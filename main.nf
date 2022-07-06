reference_fasta = file("/expanse/lustre/scratch/ux453059/temp_project/pb_pipeline/Homo_sapiens_assembly38.fasta", type: "file", checkIfExists: true)
reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)

joint_called_vcf = file("/expanse/projects/sebat1/genomicsdataanalysis/REACH_JG/concatenated_vcf/reach.sorted.norm.annotated.vcf.gz", type: "file", checkIfExists: true)
joint_called_vcf_tbi = file("${joint_called_vcf}.tbi", type: "file", checkIfExists: true)

ped = Channel
            .fromPath("REACH.2022_01_07.psam", type: "file", checkIfExists: true)
            .splitCsv(sep: "\t", header: ["family_id", "sample_id", "father_id", "mother_id", "sex", "phenotype"])
            .map { row -> tuple(row.sample_id, row.family_id) }

// The data channel will output tuples, each of which contains a sample ID, the path to the sample's input bam file, and the sample's family ID.
sample_fastq_tuple = Channel
            .fromPath("ont_fastqs.tsv", type: "file", checkIfExists: true)
            .splitCsv(sep: "\t", header: ["sample_id", "reads_file"])
            .map { row -> tuple(row.sample_id, row.reads_file) }
            .join(ped)
            .map { it -> tuple(it[2], it[0], it[1]) }


process MINIMAP2 {
    input:
    tuple val(family_id), val(sample_id), path(fastq_gz)
    output:
    tuple val(family_id), val(sample_id), path("${fastq_gz.simpleName}.sam")
    script:
    """ 
    minimap2 -t 4 -H -x map-hifi -a --MD -Y -o ${fastq_gz.simpleName}.sam $reference_fasta $fastq_gz
    """
    stub:
    """
    touch ${fastq_gz.simpleName}.sam
    """
}

process SORT {
    input:
    tuple val(family_id), val(sample_id), path(sam)
    output:
    tuple val(family_id), val(sample_id), path("${sam.simpleName}.bam")
    script:
    """
    samtools sort --reference $reference_fasta -@ ${task.cpus} -o ${sam.simpleName}.bam $sam
    """
    stub:
    """
    touch ${sam.simpleName}.bam
    """
}

process ADD_READGROUP {
    input:
    tuple val(family_id), val(sample_id), path(bam)
    output:
    tuple val(family_id), val(sample_id), path("${sample_id}.bam")
    script:
    """
    samtools addreplacerg -r "@RG\tID:$sample_id\tSM:$sample_id" -o ${sample_id}.RG.bam $bam
    """
    stub:
    """
    touch ${sample_id}.RG.bam
    """
} 

process INDEX {
    input:
    tuple val(family_id), val(sample_id), path(bam)
    output:
    tuple val(family_id), val(sample_id), path(bam), path("${bam}.bai")
    script:
    """
    samtools index -@ ${task.cpus} $bam
    """
    stub:
    """
    touch ${bam}.bai
    """
}

process SUBSET_VCF {
    input:
    tuple val(family_id), val(sample_names), path(bams), path(bais)
    path(joint_genotyped_vcf)
    output:
    tuple val(family_id), val(sample_names), path(bams), path(bais), path("${family_id}.vcf.gz")
    script:
    """
    bcftools --threads ${task.cpus} --samples ${sample_names.join(" ")} --output-type z --output ${family_id}.vcf.gz $joint_genotyped_vcf
    """
    stub:
    """
    touch "${family_id}.vcf.gz"
    """
}

process PHASE {
    input:
    tuple val(family_id), val(sample_names), path(bams), path(bais), path(family_vcf)
    output:
    tuple val(family_id), path("${family_id}.phased.vcf.gz")
    script:
    """
    whatshap phase --reference $reference_fasta --ped $ped --indels --tag PS --merge-reads --ignore-read-groups --output ${family_id}.phased.vcf.gz $family_vcf $bams 
    """
    stub:
    """
    touch ${family_id}.phased.vcf.gz
    """
}

process MERGE_FAMILY_VCFS {
    input:
    tuple val(family_id), path(vcfs, stageAs: "?.vcf.gz") 
    output:
    tuple val(family_id), path("${family_id}.phased.vcf.gz")
    script:
    """
    bcftools merge --force-samples --output "${family_id}.phased.vcf.gz" --output-type z $vcfs
    """
    stub:
    """
    echo $vcfs
    touch ${family_id}.phased.vcf.gz
    """
}



process HAPLOTAG {
    input:
    tuple val(family_id), val(sample_name), path(bam), path(bai), path(phased_family_vcf)
    output:
    tuple val(sample_name), path("${bam.simpleName}.haplotag.bam")
    script:
    """
    whatshap haplotag --reference $reference_fasta --sample $sample_name --ignore-read-groups --tag-supplementary --output-haplotag-list ${bam.simpleName}_haplotag_list.tab.gz --output ${bam.simpleName}.haplotag.bam $phased_family_vcf $bam
    """
    stub:
    """
    touch ${bam.simpleName}.haplotag.bam
    """
}


def make_pedigree_dictionary (pedigree_filepath) {
    pedigree_table = new File(pedigree_filepath).readLines().collect{ it.tokenize("\t") }
    pedigree_dictionary = [:]
    pedigree_table.each{ item ->
        if (item[2] == "0" && item[3] == "0" && item[4] == "1") { 
            pedigree_dictionary[item[1]] = "Father"
        }
        else if (item[2] == "0" && item[3] == "0" && item[4] == "2") {
            pedigree_dictionary[item[1]] = "Mother"
        }
        else {
            pedigree_dictionary[item[1]] = "Child"
        }
    }
    return pedigree_dictionary 
}

workflow {
    pedigree_dictionary = make_pedigree_dictionary("REACH.2022_01_07.psam")

    // Mapping steps (uses minimap2 for alignment to create a sam file, which is then converted to a bam file, and then index and sort the output bam file)
    sample_bam_tuple = INDEX(ADD_READGROUP(SORT(MINIMAP2(sample_fastq_tuple))))

    family_tuple = sample_bam_tuple.groupTuple()
    family_tuple
                .branch {
                    large_families: it[1].size() >= 5
                    standard_family: true
                }
                .set { families_to_be_phased }            

    families_to_be_phased.large_families
                                        .transpose()
                                        .branch {
                                            father: pedigree_dictionary[it[1]] == "Father"
                                            mother: pedigree_dictionary[it[1]] == "Mother"
                                            child: true
                                        }
                                        .set { large_family_samples }

    large_family_trios = large_family_samples.father.join(large_family_samples.mother).combine(large_family_samples.child, by: 0).map { it -> tuple( it[0], [ it[1], it[4], it[7] ], [ it[2], it[5], it[8] ], [ it[3], it[6], it[9] ] ) }

    // Now we can run the phasing on all types of family structure
    subsetted_family_tuples = SUBSET_VCF(large_family_trios.mix(families_to_be_phased.standard_family), joint_called_vcf)
    phased_family_tuples = PHASE(subsetted_family_tuples)
    phased_family_tuples
                        .groupTuple()
                        .branch {
                                families_to_merge: it[1].size() > 1 
                                standard_families: true
                        }
                        .set { families_to_merge_and_standard_families }
    merged_families = MERGE_FAMILY_VCFS(families_to_merge_and_standard_families.families_to_merge)
    all_phased_families_tuple = families_to_merge_and_standard_families.standard_families.transpose().mix(merged_families)
    // Haplotag each sample's bam file using its phased family VCF
    haplotag_bam_tuple = HAPLOTAG(sample_bam_tuple.combine(all_phased_families_tuple, by: 0))
    haplotag_bam_tuple.view()

    // // Five member families
    // families_to_be_phased.five_member_family
    //                                         .transpose()
    //                                         .branch {
    //                                             father: pedigree_dictionary[it[1]] == "Father"
    //                                             mother: pedigree_dictionary[it[1]] == "Mother"
    //                                             child: true
    //                                         }
    //                                         .set { five_member_family_samples }
    // five_member_family_trios = five_member_family_samples.father.join(five_member_family_samples.mother).combine(five_member_family_samples.child, by: 0).map { it -> tuple( it[0], [ it[1], it[4], it[7] ], [ it[2], it[5], it[8] ], [ it[3], it[6], it[9] ] ) }


    // // Now we can run the phasing on all types of family structure
    // subsetted_family_tuples = SUBSET_VCF(five_member_family_trios.mix(families_to_be_phased.standard_family), joint_called_vcf)
    // phased_family_tuples = PHASE(subsetted_family_tuples)
    // phased_family_tuples
    //                     .groupTuple()
    //                     .branch {
    //                             families_to_merge: it[1].size() > 1 
    //                             standard_families: true
    //                     }
    //                     .set { families_to_merge_and_standard_families }
    // merged_families = MERGE_FAMILY_VCFS(families_to_merge_and_standard_families.families_to_merge)

    // all_phased_families_tuple = families_to_merge_and_standard_families.standard_families.transpose().mix(merged_families)

    // // Haplotag each sample's bam file using its phased family VCF
    // haplotag_bam_tuple = HAPLOTAG(sample_bam_tuple.combine(all_phased_families_tuple, by: 0))
}
