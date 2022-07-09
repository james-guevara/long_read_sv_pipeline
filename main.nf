reference_fasta = file("/expanse/lustre/scratch/ux453059/temp_project/pb_pipeline/Homo_sapiens_assembly38.fasta", type: "file", checkIfExists: true)
reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)

joint_called_vcf = file("/expanse/projects/sebat1/genomicsdataanalysis/REACH_JG/concatenated_vcf/reach.sorted.norm.annotated.vcf.gz", type: "file", checkIfExists: true)
joint_called_vcf_tbi = file("${joint_called_vcf}.tbi", type: "file", checkIfExists: true)


ped_file = file("REACH.2022_01_07.psam", type: "file", checkIfExists: true)

ped = Channel
            .fromPath("REACH.2022_01_07.psam", type: "file", checkIfExists: true)
            .splitCsv(sep: "\t", header: ["family_id", "sample_id", "father_id", "mother_id", "sex", "phenotype"])
            .map { row -> tuple(row.sample_id, row.family_id) }

// The data channel will output tuples, each of which contains a sample ID, the path to the sample's input bam file, and the sample's family ID.
sample_fastq_tuple = Channel
            .fromPath("ont_fastqs_test1.tsv", type: "file", checkIfExists: true)
            .splitCsv(sep: "\t", header: ["sample_id", "reads_file"])
            .map { row -> tuple(row.sample_id, row.reads_file) }
            .join(ped)
            .map { it -> tuple(it[2], it[0], it[1]) }


process MINIMAP2 {
    input:
    tuple val(family_id), val(sample_id), path(fastq_gz)
    output:
    tuple val(family_id), val(sample_id), path("${sample_id}.sam")
    script:
    """ 
    minimap2 -t ${task.cpus} -H -x map-ont -a --MD -Y -o ${sample_id}.sam $reference_fasta $fastq_gz
    """
    stub:
    """
    touch ${sample_id}.sam
    """
}

process SORT {
    input:
    tuple val(family_id), val(sample_id), path(sam)
    output:
    tuple val(family_id), val(sample_id), path("${sam.simpleName}.bam")
    script:
    """
    samtools sort --threads ${task.cpus} --reference $reference_fasta  -o ${sample_id}.bam $sam
    """
    stub:
    """
    touch ${sample_id}.bam
    """
}

process ADD_READGROUP {
    input:
    tuple val(family_id), val(sample_id), path(bam)
    output:
    tuple val(family_id), val(sample_id), path("${sample_id}.RG.bam")
    script:
    """
    samtools addreplacerg --threads ${task.cpus} -r "@RG\tID:$sample_id\tSM:$sample_id" -o ${sample_id}.RG.bam $bam
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
    bcftools view --threads ${task.cpus} --samples ${sample_names.join(" ")} --output-type z --output-file ${family_id}.vcf.gz $joint_genotyped_vcf
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
    whatshap phase --reference $reference_fasta --ped $ped_file --indels --tag PS --merge-reads --output ${family_id}.phased.vcf.gz $family_vcf $bams 
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
    bcftools merge --threads ${task.cpus} --force-samples --output "${family_id}.phased.vcf.gz" --output-type z $vcfs
    """
    stub:
    """
    touch ${family_id}.phased.vcf.gz
    """
}

process INDEX_VCF {
    input:
    tuple val(family_id), path(phased_family_vcf)
    output:
    tuple val(family_id), path(phased_family_vcf), path("${phased_family_vcf}.tbi")
    script:
    """
    tabix $phased_family_vcf
    """
    stub:
    """
    touch ${phased_family_vcf}.tbi
    """
}

process HAPLOTAG {
    input:
    tuple val(family_id), val(sample_name), path(bam), path(bai), path(phased_family_vcf), path(phased_family_vcf_tbi)
    output:
    tuple val(sample_name), path("${bam.simpleName}.haplotag.bam"), path("${bam.simpleName}.haplotag.bam.bai")
    script:
    """
    whatshap haplotag --output-threads ${task.cpus} --reference $reference_fasta --sample $sample_name --ignore-read-groups --tag-supplementary --output-haplotag-list ${bam.simpleName}_haplotag_list.tab.gz --output ${bam.simpleName}.haplotag.bam $phased_family_vcf $bam
    samtools index -@ ${task.cpus} ${bam.simpleName}.haplotag.bam
    """
    stub:
    """
    touch ${bam.simpleName}.haplotag.bam
    """
}


process CUTESV {
    input:
    tuple val(sample_name), path(bam), path(bai)
    output:
    tuple val(sample_name), path("${bam.simpleName}.mm2.cutesv.s1.vcf")
    script:
    """
    mkdir work_folder
    cuteSV -t ${task.cpus} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.9 --report_readid --min_support 1 --genotype $bam $reference_fasta ${bam.simpleName}.mm2.cutesv.s1.vcf work_folder 
    """
    stub:
    """
    touch ${bam.simpleName}.mm2.cutesv.s1.vcf
    """
}

process POSTPROCESS_CUTESV {
    input:
    tuple val(sample_name), path(vcf)
    output:
    tuple val(sample_name), path("${vcf.simpleName}.mm2.cutesv.s1.vcf.gz")
    script:
    """
    bcftools view --threads ${task.cpus} --include 'POS>0' $vcf | bcftools sort --output-type z --output ${vcf.simpleName}.mm2.cutesv.s1.fixed.vcf.gz
    """
    stub:
    """
    touch ${vcf.simpleName}.mm2.cutesv.s1.vcf.gz
    """
}

process SNIFFLES {
    input:
    tuple val(sample_name), path(bam) 
    output:
    tuple val(sample_name), path("${bam.simpleName}.mm2.sniffles.s1.vcf")
    script:
    """
    sniffles -t ${task.cpus} --num_reads_report -1 --tmp_file tmp1 --min_support 1 --cluster -m $map_sort_bam_sniffles -v ${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf
    """
    stub:
    """
    touch ${bam.simpleName}.mm2.sniffles.s1.vcf
    """
}

process POSTPROCESS_SNIFFLES {
    input:
    tuple val(sample_name), path(vcf)
    output:
    tuple val(sample_name), path("${vcf.simpleName}.mm2.sniffles.s1.vcf.gz")
    script:
    """
    bash postprocess_sniffles.sh $vcf
    """
    stub:
    """
    touch ${vcf.simpleName}.mm2.sniffles.s1.vcf.gz
    """
}


process SVIM {
    input:
    tuple val(sample_name), path(bam)
    output:
    tuple val(sample_name), path("${bam.simpleName}.mm2.svim.vcf")
    script:
    """
    svim alignment --read_names --zmws out_svim/  $map_sort_bam_svim $reference_fasta
    mv out_svim/variants.vcf ${bam.simpleName}.mm2.svim.vcf
    """
    stub:
    """
    touch ${bam.simpleName}.mm2.svim.vcf
    """
}

process POSTPROCESS_SVIM {
    input:
    tuple val(sample_name), path(vcf)
    output:
    tuple val(sample_name), path("${vcf.simpleName}.mm2.svim.s1.vcf.gz")
    script:
    """
    bash postprocess_svim.sh $vcf
    """
    stub:
    """
    touch ${vcf.simpleName}.mm2.svim.s1.vcf.gz
    """
}

process JASMINE {
    input:
    tuple val(sample_name), path(bam), path(cutesv_vcf), path(sniffles_vcf), path(svim_vcf)
    output:
    tuple val(sample_name), path("${sample_name}_jasmine_iris.vcf"), path(cutesv_vcf), path(sniffles_vcf), path(svim_vcf)
    script:
    """
    echo $cutesv_vcf    > vcf_paths.txt
    echo $sniffles_vcf >> vcf_paths.txt
    echo $svim_vcf     >> vcf_paths.txt
    VCF_FILE_LIST=vcf_paths.txt

    echo $bam > bam_path.txt
    BAM_FILE_LIST=bam_path.txt

    jasmine file_list=\$VCF_FILE_LIST out_file=${sample_name}_jasmine_iris.vcf genome_file=$reference_fasta bam_list=\$BAM_FILE_LIST out_dir=. --output_genotypes --dup_to_ins --normalize_type --ignore_strand threads=${task.cpus} --run_iris iris_args=--keep_long_variants,threads=${task.cpus}
    """
    stub:
    """
    touch ${sample_name}_jasmine_iris.vcf
    """
}

process POSTPROCESS_JASMINE {
    input:
    tuple val(sample_name), path(jasmine_vcf), path(cutesv_vcf), path(sniffles_vcf), path(svim_vcf)
    output:
    tuple val(sample_name), path("${jasmine_vcf.simpleName}.s1.vcf.gz"), path("${jasmine_vcf.simpleName}.s1.vcf.gz.tbi")
    script:
    """
    bash postprocess_jasmine.sh $jasmine_vcf $cutesv_vcf $sniffles_vcf $svim_vcf
    """
    stub:
    """
    touch ${jasmine_vcf.simpleName}.s1.vcf.gz
    touch ${jasmine_vcf.simpleName}.s1.vcf.gz.tbi
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

    // Create the family tuples that will be used for phasing
    family_tuple = sample_bam_tuple.groupTuple()
    family_tuple
                .branch {
                    large_families: it[1].size() >= 5
                    standard_family: true
                }
                .set { families_to_be_phased }            
    // Deal with large families by splitting them up into trios (and then merging them afterward)
    families_to_be_phased.large_families
                                        .transpose()
                                        .branch {
                                            father: pedigree_dictionary[it[1]] == "Father"
                                            mother: pedigree_dictionary[it[1]] == "Mother"
                                            child: true
                                        }
                                        .set { large_family_samples }
    large_family_trios = large_family_samples.father.join(large_family_samples.mother).combine(large_family_samples.child, by: 0).map { it -> tuple( it[0], [ it[1], it[4], it[7] ], [ it[2], it[5], it[8] ], [ it[3], it[6], it[9] ] ) }

    // Phasing (subset the joint-genotyped SNP VCF beforehand)
    subsetted_family_tuples = SUBSET_VCF(large_family_trios.mix(families_to_be_phased.standard_family), joint_called_vcf)
    phased_family_tuples = PHASE(subsetted_family_tuples)
    phased_family_tuples
                        .groupTuple()
                        .branch {
                                families_to_merge: it[1].size() > 1 
                                standard_families: true
                        }
                        .set { families_to_merge_and_standard_families }
    // Merging the large family trio VCFs (and then use the mix operator to put all the family VCFs into one channel, all_phased_families_tuple) 
    merged_families = MERGE_FAMILY_VCFS(families_to_merge_and_standard_families.families_to_merge)
    all_phased_families_tuple_with_index = INDEX_VCF(families_to_merge_and_standard_families.standard_families.transpose().mix(merged_families))


    // // Haplotag each sample's bam file using its phased family VCF
    haplotag_bam_tuple = HAPLOTAG(sample_bam_tuple.combine(all_phased_families_tuple_with_index, by: 0))

    // // Run variant callers
    cutesv_tuple = POSTPROCESS_CUTESV(CUTESV(haplotag_bam_tuple))
    // sniffles_tuple = POSTPROCESS_SNIFFLES(SNIFFLES(haplotag_bam_tuple))
    // svim_tuple = POSTPROCESS_SVIM(SVIM(haplotag_bam_tuple))

    // // Merge variant callers (within each sample) using Jasmine
    // jasmine_tuple = POSTPROCESS_JASMINE(JASMINE(haplotag_bam_tuple.join(cutesv_tuple).join(sniffles_tuple).join(svim_tuple)))
}
