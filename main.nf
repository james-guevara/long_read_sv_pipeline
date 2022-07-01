nextflow.enable.dsl = 2

reference_fasta = file("test_data/Homo_sapiens_assembly38.fasta", type: "file", checkIfExists: true)
reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)

joint_called_vcf = file("/expanse/projects/sebat1/genomicsdataanalysis/REACH_JG/concatenated_vcf/reach.sorted.norm.annotated.vcf.gz", type: "file", checkIfExists: true)
joint_called_vcf_tbi = file("${joint_called_vcf}.tbi", type: "file", checkIfExists: true)

ped = file("REACH.2022_01_07.20.psam", type: "file", checkIfExists: true)

process BAM_TO_FASTQ {
    input:
    tuple val(sample_id), path(bam), val(family_id)
    output:
    tuple val(family_id), val(sample_id), path("${bam.simpleName}.fastq")
    script:
    """
    samtools bam2fq -@ ${task.cpus} $bam > ${bam.simpleName}.fastq
    """
    stub:
    """
    touch ${bam.simpleName}.fastq 
    """
}

process BGZIP {
    input:
    tuple val(family_id), val(sample_id), path(fastq)
    output:
    tuple val(family_id), val(sample_id), path("${fastq}.gz")
    script:
    """
    bgzip $fastq
    """
    stub:
    """
    touch ${fastq}.gz
    """
}

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
    tuple val(family_id), val(sample_id), path("${bam.simpleName}.${sample_id}.bam")
    script:
    """
    samtools addreplacerg -r "@RG\tID:$sample_id\tSM:$sample_id" -o ${bam.simpleName}.${sample_id}.bam $bam
    """
    stub:
    """
    touch ${bam.simpleName}.${sample_id}.bam
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
    tuple val(family_id), path("${family_vcf.simpleName}.phased.vcf.gz")
    script:
    """
    whatshap phase --reference $reference_fasta --ped $ped --indels --tag PS --merge-reads --ignore-read-groups --output ${family_vcf.simpleName}.phased.vcf.gz $family_vcf $bams 
    """
    stub:
    """
    touch ${family_vcf.simpleName}.phased.vcf.gz
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

process CUTESV {
    input:
    tuple val(sample_name), path(bam)
    output:
    tuple val(sample_name), path("${bam.simpleName}.cutesv.vcf")
    script:
    """
    mkdir work_folder 
    cuteSV --threads ${task.cpu} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.9 --report_readid --min_support 1 --genotype $bam $reference_fasta ${bam.simpleName}.cutesv.vcf work_folder 
    """
    stub:
    """
    touch ${bam.simpleName}.cutesv.vcf
    """
}

process FIX_CUTESV {
    input:
    tuple val(sample_name), path(vcf)
    output:
    tuple val(sample_name), path("${vcf.simpleName}.cutesv.fixed.sorted.vcf")
    script:
    """
    bcftools view --include 'POS>0' $vcf | bcftools sort --output ${vcf.simpleName}.cutesv.fixed.sorted.vcf
    """    
    stub:
    """
    touch ${vcf.simpleName}.cutesv.fixed.sorted.vcf
    """
}

process SNIFFLES {
    input:
    tuple val(sample_name), path(bam)
    output:
    tuple val(sample_name), file("${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf") into sniffles_output_raw
    script:
    """
    sniffles --ccs_reads -t 32 --num_reads_report -1 --tmp_file tmp1 --min_support 1 --cluster -m $map_sort_bam_sniffles -v ${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf
    """
    stub:
    """
    touch ${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf
    """
}

process FIX_SNIFFLES {
    input:
    tuple val(sample_name), path(vcf)
    output:
    tuple val(sampleName),file("${sniffles_output_raw.simpleName}.mm2.sniffles.s1.fixed.vcf") into sniffles_output
    shell:
    '''
    if grep -Fq STRANDBIAS !{sniffles_output_raw}; then
        awk 'BEGIN{FS="\t";OFS="\t"} \
        { \
            if ($1=="#CHROM") { \
                print("##FILTER=<ID=STRANDBIAS,Description=\\"Strand Bias\\">"); \
            } \
            print $0; \
        }' !{sniffles_output_raw} | \
        bcftools reheader -s <(echo !{sampleName}) | \
        bcftools sort -o !{sniffles_output_raw.simpleName}.mm2.sniffles.s1.fixed.vcf 
    else
        bcftools reheader -s <(echo !{sampleName}) !{sniffles_output_raw} | \
        bcftools sort -o !{sniffles_output_raw.simpleName}.mm2.sniffles.s1.fixed.vcf

    fi
    '''
    stub:
    """
    touch ${sniffles_output_raw.simpleName}.mm2.sniffles.s1.fixed.vcf
    """
}

process SVIM {
    input:
    tuple val(sample_name), path(bam)
    output:
    tuple val(sample_name), path("${bam.simpleName}.svim.vcf")
    script:
    """
    svim alignment --read_names --zmws out/ $bam $reference_fasta
    mv out_svim/variants.vcf ${bam.simpleName}.svim.vcf
    """
    stub:
    """
    touch ${bam.simpleName}.svim.vcf
    """    
}

process FIX_SVIM {
    input:
    tuple val(sample_name), path(vcf)
    output:
    tuple val(sample_name), path("${vcf.simpleName}.svim.fixed.vcf")
    shell:
    '''
    awk 'BEGIN{FS="\t";OFS="\t"} \
        { \
            if ($1=="#CHROM") { \
                print("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\\"# high-quality variant reads\\">"); \
            } \
            gsub("DUP:INT","DUP_INT",$0); \
            gsub("DUP:TANDEM","DUP",$0); \
            gsub("READS","RNAMES",$0); \
            DV="."; \
            SVTYPE_INV = "F"; \
            if ($8~"SVTYPE=INV;") { \
                SVTYPE_INV = "T"; \
            } \
            n=split($8,p_info,";"); \
            for (i=1;i<=n;i++) { \
                split(p_info[i],pp,"="); \
                if ((SVTYPE_INV=="T") && (pp[1]=="END")) { \
                    svlen = pp[2] - $2 + 1; \
                    $8=$8";SVLEN="svlen; \
                } \
                if (pp[1]=="SUPPORT") { \
                    DV = pp[2]; \
                    $9=$9":DV"; \
                    $10=$10":"DV; \
                } \
            } \
            print $0; \
        }' !{svim_output_raw} | \
     bcftools view -i 'SVTYPE=="DEL" || SVTYPE=="INS" || SVTYPE=="DUP" || SVTYPE=="INV" || SVTYPE=="BND"' | \
     bcftools view -i 'SUPPORT>=1' | bcftools sort -o !{svim_output_raw.simpleName}.mm2.svim.s1.fixed.vcf
    '''
    stub:
    """
    touch ${vcf.simpleName}.svim.fixed.vcf
    """
}

process PBSV {
    input:
    tuple val(sample_name), path(bam)
    output:
    tuple val(sample_name), path("${bam.simpleName}.pbsv.vcf")
    script:
    """
    pbsv discover --sample $sample_name --tandem-repeats $trBedFile $bam ${bam.simpleName}.svsig.gz
    pbsv call --ccs -O 1 $reference_fasta ${bam.simpleName}.svsig.gz ${bam.simpleName}.pbsv.vcf
    """
    stub:
    """
    touch ${bam.simpleName}.pbsv.vcf
    """
}

// process FIX_PBSV {
//     input:
//     tuple val(sample_name), path(vcf)
//     output:
//     tuple val(sample_name), path("${vcf.simpleName}.pbsv.fixed.vcf")
// 
//     shell:
//     '''
//     header1='##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">'
//     header2='##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">'
// 
//     cat <(bcftools view -h !{pbsv_output_raw} | head -n -1) \
//     <(echo -e $header1) <(echo -e $header2) \
//     <(bcftools view -h !{pbsv_output_raw} | tail -n 1) \
//     <(bcftools view -H !{pbsv_output_raw} | \
//     awk 'BEGIN{FS="\t";OFS="\t"} \
//     { \
//         split($NF,p,":"); \
//         AD = p[2]; \
//         split(AD, p_ad, ","); \
//         DR = p_ad[1]; \
//         DV = p_ad[2]; \
//         $(NF-1) = $(NF-1)":DR:DV"; \
//         $NF = $NF":"DR":"DV; \
//         print $0; \
//     }') | bcftools sort -o !{pbsv_output_raw.simpleName}.mm2.pbsv.s1.fixed.vcf
//     '''
//     stub:
//     """
//     touch ${pbsv_output_raw.simpleName}.mm2.pbsv.s1.fixed.vcf
//     """
// }

process MERGE {
    input:
    tuple val(sampleName),file(map_sort_bam),file(map_sort_bam_index),file(vcfs_sn),file(vcfs_cu),file(vcfs_pb),file(vcfs_sv) from map_sort_bam.join(sniffles_output).join(cutesv_output).join(pbsv_output).join(svim_output)         
 
    output:
    tuple val(sampleName),file("Jasmine_Iris_merge_scps_haptag_${sampleName}.vcf") into jasmine_output_raw
 
    shell:
    '''
    OUT_FILE=Jasmine_Iris_merge_scps_haptag_!{sampleName}.vcf
    echo "!{map_sort_bam.name}\n!{map_sort_bam.name}\n!{map_sort_bam.name}\n!{map_sort_bam.name}" > bam_files_scps.txt
    echo "!{vcfs_sn.name}\n!{vcfs_cu.name}\n!{vcfs_pb.name}\n!{vcfs_sv.name}" > vcf_files_scps.txt
    BAMS_LIST=bam_files_scps.txt
    VCFS_FILE_LIST=vcf_files_scps.txt
    TMP_DIR=/scratch/$USER/job_$SLURM_JOBID

    jasmine file_list=$VCFS_FILE_LIST out_file=$OUT_FILE genome_file=!{fasta} bam_list=$BAMS_LIST out_dir=$TMP_DIR --output_genotypes \
      --dup_to_ins --normalize_type --ignore_strand threads=120 --run_iris iris_args=--keep_long_variants,threads=124
    '''

    stub:
    """
    echo $vcfs_cu.baseName > Z2.txt
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.vcf
    """ 
}


workflow {
    ped = channel
            .fromPath("REACH.2022_01_07.20.psam", type: "file", checkIfExists: true)
            .splitCsv(sep: "\t", header: ["family_id", "sample_id", "father_id", "mother_id", "sex", "phenotype"])
            .map { row -> tuple(row.sample_id, row.family_id) }

    // The data channel will output tuples, each of which contains a sample ID, the path to the sample's input bam file, and the sample's family ID.
    data = channel
            .fromPath("test_data/samples_files.tsv", type: "file", checkIfExists: true)
            .splitCsv(sep: "\t", header: ["sample_id", "reads_file"])
            .map { row -> tuple(row.sample_id, row.reads_file) }
            .join(ped)

    // Preprocessing steps (converts the unaligned reads bam file to a fastq.gz file)
    sample_fastq_tuple = BGZIP(BAM_TO_FASTQ(data))
    
    // Mapping steps (uses minimap2 for alignment to create a sam file, which is then converted to a bam file, and then index and sort the output bam file)
    sample_bam_tuple = INDEX(ADD_READGROUP(SORT(MINIMAP2(sample_fastq_tuple))))

    // Subset a family VCF and then phase it using the bam files in the family
    family_tuple = SUBSET_VCF(sample_bam_tuple.groupTuple(), joint_called_vcf)
    phased_family_tuple = PHASE(family_tuple) 

    // Haplotag each sample's bam file using its phased family VCF
    haplotag_bam_tuple = HAPLOTAG(sample_bam_tuple.join(phased_family_tuple))
}
