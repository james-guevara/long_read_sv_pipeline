params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    Long-read structural variation pipeline (mapping, variant calling, genotyping)
    ============================
    Required:
    --sample_fastq_table_tsv    Tab-delimited file, where the first column corresponds to the sample name and the second column corresponds to the file path of the fastq reads input.
    --reference_fasta           Reference genome to which the reads are aligned.
    --pedigree_file             Pedigree file (.fam, .ped, or .psam) containing all the samples in the cohort.
    --sequencing_mode           The sequencing technology used to generate the input reads. Must be one of: ccs_hifi, ccs_subread_fallback, or ont
    --tandem_repeats_bed        File path of tandem repeats bed file (used for Sniffles and ...).

    One of these modes are required for phasing:
    --cohort_vcf                File path of the cohort VCF (that contains all the samples in the cohort).
    --family_vcfs               Tab-delimited file, where the first column corresponds to the family names and the second column corresponds to the family VCF filepath.
    """.stripIndent()
    exit 0
}

params.reference_fasta = file("resources/Homo_sapiens_assembly38.fasta", type: "file", checkIfExists: true)
params.reference_fasta_fai = file("${params.reference_fasta}.fai", type: "file", checkIfExists: true)
params.sample_fastq_table_tsv = file("testing_data/sample_files_Ashkenazi_tests.tsv", type: "file", checkIfExists: true)
params.pedigree_file = file("testing_data/AshkenazimTrio.ped", type: "file", checkIfExists: true)
params.sequencing_mode = ""
sequencing_modes = ["ccs_hifi", "ccs_subread_fallback", "ont"]
if (sequencing_modes.contains(params.sequencing_mode) == false) {
    log.info "Sequencing mode parameter must be one of: ccs_hifi, ccs_subread_fallback, or ont"
    exit(0)
}
params.tandem_repeats_bed =  file("resources/human_GRCh38_no_alt_analysis_set.trf.bed", type: "file", checkIfExists: true)
params.TR_file_TRF_target = file("/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Simple_Repeats_TRF_annot.bed", type: "file", checkIfExists: true)
params.TR_file_TRF_all =  file("/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Simple_Repeats_TRF_annot.bed", type: "file", checkIfExists: true)
params.TR_file_RM_Simpe =  file("/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Repeats_Masker_Simple_repeat.bed", type: "file", checkIfExists: true)
params.TR_file_MG =  file("/home/smmortazavi/HUMAN_DATA/REF/REPEATS/hg38_ver13.bed", type: "file", checkIfExists: true)
// Choose between these 2 types of files
params.cohort_vcf = ""
params.family_vcfs = ""

chromIndices = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]

log.info """
Running with parameters:
--sample_fastq_table_tsv : ${params.sample_fastq_table_tsv}
--reference_fasta : ${params.reference_fasta}
--pedigree_file : ${params.pedigree_file}
--sequencing_mode : ${params.sequencing_mode}
--tandem_repeats_bed : ${params.tandem_repeats_bed}
--cohort_vcf : ${params.cohort_vcf}
--family_vcfs : ${params.family_vcfs}
""".stripIndent()

process MAP {
    input: 
    tuple val(sample_name), path(fastq_gz)
    output: 
    tuple val(sample_name), path("${sample_name}.bam")
 
    script:
    if (params.sequencing_mode == "ccs_hifi") 
    """ 
    minimap2 -t ${task.cpus} -a -x map-hifi --MD -Y -R '@RG\\tID:${sample_name}\\tSM:${sample_name}'  ${params.reference_fasta} $fastq_gz | samtools sort --reference ${params.reference_fasta} --threads ${task.cpus} --output-fmt BAM -o ${sample_name}.bam
    """
    else if (params.sequencing_mode == "ccs_subread_fallback") // I might have to change the -Q parameter (though I don't see it in the most recent version of minimap2)
    """
    minimap2 -t ${task.cpus} -a -x map-pb --MD -Y -R '@RG\\tID:${sample_name}\\tSM:${sample_name}' ${params.reference_fasta} $fastq_gz | samtools sort --reference ${params.reference_fasta} --threads ${task.cpus} --output-fmt BAM -o ${sample_name}.bam
    """
    else if (params.sequencing_mode == "ont")
    """
    minimap2 -t ${task.cpus} -a -x map-ont  --MD -Y -R '@RG\\tID:${sample_name}\\tSM:${sample_name}'  ${params.reference_fasta} $fastq_gz | samtools sort --reference ${params.reference_fasta} --threads ${task.cpus} --output-fmt BAM -o ${sample_name}.bam
    """
 
    stub: 
    """
    touch ${sample_name}.bam 
    """ 
}

process INDEX_BAM {
    input:
    tuple val(sample_name), path(bam)
    output:
    tuple val(sample_name), path(bam), path("${bam}.bai")

    script:
    """
    samtools index -@ ${task.cpus} $bam
    """
    
    stub:
    """
    touch ${bam}.bai
    """
}

process COVERAGE {
    publishDir "results/$sample_name", mode: "copy"

    input:
    tuple val(sample_name), path(bam), path(bai)
    output:
    tuple path("${sample_name}.mosdepth.global.dist.txt"), path("${sample_name}.mosdepth.summary.txt")

    script:
    """
    mosdepth --threads ${task.cpus} --no-per-base --fasta ${params.reference_fasta} --fast-mode $sample_name $bam
    """

    stub:
    """
    touch ${sample_name}.mosdepth.global.dist.txt
    touch ${sample_name}.mosdepth.summary.txt
    """
}

process SUBSET_VCF {
    input:
    tuple val(family_name), val(sample_names), path(bams), path(bais)
    path(cohort_vcf)
    output:
    tuple val(family_name), val(sample_names), path(bams), path(bais), path("${family_name}.vcf.gz")

    script:
    """
    bcftools view --threads ${task.cpus} --samples ${sample_names.join(",")} --output-type z --output-file ${family_name}.vcf.gz $cohort_vcf
    """

    stub:
    """
    touch ${family_name}.vcf.gz
    """
}

process INDEX_VCF {
    input:
    tuple val(family_name), val(sample_names), path(bams), path(bais), path(family_vcf)
    output:
    tuple val(family_name), val(sample_names), path(bams), path(bais), path(family_vcf), path("${family_vcf}.tbi")

    script:
    """
    tabix $family_vcf 
    """

    stub:
    """
    touch ${family_vcf}.tbi
    """
}

process PHASE {
    input:
    tuple val(family_name), val(sample_names), path(bams), path(bais), path(family_vcf), path(family_vcf_tbi)
    output:
    tuple val(family_name), path("${family_name}.phased.vcf.gz"), path("${family_name}.phased.vcf.gz.tbi")

    script:
    """
    whatshap phase --reference ${params.reference_fasta} --ped ${params.pedigree_file} --indels --tag PS --output ${family_name}.phased.vcf.gz $family_vcf $bams
    tabix ${family_name}.phased.vcf.gz
    """

    stub:
    """
    touch ${family_name}.phased.vcf.gz
    touch ${family_name}.phased.vcf.gz.tbi
    """
}

process MERGE_FAMILY_VCFS {
    input:
    tuple val(family_name), path(vcfs, stageAs: "?.vcf.gz"), path(vcf_tbis, stageAs: "?.vcf.gz.tbi")
    output:
    tuple val(family_name), path("${family_name}.phased.vcf.gz"), path("${family_name}.phased.vcf.gz.tbi")

    script:
    """
    bcftools merge --threads ${task.cpus} --force-samples --output ${family_name}.phased.vcf.gz --output-type z $vcfs
    tabix ${family_name}.phased.vcf.gz
    """

    stub:
    """
    touch ${family_name}.phased.vcf.gz
    touch ${family_name}.phased.vcf.gz.tbi
    """
}

process HAPLOTAG {
    input:
    tuple val(family_name), val(sample_name), path(bam), path(bai), path(phased_family_vcf), path(phased_family_vcf_tbi)
    output:
    tuple val(sample_name), path("${sample_name}.haplotag.bam")

    script:
    """
    whatshap haplotag --reference ${params.reference_fasta} --sample $sample_name --skip-missing-contigs --ignore-read-groups --tag-supplementary --output-haplotag-list ${sample_name}_haplotag_list.tab.gz --output ${sample_name}.haplotag.bam $phased_family_vcf $bam
    """

    stub:
    """
    touch ${sample_name}.haplotag.bam
    """
}

process INDEX_BAM2 {
    publishDir "results/$sample_name", mode: "copy"

    input:
    tuple val(sample_name), path(bam)
    output:
    tuple val(sample_name), path(bam), path("${bam}.bai")

    script:
    """
    samtools index -@ ${task.cpus} $bam
    """

    stub:
    """
    touch ${bam}.bai
    """
}

process CUTESV {
    publishDir "results/$sample_name", mode: "copy"

    input:
    tuple val(sample_name), path(bam), path(bai)
    output:
    tuple val(sample_name), path("${sample_name}.cutesv.vcf")

    script:
    if (params.sequencing_mode == "ccs_hifi" || params.sequencing_mode == "ccs_subread_fallback")
    """
    mkdir tmp
    cuteSV -t ${task.cpus} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.9 --report_readid --min_support 1 --genotype $bam ${params.reference_fasta} ${sample_name}.cutesv.vcf tmp
    """
    else if (params.sequencing_mode == "ont")
    """
    mkdir tmp
    cuteSV -t ${task.cpus} --max_cluster_bias_INS 100  --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100  --diff_ratio_merging_DEL 0.3 --report_readid --min_support 1 --genotype $bam ${params.reference_fasta} ${sample_name}.cutesv.vcf tmp
    """
    
    stub:
    """
    touch ${sample_name}.cutesv.vcf
    """
}

process POSTPROCESS_CUTESV {
    input:
    tuple val(sample_name), path(vcf, stageAs: "input.vcf")
    output:
    tuple val(sample_name), path("${sample_name}.cutesv.vcf")

    script:
    """
    mkdir tmp
    bcftools view --include 'POS>0' $vcf | bcftools sort --temp-dir tmp --output-file ${sample_name}.cutesv.vcf
    """

    stub:
    """
    touch ${sample_name}.cutesv.vcf
    """
}

process SNIFFLES {
    publishDir "results/$sample_name", mode: "copy"

    input:
    tuple val(sample_name), path(bam), path(bai)
    output:
    tuple val(sample_name), path("${sample_name}.sniffles.vcf")

    script:
    """
    sniffles --threads ${task.cpus} --minsupport 1 --output-rnames --reference ${params.reference_fasta} --tandem-repeats ${params.tandem_repeats_bed} --input $bam --vcf ${sample_name}.sniffles.vcf
    """

    stub:
    """
    touch ${sample_name}.sniffles.vcf
    """
}

process POSTPROCESS_SNIFFLES {
    input:
    tuple val(sample_name), path(vcf, stageAs: "input.vcf")
    output:
    tuple val(sample_name), path("${sample_name}.sniffles.vcf")

    shell:
    '''
    mkdir tmp
    if grep -Fq STRANDBIAS !{vcf}; then
        awk '
            BEGIN {
                FS = "\t"
                OFS = "\t"
            }
            {
                if ($1 == "#CHROM") {
                    print("##FILTER=<ID=STRANDBIAS,Description=\\"Strand Bias\\">");
                }
                print $0;
            }
            END {
            }' !{vcf} | bcftools reheader -s <(echo !{sample_name}) | bcftools sort --temp-dir tmp --output !{sample_name}.sniffles.vcf 
    else
        bcftools reheader --samples <(echo !{sample_name}) !{vcf} | bcftools sort --temp-dir tmp --output !{sample_name}.sniffles.vcf
    fi
    '''

    stub:
    """
    touch ${sample_name}.sniffles.vcf
    """
}

process SVIM {
    publishDir "results/$sample_name", mode: "copy"

    input:
    tuple val(sample_name), path(bam), path(bai)
    output:
    tuple val(sample_name), path("${sample_name}.svim.vcf")

    script:
    """
    svim alignment --read_names --zmws out_svim/ $bam ${params.reference_fasta}
    mv out_svim/variants.vcf ${sample_name}.svim.vcf
    """

    stub:
    """
    touch ${sample_name}.svim.vcf
    """
}

process POSTPROCESS_SVIM {
    input:
    tuple val(sample_name), path(vcf, stageAs: "input.vcf") 
    output:
    tuple val(sample_name), path("${sample_name}.svim.vcf")

    shell:
    '''
    mkdir tmp
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
        }' !{vcf} | \
     bcftools view -i 'SVTYPE=="DEL" || SVTYPE=="INS" || SVTYPE=="DUP" || SVTYPE=="INV" || SVTYPE=="BND"' | \
     bcftools view -i 'SUPPORT>=1' | bcftools sort --temp-dir tmp -o !{sample_name}.svim.vcf
    '''

    stub:
    """
    touch ${sample_name}.svim.vcf
    """
}

process PBSV {
    publishDir "results/$sample_name", mode: "copy"

    input:
    tuple val(sample_name), path(bam), path(bai)
    output:
    tuple val(sample_name), path("${sample_name}.pbsv.vcf")
 
    script:
    if (params.sequencing_mode == "ccs_hifi" || params.sequencing_mode == "ccs_subread_fallback")
    """
    pbsv discover --sample $sample_name --tandem-repeats ${params.tandem_repeats_bed} $bam ${sample_name}.svsig.gz
    pbsv call --ccs --call-min-reads-one-sample 1 ${params.reference_fasta} ${sample_name}.svsig.gz ${sample_name}.pbsv.vcf
    """
    else
    """
    pbsv discover --sample $sample_name --tandem-repeats ${params.tandem_repeats_bed} $bam ${sample_name}.svsig.gz
    pbsv call --call-min-reads-one-sample 1 ${params.reference_fasta} ${sample_name}.svsig.gz ${sample_name}.pbsv.vcf
    """

    stub:
    """
    touch ${sample_name}.pbsv.vcf
    """
}

process POSTPROCESS_PBSV {
    input:
    tuple val(sample_name), path(vcf, stageAs: "input.vcf")
    output:
    tuple val(sample_name), path("${sample_name}.pbsv.vcf")

    shell:
    '''
    mkdir tmp
    header1='##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">'
    header2='##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">'

    cat <(bcftools view -h !{vcf} | head -n -1) \
    <(echo -e $header1) <(echo -e $header2) \
    <(bcftools view -h !{vcf} | tail -n 1) \
    <(bcftools view -H !{vcf} | \
    awk 'BEGIN{FS="\t";OFS="\t"} \
    { \
        split($NF,p,":"); \
        AD = p[2]; \
        split(AD, p_ad, ","); \
        DR = p_ad[1]; \
        DV = p_ad[2]; \
        $(NF-1) = $(NF-1)":DR:DV"; \
        $NF = $NF":"DR":"DV; \
        print $0; \
    }') | bcftools sort --temp-dir tmp -o !{sample_name}.pbsv.vcf
    '''

    stub:
    """
    touch ${sample_name}.pbsv.vcf
    """
}

process JASMINE {
    publishDir "results/$sample_name", mode: "copy", pattern: "*.merged.vcf"

    input:
    tuple val(sample_name), path(bam), path(bai), path(sniffles_vcf), path(cutesv_vcf), path(svim_vcf)
    output:
    tuple val(sample_name), path("${sample_name}.merged.vcf"), path(sniffles_vcf), path(cutesv_vcf), path(svim_vcf)
 
    shell:
    '''
    echo "!{bam.name}\n!{bam.name}\n!{bam.name}" > bam_files.txt
    echo "!{sniffles_vcf.name}\n!{cutesv_vcf.name}\n!{svim_vcf.name}" > vcf_files.txt
    mkdir tmp
    jasmine file_list=vcf_files.txt out_file=!{sample_name}.merged.vcf genome_file=!{params.reference_fasta} bam_list=bam_files.txt out_dir=tmp --output_genotypes --dup_to_ins --normalize_type --ignore_strand threads=!{task.cpus} --run_iris iris_args=--keep_long_variants,threads=!{task.cpus}
    '''

    stub:
    """
    touch ${sample_name}.merged.vcf
    """ 
}

process POSTPROCESS_JASMINE {
    publishDir "results/$sample_name", mode: "copy"

    input:
    tuple val(sample_name), path(jasmine_vcf, stageAs: "input.vcf"), path(sniffles_vcf), path(cutesv_vcf), path(svim_vcf)
    output:
    tuple val(sample_name), path("${sample_name}.merged.vcf.gz"), path("${sample_name}.merged.vcf.gz.tbi")

    shell:
    '''
    HEADER_KEYS_CUTESV="ID=CIPOS ID=CILEN ID=q5 ID=STRAND ID=RE"
    HEADER_KEYS_SVIM="ID=SUPPORT ID=STD_SPAN ID=STD_POS ID=STD_POS1 ID=STD_POS2 ID=hom_ref ID=ZMWS ID=not_fully_covered"

    bcftools view --header-lines !{jasmine_vcf} | head -n -1 > vcf_header.txt

    for key in $HEADER_KEYS_CUTESV; do
        bcftools view --header-lines !{cutesv_vcf} | grep --word-regexp $key >> vcf_header.txt
    done

    for key in $HEADER_KEYS_SVIM; do
        bcftools view --header-lines !{svim_vcf} | grep --word-regexp $key >> vcf_header.txt
    done

    bcftools view --header-lines !{jasmine_vcf} | tail -n 1 >> vcf_header.txt

    bcftools view --no-header --exclude "INFO/SVLEN==0 && INFO/SVTYPE=='INV'" !{jasmine_vcf} > variants_minus_svim_invs.txt
    bcftools view --no-header --include "INFO/SVLEN==0 && INFO/SVTYPE=='INV'" !{jasmine_vcf} > svim_invs.txt

    awk '
        BEGIN {
            FS="\t"
        }
        {
            TRUE_SVLEN = length($5)
            # print SVLEN
            gsub("SVLEN=0;", "SVLEN="TRUE_SVLEN";")
            print
        }
        END {
        }
        ' svim_invs.txt > svim_invs_with_updated_svlens.txt

    cat vcf_header.txt variants_minus_svim_invs.txt svim_invs_with_updated_svlens.txt > jasmine_merged_unsorted.vcf

    mkdir tmp
    bcftools sort --temp-dir tmp jasmine_merged_unsorted.vcf | bgzip > !{sample_name}.merged.vcf.gz
    tabix !{sample_name}.merged.vcf.gz
    '''
    // Remember to rename the samples in the VCF! (${sample_name}_cutesv, ${sample_name}_sniffles, etc.)
    stub:
    """
    touch ${sample_name}.merged.vcf.gz
    touch ${sample_name}.merged.vcf.gz.tbi
    """
}

process SPLIT_VCF_BY_CHROMOSOME {
    input:
    tuple val(sample_name), path(vcf), path(tbi)
    output:
    tuple val(sample_name), path("${sample_name}.chr*.vcf")

    script:
    """
    for n in {1..22}; do
        bcftools view --regions chr\$n $vcf > ${sample_name}.chr\$n.vcf
    done
    bcftools view --regions chrX $vcf > ${sample_name}.chrX.vcf
    bcftools view --regions chrY $vcf > ${sample_name}.chrY.vcf
    """
}

process GENOTYPE_NONTR {
    input:
    tuple val(sample_name), path(bam), path(bai), path(vcfs)
    each chromIndex

    output:
    tuple val(sample_name), path("${vcfs[chromIndex].baseName}.genotyped_nontr.vcf")

    script:
    """
    VCF_SAMPLE_NAME="\$(bcftools query -l ${vcfs[chromIndex]} | head -n +1)"
    echo "\$VCF_SAMPLE_NAME\t${bam.name}" > sample_bam.txt
    python3.8 /genSV/genotype.py nontr -v ${vcfs[chromIndex]} -o ${vcfs[chromIndex].baseName}.genotyped_nontr.vcf -s sample_bam.txt --TR_file_TRF_target ${params.TR_file_TRF_target} --TR_file_TRF_all ${params.TR_file_TRF_all} --TR_file_RM_Simpe ${params.TR_file_RM_Simpe} --TR_file_MG ${params.TR_file_MG}
    """
}

process GENOTYPE_TR {
    input:
    tuple val(sample_name), path(bam), path(bai), path(vcfs)
    each chromIndex

    output:
    tuple val(sample_name), path("${vcfs[chromIndex].baseName}.genotyped_tr.vcf")

    script:
    """
    VCF_SAMPLE_NAME="\$(bcftools query -l ${vcfs[chromIndex]} | head -n +1)"
    echo "\$VCF_SAMPLE_NAME\t${bam.name}" > sample_bam.txt
    python3.8 /genSV/genotype.py tr -t ${params.TR_file_TRF_target} -v ${vcfs[chromIndex]} -o ${vcfs[chromIndex].baseName}.genotyped_tr.vcf -s sample_bam.txt
    """
}

process MERGE_GENOTYPE_NONTR_VCFS {
    input:
    tuple val(sample_name), path(genotyped_nontr_vcfs)
    output:
    path("${sample_name}.genotyped_nontr.vcf.gz")

    script:
    """
    for n in {1..22}; do
        echo "${sample_name}.chr\$n.genotyped_nontr.vcf" >> chromosome_file_names.txt
    done
    echo "${sample_name}.chrX.genotyped_nontr.vcf" >> chromosome_file_names.txt
    echo "${sample_name}.chrY.genotyped_nontr.vcf" >> chromosome_file_names.txt

    bcftools concat --output ${sample_name}.genotyped_nontr.vcf.gz --output-type z --threads ${task.cpus} --file-list chromosome_file_names.txt 
    """
}

process SAVE_PHASED_VCFS {
    publishDir "results/$family_name", mode: "copy"
    input:
    tuple val(family_name), path(vcf), path(vcf_tbi)
    output:
    tuple val(family_name), path(vcf), path(vcf_tbi)

    script:
    """
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
    pedigree_dictionary = make_pedigree_dictionary(params.pedigree_file)

    // Create channel for sample_fastq_tuples. Each tuple is of the form: (sample_name, fastq_file_path).
    sample_fastq_tuple_channel = Channel
                                        .fromPath(params.sample_fastq_table_tsv, type: "file", checkIfExists: true)
                                        .splitCsv(sep: "\t", header: ["sample_name", "fastq_file_path"])
                                        .map { row -> tuple(row.sample_name, row.fastq_file_path) }

    // The MAP process uses minimap2 for mapping (and adding read group information) and uses samtools for sorting the output BAM files.
    // The output BAM file from MAP is then indexed in the INDEX_BAM process (using samtools index).
    sample_bam_tuple_channel = INDEX_BAM( ( MAP(sample_fastq_tuple_channel) ) )

    // Get coverage using mosdepth
    COVERAGE(sample_bam_tuple_channel)
    sample_bam_tuple_channel.view()
    // Phasing step uses pedigree structure.
    pedigree_tuple_channel = Channel
                                    .fromPath(params.pedigree_file, type: "file", checkIfExists: true)
                                    .splitCsv(sep: "\t", header: ["family_name", "sample_name", "father_name", "mother_name", "sex", "phenotype"])
                                    .map { row -> tuple(row.sample_name, row.family_name) }


    // Join the sample bam tuples with their family names, and then group samples by family (the output tuples will look like this: [ family_name, [sample_name_0, sample_name_1, sample_name_2], [bam_0, bam_1, bam_2], [bai_0, bai_1, bai_2] ]
    family_sample_bam_tuple_channel = sample_bam_tuple_channel
                                                              .join(pedigree_tuple_channel)
                                                              .map { it -> tuple(it[3], it[0], it[1], it[2]) }
    family_sample_bam_tuple_channel
                                    .groupTuple()
                                    .branch {
                                        large_families: it[1].size() >= 5
                                        standard_families: true
                                    }
                                    .set { families_to_be_phased } 

    // Deal with large families (too large to phased together) by splitting them up into trios (and they will be merged after phasing).
    families_to_be_phased.large_families
                                        .transpose()
                                        .branch {
                                            father: pedigree_dictionary[it[1]] == "Father"
                                            mother: pedigree_dictionary[it[1]] == "Mother"
                                            child: true
                                        }
                                        .set { large_family_samples}
    large_family_trios = large_family_samples.father.join(large_family_samples.mother).combine(large_family_samples.child, by: 0).map { it -> tuple( it[0], [ it[1], it[4], it[7] ], [ it[2], it[5], it[8] ], [ it[3], it[6], it[9] ] ) }
    // For phasing, if we're in joint-called SNV VCF mode, then we have to subset the VCF. Otherwise, we can skip the SUBSET_VCF process and find the family VCF corresponding to this particular family.
    // But, if we're phasing a large family, we need to subset the VCF regardless if it is joint-called or family VCF because it's too big no matter what.
    // Focus on standard families
    if (params.cohort_vcf) {
        cohort_vcf_tbi = file("${params.cohort_vcf}.tbi", checkIfExists: true)
        subsetted_family_tuples = SUBSET_VCF( large_family_trios.mix( families_to_be_phased.standard_families ), params.cohort_vcf )
    }
    else if (params.family_vcfs) {
        // Join the family VCF to the standard_families tuple 
        family_vcf_tuple_channel = Channel
                                            .fromPath(params.family_vcfs, type: "file", checkIfExists: true)
                                            .splitCsv(sep: "\t", header: ["family_name", "family_vcf_file_path"])
                                            .map { row -> tuple(row.family_name, row.family_vcf_file_path) }
        standard_subsetted_family_tuples = families_to_be_phased.standard_families.join(family_vcf_tuple_channel)
        // Still need to use SUBSET_VCF for the large families, but not sure how... (need to test this part).
        large_family_trios
                                                        .join(family_vcf_tuple_channel)
                                                        .multiMap { it ->
                                                                    bam_channel: tuple(it[0], it[1], it[2], it[3])
                                                                    vcf_channel: it[-1]
                                                                  } 
                                                        .set { large_family_trios_with_vcf }
        subsetted_large_family_trios_tuples = SUBSET_VCF( large_family_trios_with_vcf.bam_channel, large_family_trios_with_vcf.vcf_channel  )
        subsetted_family_tuples = subsetted_large_family_trios_tuples.mix( standard_subsetted_family_tuples )
    }
    else {
        log.info "Specify one of the following parameters: cohort_vcf or family_vcfs"
        exit(0)
    }
    subsetted_indexed_family_tuples = INDEX_VCF( subsetted_family_tuples )

    // Phasing
    phased_family_tuples = PHASE( subsetted_indexed_family_tuples)
    phased_family_tuples
                        .groupTuple()
                        .branch {
                            families_to_merge: it[1].size() > 1
                            standard_families: true
                        }
                        .set { families_to_merge_and_standard_families }

    // Merging the large family trio VCFs (and then use the mix operator to put all the family VCFs into one channel, all_phased_families_tuple) 
    merged_families = MERGE_FAMILY_VCFS( families_to_merge_and_standard_families.families_to_merge )
    all_phased_families_tuple_with_index = families_to_merge_and_standard_families.standard_families.transpose().mix(merged_families)

    SAVE_PHASED_VCFS( all_phased_families_tuple_with_index)

    // Haplotag each sample's BAM using its phased family VCF 
    haplotag_bam_tuple = INDEX_BAM2( HAPLOTAG( family_sample_bam_tuple_channel.combine( all_phased_families_tuple_with_index, by: 0 ) ) )

    // Run variant callers
    cutesv_tuple = POSTPROCESS_CUTESV( CUTESV(haplotag_bam_tuple) )
    sniffles_tuple = POSTPROCESS_SNIFFLES( SNIFFLES(haplotag_bam_tuple) )
    svim_tuple = POSTPROCESS_SVIM( SVIM(haplotag_bam_tuple) )
    // pbsv_tuple = POSTPROCESS_PBSV( PBSV(haplotag_bam_tuple) )

    // Merge variant calls (within each sample) using Jasmine
    jasmine_vcf_tuple = POSTPROCESS_JASMINE( JASMINE( haplotag_bam_tuple.join( sniffles_tuple ).join( cutesv_tuple ).join( svim_tuple ) ) )

    // Split each sample's VCF by chromosome
    split_vcfs_tuple = SPLIT_VCF_BY_CHROMOSOME(jasmine_vcf_tuple)

    // Genotyping variants (per chromosome)
    sample_genotype_nontr_vcf_tuple = GENOTYPE_NONTR(haplotag_bam_tuple.join(split_vcfs_tuple), chromIndices)
    sample_genotype_tr_vcf_tuple = GENOTYPE_TR(haplotag_bam_tuple.join(split_vcfs_tuple), 1)
    // Merging non-TR VCFs per sample
    MERGE_GENOTYPE_NONTR_VCFS( sample_genotype_nontr_vcf_tuple.groupTuple() )

    // To do:
    // 1. awk and bash commands should be cleaner
    // 2. Comment some code
    // 3. Add code for PBSV (with or without)
}
