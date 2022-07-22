params.reference_fasta = file("/expanse/projects/sebat1/tjena/LongReadAnalysisPipeline/Homo_sapiens_assembly38.fasta", type: "file", checkIfExists: true)
params.reference_fasta_fai = file("${params.reference_fasta}.fai", checkIfExists: true)

params.joint_called_vcf = file("/expanse/projects/sebat1/tjena/LongReadAnalysisPipeline/Ashkenazi/HG.vcf.gz", type: "file", checkIfExists: true)
params.joint_called_vcf_tbi = file("${params.joint_called_vcf}.tbi", checkIfExists: true)

params.trBedFile = file("/expanse/projects/sebat1/tjena/LongReadAnalysisPipeline/human_GRCh38_no_alt_analysis_set.trf.bed", type: "file", checkIfExists: true)
params.genSVtrBedFile = file("/expanse/projects/sebat1/tjena/LongReadAnalysisPipeline/Simple_Repeats_TRF_annot.bed", type: "file")
params.pedFile = file("AshkenazimTrio.ped", type: "file", checkIfExists: true)
params.sample_files = file("sample_files_Ashkenazi.tsv", type: "file", checkIfExists: true)
//params.TMP_DIR = "/scratch/$USER/job_$SLURM_JOBID"
params.TMP_DIR = "work_dir"

if (params.reference_fasta) { reference_fasta = file(params.reference_fasta, type: "file", checkIfExists: true) }
if (params.joint_called_vcf) { joint_called_vcf = file(params.joint_called_vcf, type: "file", checkIfExists: true) }
if (params.genSVtrBedFile) { genSVtrBedFile = file(params.genSVtrBedFile, type: "file", checkIfExists: true) }
if (params.trBedFile) { trBedFile = file(params.trBedFile, type: "file", checkIfExists: true) }
if (params.pedFile) { pedFile = file(params.pedFile, type: "file", checkIfExists: true) }
if (params.sample_files) { sample_files = file(params.sample_files, type: "file", checkIfExists: true) }
if (params.TMP_DIR) { TMP_DIR = params.TMP_DIR }

println "Project : $workflow.projectDir"
params.outdir = '/expanse/projects/sebat1/tjena/LongReadAnalysisPipeline/Ashkenazi/pipeline_output'
outdir= params.outdir ?: '/expanse/projects/sebat1/tjena/LongReadAnalysisPipeline/Ashkenazi/pipeline_output'
println "output : $params.outdir"

params.useSnifflesV1 = true
if (params.useSnifflesV1) { useSnifflesV1 = params.useSnifflesV1 }
params.hifi = false
if (params.hifi) { hifi = params.hifi }

params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    testing_sv_workflow_for_long_reads: a SV calling workflow
    ============================
    --referencs_fasta	reference file, input to minimap2.
    --joint_called_vcf	file to vcf file to be used to phase with long read bam files
    --genSVtrBedFile	TR bed file, in a 6 column format as expected by genotype.py for TR analysis
    --trBedFile		TR Bed file, input to pbsv discover
    --pedFile		ped file to use to phase
    --sample_files	file listing sample names and path to where the read fastq.gz files are located
    Other arguments:
    ----------------
    --outdir		path to where the pipeline outputs will be saved.
    --hifi		true if input files are PB hifi reads. Parameters for mapping and SV calling depend on this
    --useSnifflesV1     use sniffles v1 (as opposed to sniffles v2)
    Execution environment dependency:
    ---------------------------------
    process CUTESV and JASMINE assume a SDSC slurm environment to create a temporary TMP_DIR directory in the 
    execution directory. This is preferred but not necessary; it is used to read/write temporary files which may
    not be allowed at high frequencies unless at a designated location. If not in SDSC expanse, you should modify 
    the path to a location recommended by the system administrators or simply delete it in the script. This will put
    the temporary files in the execution directory.
    -----------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

chromNum = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]

ped = channel
            .fromPath("$pedFile", type: "file", checkIfExists: true)
            .splitCsv(sep: "\t", header: ["family_id", "sample_id", "father_id", "mother_id", "sex", "phenotype"])
            .map { row -> tuple(row.sample_id, row.family_id) }

// The data channel will output tuples, each contains a sample ID, the path to the sample's input bam file, and the sample's family ID. 
sample_fastq_tuple = channel
	    .fromPath("$sample_files", type: "file", checkIfExists: true) 
            .splitCsv(sep: "\t", header: ["sample_id", "reads_file"])
            .map { row -> tuple(row.sample_id, row.reads_file) } 
            .join(ped)
            .map { it -> tuple(it[2], it[0], it[1]) }

process MINIMAP2 {
    conda '/home/ux453059/miniconda3/envs/pbenv'
 
    input: 
    tuple val(family_id), val(sample_id), path(fastq_gz)
 
    output: 
    tuple val(family_id), val(sample_id), path("${fastq_gz.simpleName}.sam")
 
    script:
      if(params.hifi) 
      """ 
      minimap2 -t ${task.cpus} -H -x map-hifi -a --MD -Y -o ${fastq_gz.simpleName}.sam $reference_fasta $fastq_gz
      """
    else 
      """
      minimap2 -t ${task.cpus} -x map-ont -a --MD -Y -o ${fastq_gz.simpleName}.sam $reference_fasta $fastq_gz
      """
 
    stub: 
    """
    touch ${fastq_gz.simpleName}.sam 
    """ 
} 

process SORT {
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(family_id), val(sample_id), path(sam)

    output:
    tuple val(family_id), val(sample_id), path("${sam.baseName}.bam")

    script:
    """
    samtools sort --reference $reference_fasta -@ ${task.cpus} -O bam -o ${sam.baseName}.bam $sam
    """
    stub:
    """
    touch ${sam.baseName}.bam
    """
}

process ADD_READGROUP {
    publishDir params.outdir, mode:'copy', overwrite:true
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(family_id), val(sample_id), path(bam)

    output:
    tuple val(family_id), val(sample_id), path("${bam.baseName}.withSM.bam")

    script:
    """
    samtools addreplacerg -@ ${task.cpus} -r "@RG\tID:$sample_id\tSM:$sample_id" -O bam -o ${bam.baseName}.withSM.bam $bam
    """

    stub:
    """
    touch ${bam.baseName}.withSM.bam
    """
} 

process INDEX {
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

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
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(family_id), val(sample_names), path(bams), path(bais)
    path(joint_called_vcf)

    output:
    tuple val(family_id), val(sample_names), path(bams), path(bais), path("${family_id}.vcf.gz"), path("${family_id}.vcf.gz.csi")

    script:
    """
    bcftools view --threads ${task.cpus} --samples ${sample_names.join(",")} --output-type z --output ${family_id}.vcf.gz $joint_called_vcf
    bcftools index --threads ${task.cpus} ${family_id}.vcf.gz
    """

    stub:
    """
    touch "${family_id}.vcf.gz"
    touch "${family_id}.vcf.gz.csi"
    """
}

process PHASE {
    conda '/home/ux453059/miniconda3/envs/pbenv'

    input:
    tuple val(family_id), val(sample_names), path(bams), path(bais), path(family_vcf), path(family_vcf_csi)

    output:
    tuple val(family_id), path("${family_id}.phased.vcf.gz")

    script:
    """
    whatshap phase --reference $reference_fasta --ped $pedFile --indels --tag PS --output ${family_id}.phased.vcf.gz $family_vcf $bams
    """

    stub:
    """
    touch ${family_id}.phased.vcf.gz
    touch ${family_id}.phased.vcf.gz.tbi
    """
}

process INDEX_PHASE { 
    conda '/home/ux453059/miniconda3/envs/sambcfenv' 
 
    input: 
    tuple val(family_id), path(family_vcf) 
 
    output: 
    tuple val(family_id), path(family_vcf), path("${family_vcf}.csi") 
 
    script: 
    """
    bcftools index --threads ${task.cpus} $family_vcf
    """ 
 
    stub: 
    """ 
    touch ${family_vcf} 
    touch "${family_vcf}.csi" 
    """ 
}

process MERGE_FAMILY_VCFS {
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(family_id), path(vcfs, stageAs: "?.vcf.gz"), path(vcf_csis, stageAs: "?.vcf.gz.csi")
    output:
    tuple val(family_id), path("${family_id}.phased.vcf.gz"), path("${family_id}.phased.vcf.gz.csi")
    script:
    """
    bcftools merge --threads ${task.cpus} --force-samples --output "${family_id}.phased.vcf.gz" --output-type z $vcfs
    bcftools index --threads ${task.cpus} "${family_id}.phased.vcf.gz"
    """
    stub:
    """
    touch ${family_id}.phased.vcf.gz
    """
}

process HAPLOTAG {
    publishDir params.outdir, mode:'copy', overwrite:true
    conda '/home/ux453059/miniconda3/envs/pbenv'

    input:
    tuple val(family_id), val(sample_name), path(bam), path(bai), path(phased_family_vcf), path(phased_family_vcf_csi)

    output:
    tuple val(sample_name), path("${bam.simpleName}.haplotag.bam"), path("${bam.simpleName}.haplotag.bam.bai")

    script:
    """
    whatshap haplotag --reference $reference_fasta --sample $sample_name --ignore-read-groups --tag-supplementary --output-haplotag-list ${bam.simpleName}_haplotag_list.tab.gz --output ${bam.simpleName}.haplotag.bam $phased_family_vcf $bam
    samtools index -@ ${task.cpus} ${bam.simpleName}.haplotag.bam
    """

    stub:
    """
    touch ${bam.simpleName}.haplotag.bam
    touch ${bam.simpleName}.haplotag.bam.bai
    """
}

process CUTESV {
    conda '/home/ux453059/miniconda3/envs/swimenv'

    input:
    tuple val(sampleName),path(map_sort_bam_cutesv),path(map_sort_bam_index_cutesv) 

    output:
    tuple val(sampleName),path("${map_sort_bam_cutesv.simpleName}.mm2.cutesv.s1.vcf")

    script:
    if(params.hifi)
      """
      mkdir tmp
      cuteSV -t ${task.cpus} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.9 --report_readid \
        --min_support 1 --genotype ${map_sort_bam_cutesv} ${reference_fasta} ${map_sort_bam_cutesv.simpleName}.mm2.cutesv.s1.vcf tmp
      """
    else
      """
      mkdir tmp
      cuteSV -t ${task.cpus} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --report_readid \
        --min_support 1 --genotype ${map_sort_bam_cutesv} ${reference_fasta} ${map_sort_bam_cutesv.simpleName}.mm2.cutesv.s1.vcf tmp
      """

    stub:
    """
    touch ${map_sort_bam_cutesv.simpleName}.mm2.cutesv.s1.vcf
    """
}

process POSTPROCESS_CUTESV {
    publishDir params.outdir, mode:'copy', overwrite:true
    conda '/home/ux453059/miniconda3/envs/sambcfenv'
    
    input:
    tuple val(sampleName),path(cutesv_output_raw)
    
    output:
    tuple val(sampleName),path("${cutesv_output_raw.simpleName}.mm2.cutesv.s1.fixed.vcf")
    
    script:
    """
    bcftools view -i 'POS>0' ${cutesv_output_raw} | bcftools sort -o ${cutesv_output_raw.simpleName}.mm2.cutesv.s1.fixed.vcf
    """    
    
    stub:
    """
    touch ${cutesv_output_raw.simpleName}.mm2.cutesv.s1.fixed.vcf
    """
}

process SNIFFLES {
    conda '/home/ux453059/miniconda3/envs/pbenv'
 
    input:
    tuple val(sampleName),path(map_sort_bam_sniffles),path(map_sort_bam_index_sniffles) 
 
    output:
    tuple val(sampleName),path("${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf")
 
    script:
    if(params.hifi && params.useSnifflesV1)
      """
      sniffles --ccs_reads -t ${task.cpus} --num_reads_report -1 --tmp_file tmp1 --min_support 1 --cluster -m $map_sort_bam_sniffles \
        -v ${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf
      """
    else if(!params.hifi && params.useSnifflesV1)
      """
      sniffles -t ${task.cpus} --num_reads_report -1 --tmp_file tmp1 --min_support 1 --cluster -m $map_sort_bam_sniffles \
        -v ${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf
      """

    stub:
    """
    touch ${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf
    """
}

process SNIFFLES_V2 {
    conda '/home/ux453059/miniconda3/envs/swimenv'

    input:
    tuple val(sampleName),path(map_sort_bam_sniffles),path(map_sort_bam_index_sniffles)

    output:
    tuple val(sampleName),path("${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf")

    script:
    """
    sniffles -t ${task.cpus} --minsupport 1 --reference $reference_fasta --tandem-repeats $trBedFile -i $map_sort_bam_sniffles -v ${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf
    """

    stub:
    """
    touch ${map_sort_bam_sniffles.simpleName}.mm2.sniffles.s1.vcf
    """
}

process POSTPROCESS_SNIFFLES {
    publishDir params.outdir, mode:'copy', overwrite:true
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(sampleName),path(sniffles_output_raw)

    output:
    tuple val(sampleName),path("${sniffles_output_raw.simpleName}.mm2.sniffles.s1.fixed.vcf") 

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
    conda '/home/ux453059/miniconda3/envs/swimenv'
 
    input:
    tuple val(sampleName),path(map_sort_bam_svim),path(map_sort_bam_index_svim) 
 
    output:
    tuple val(sampleName),path("${map_sort_bam_svim.simpleName}.mm2.svim.vcf") 
 
    script:
    """
    svim alignment --read_names --zmws out_svim/  $map_sort_bam_svim $reference_fasta
    mv out_svim/variants.vcf ${map_sort_bam_svim.simpleName}.mm2.svim.vcf
    """

    stub:
    """
    touch ${map_sort_bam_svim.simpleName}.mm2.svim.vcf
    """    
}

process POSTPROCESS_SVIM {
    publishDir params.outdir, mode:'copy', overwrite:true
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(sampleName),path(svim_output_raw) 

    output:
    tuple val(sampleName),path("${svim_output_raw.simpleName}.mm2.svim.s1.fixed.vcf") 

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
    touch ${svim_output_raw.simpleName}.mm2.svim.s1.fixed.vcf
    """
}

process PBSV {
    conda '/home/ux453059/miniconda3/envs/swimenv'
 
    input:
    tuple val(sampleName),path(map_sort_bam_pbsv),path(map_sort_bam_index_pbsv) 
 
    output:
    tuple val(sampleName),path("${map_sort_bam_pbsv.simpleName}.mm2.pbsv.s1.vcf") 
    
    script:
    if(params.hifi)
      """
      pbsv discover --sample $sampleName --tandem-repeats $trBedFile $map_sort_bam_pbsv \
        ${map_sort_bam_pbsv.simpleName}.svsig.gz
      pbsv call --ccs -O 1 $reference_fasta ${map_sort_bam_pbsv.simpleName}.svsig.gz ${map_sort_bam_pbsv.simpleName}.mm2.pbsv.s1.vcf
      """
    else
      """
      pbsv discover --sample $sampleName --tandem-repeats $trBedFile $map_sort_bam_pbsv \
        ${map_sort_bam_pbsv.simpleName}.svsig.gz
      pbsv call -O 1 $reference_fasta ${map_sort_bam_pbsv.simpleName}.svsig.gz ${map_sort_bam_pbsv.simpleName}.mm2.pbsv.s1.vcf
      """

    stub:
    """
    touch ${map_sort_bam_pbsv.simpleName}.mm2.pbsv.s1.vcf
    """
}

process POSTPROCESS_PBSV {
    publishDir params.outdir, mode:'copy', overwrite:true
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(sampleName),path(pbsv_output_raw) 

    output:
    tuple val(sampleName),path("${pbsv_output_raw.simpleName}.mm2.pbsv.s1.fixed.vcf") 

    shell:
    '''
    header1='##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">'
    header2='##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">'

    cat <(bcftools view -h !{pbsv_output_raw} | head -n -1) \
    <(echo -e $header1) <(echo -e $header2) \
    <(bcftools view -h !{pbsv_output_raw} | tail -n 1) \
    <(bcftools view -H !{pbsv_output_raw} | \
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
    }') | bcftools sort -o !{pbsv_output_raw.simpleName}.mm2.pbsv.s1.fixed.vcf
    '''

    stub:
    """
    touch ${pbsv_output_raw.simpleName}.mm2.pbsv.s1.fixed.vcf
    """
}

process JASMINE {
    conda '/home/ux453059/miniconda3/envs/swimenv'

    input:
    tuple val(sampleName),path(map_sort_bam),path(map_sort_bam_index),path(vcfs_sn),path(vcfs_cu),path(vcfs_pb),path(vcfs_sv)

    output:
    tuple val(sampleName),path("Jasmine_Iris_merge_scps_haptag_${sampleName}.vcf") 
 
    shell:
    '''
    OUT_FILE=Jasmine_Iris_merge_scps_haptag_!{sampleName}.vcf
    echo "!{map_sort_bam.name}\n!{map_sort_bam.name}\n!{map_sort_bam.name}\n!{map_sort_bam.name}" > bam_files_scps.txt
    echo "!{vcfs_sn.name}\n!{vcfs_cu.name}\n!{vcfs_pb.name}\n!{vcfs_sv.name}" > vcf_files_scps.txt
    BAMS_LIST=bam_files_scps.txt
    VCFS_FILE_LIST=vcf_files_scps.txt

    jasmine file_list=$VCFS_FILE_LIST out_file=$OUT_FILE genome_file=!{reference_fasta} bam_list=$BAMS_LIST out_dir=!TMP_DIR --output_genotypes \
      --dup_to_ins --normalize_type --ignore_strand threads=!{task.cpus} --run_iris iris_args=--keep_long_variants,threads=!{task.cpus}
    '''

    stub:
    """
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.vcf
    """ 
}

process POSTPROCESS_JASMINE {
    publishDir params.outdir, mode:'copy', overwrite:true
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(sampleName),path(jasmine_output_raw),path(cutesv_vcf),path(pbsv_vcf),path(svim_vcf)

    output:
    tuple val(sampleName),path("${jasmine_output_raw.simpleName}.s1.vcf.gz"),path("${jasmine_output_raw.simpleName}.s1.vcf.gz.tbi") 

    shell:
    '''
    #!/bin/bash

    HEADER_KEYS_CSV="ID=CIPOS ID=CILEN ID=q5 ID=STRAND"
    HEADER_KEYS_PBSV="ID=SVANN ID=MATEID ID=MATEDIST ID=NearReferenceGap ID=NearContigEnd ID=Decoy ID=InsufficientStrandEvidence"
    HEADER_KEYS_SVIM="ID=SUPPORT ID=STD_SPAN ID=STD_POS ID=STD_POS1 ID=STD_POS2 ID=SEQS ID=hom_ref ID=ZMWS ID=not_fully_covered"

    vcf_temp=temp.vcf

    bcftools view -h !{jasmine_output_raw} | head -n -1 > $vcf_temp
    for key in $HEADER_KEYS_CSV; do
        bcftools view -h !{cutesv_vcf} | grep -w $key >> $vcf_temp
    done
    for key in $HEADER_KEYS_PBSV; do
        bcftools view -h !{pbsv_vcf} | grep -w $key >> $vcf_temp
    done
    for key in $HEADER_KEYS_SVIM; do
        bcftools view -h !{svim_vcf} | grep -w $key >> $vcf_temp
    done
    bcftools view -h !{jasmine_output_raw} | tail -n 1 >> $vcf_temp
    bcftools view -H !{jasmine_output_raw} >> $vcf_temp

    bcftools sort $vcf_temp | bgzip > !{jasmine_output_raw.simpleName}.s1.vcf.gz
    tabix !{jasmine_output_raw.simpleName}.s1.vcf.gz
    '''

    stub:
    """
    touch ${jasmine_output_raw.simpleName}.s1.vcf.gz
    touch ${jasmine_output_raw.simpleName}.s1.vcf.gz.tbi
    """
}

process SEP_CHROMS_PRE_GENOTYPE {
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(sampleName),path(jasmine_output),path(jasmine_output_index) 

    output:
    tuple val(sampleName),path("${jasmine_output.simpleName}.chr*.vcf") 

    shell:
    '''
    IN_VCF=!{jasmine_output}
    IN_FILE=!{jasmine_output.simpleName}

    for n in {1..22}; do
        CHR="chr$n"
        OUT_VCF=$IN_FILE"."$CHR".vcf"
        bcftools view --threads !{task.cpus} -r $CHR $IN_VCF > $OUT_VCF
    done


    CHR="chrX"
    OUT_VCF=$IN_FILE"."$CHR".vcf"
    bcftools view --threads !{task.cpus} -r $CHR $IN_VCF > $OUT_VCF

    CHR="chrY"
    OUT_VCF=$IN_FILE"."$CHR".vcf"
    bcftools view --threads !{task.cpus} -r $CHR $IN_VCF > $OUT_VCF
    '''

    stub:
    """
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.chr1.vcf
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.chr2.vcf
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.chr3.vcf
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.chr10.vcf
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.chr11.vcf
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.chr20.vcf
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.chr21.vcf
    touch Jasmine_Iris_merge_scps_haptag_${sampleName}.chrX.vcf
    """
}

process GENOTYPE_NONTR {
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(sampleName),path(sep_chrs),path(bamsample),path(bamsample_index) 
    each chrmNum

    output:
    tuple val(sampleName),path("Jasmine_Iris_merge*chr*genotyped.nonTR.vcf") 

    shell:
    '''
    echo "0_!{sampleName}\t!{bamsample.name}" > sample_bam.txt
    python /home/ux453059/genSV/genotype.py nontr -v !{sep_chrs[chrmNum]} -o !{sep_chrs[chrmNum].baseName}.genotyped.nonTR.vcf -s sample_bam.txt
    '''

    stub:
    """
    touch ${sep_chrs[chrmNum].baseName}.genotyped.nonTR.vcf
    """
}

process GENOTYPE_TR {
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(sampleName),path(sep_chrs),path(bamsample),path(bamsample_index)
    each chrmNum

    output:
    tuple val(sampleName),path("Jasmine_Iris_merge*chr*genotyped.TR.vcf")

    shell:
    '''
    echo "0_!{sampleName}\t!{bamsample.name}" > sample_bam.txt
    python /home/ux453059/genSV/genotype.py tr -t !{genSVtrBedFile} -v !{sep_chrs[chrmNum]} -o !{sep_chrs[chrmNum].baseName}.genotyped.TR.vcf -s sample_bam.txt
    '''

    stub:
    """
    touch ${sep_chrs[chrmNum].baseName}.genotyped.TR.vcf
    """
}

process MERGE_GENOTYPE_CHROMS_NONTR {
    publishDir params.outdir, mode:'copy', overwrite:true 
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(sampleName),path(gtyped_chrs) 
    output:
    path("${gtyped_chrs[0].simpleName}.s1.genotyped.nonTR.vcf.gz") 

    shell:
    '''
    #!/bin/bash

    mkdir tmp
    TMP_DIR=tmp
    bcftools view --threads !{task.cpus} -h !{gtyped_chrs[0]} > !{gtyped_chrs[0].simpleName}.s1.genotyped.nonTR.vcf

    VCF_FILES=Jasmine_Iris_merge*chr*genotyped.nonTR.vcf

    for vcf_file in $VCF_FILES; do
            bcftools view --threads !{task.cpus} -H $vcf_file >> !{gtyped_chrs[0].simpleName}.s1.genotyped.nonTR.vcf
    done

    bcftools sort -T $TMP_DIR !{gtyped_chrs[0].simpleName}.s1.genotyped.nonTR.vcf | bgzip > !{gtyped_chrs[0].simpleName}.s1.genotyped.nonTR.vcf.gz
    tabix !{gtyped_chrs[0].simpleName}.s1.genotyped.nonTR.vcf.gz
    '''

    stub:
    """
    touch ${gtyped_chrs[0].simpleName}.s1.genotyped.nonTR.vcf.gz
    """
}

process MERGE_GENOTYPE_CHROMS_TR {
    publishDir params.outdir, mode:'copy', overwrite:true
    conda '/home/ux453059/miniconda3/envs/sambcfenv'

    input:
    tuple val(sampleName),path(gtyped_chrs)
    output:
    path("${gtyped_chrs[0].simpleName}.s1.genotyped.TR.vcf.gz")

    shell:
    '''
    #!/bin/bash

    mkdir tmp
    TMP_DIR=tmp
    bcftools view --threads !{task.cpus} -h !{gtyped_chrs[0]} > !{gtyped_chrs[0].simpleName}.s1.genotyped.TR.vcf

    VCF_FILES=Jasmine_Iris_merge*chr*genotyped.TR.vcf

    for vcf_file in $VCF_FILES; do
            bcftools view --threads !{task.cpus} -H $vcf_file >> !{gtyped_chrs[0].simpleName}.s1.genotyped.TR.vcf
    done

    bcftools sort -T $TMP_DIR !{gtyped_chrs[0].simpleName}.s1.genotyped.TR.vcf | bgzip > !{gtyped_chrs[0].simpleName}.s1.genotyped.TR.vcf.gz
    tabix !{gtyped_chrs[0].simpleName}.s1.genotyped.TR.vcf.gz
    '''

    stub:
    """
    touch ${gtyped_chrs[0].simpleName}.s1.genotyped.TR.vcf.gz
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
    pedigree_dictionary = make_pedigree_dictionary("$pedFile")
    // Mapping steps (uses minimap2 for alignment to create a sam file, which is then converted to a bam file, and then index and sort the output bam file)
    sample_bam_tuple = INDEX(ADD_READGROUP(SORT(MINIMAP2( sample_fastq_tuple ))))

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
    subsetted_family_tuples = SUBSET_VCF( large_family_trios.mix( families_to_be_phased.standard_family ), joint_called_vcf )
    phased_family_tuples = INDEX_PHASE(PHASE( subsetted_family_tuples ))
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

    // Haplotag each sample's bam file using its phased family VCF
    haplotag_bam_tuple = HAPLOTAG( sample_bam_tuple.combine( all_phased_families_tuple_with_index, by: 0 ) )

    // Run variant callers
    cutesv_tuple = POSTPROCESS_CUTESV(CUTESV( haplotag_bam_tuple ))
    if(params.useSnifflesV1) {
	sniffles_tuple = POSTPROCESS_SNIFFLES(SNIFFLES( haplotag_bam_tuple ))
    } else {
	sniffles_tuple = POSTPROCESS_SNIFFLES(SNIFFLES_V2( haplotag_bam_tuple ))
    }
    svim_tuple = POSTPROCESS_SVIM(SVIM( haplotag_bam_tuple ))
    pbsv_tuple = POSTPROCESS_PBSV( PBSV(haplotag_bam_tuple) )

    // Merge variant callers (within each sample) using Jasmine
    jasmine_vcf = JASMINE( haplotag_bam_tuple.join( sniffles_tuple ).join( cutesv_tuple ).join( pbsv_tuple ).join( svim_tuple ) )
    jasmine_merged_tuple = POSTPROCESS_JASMINE( jasmine_vcf.join(cutesv_tuple).join(pbsv_tuple).join(svim_tuple) )

    // Extract chromosomes, genotype them separately then merge to get final genotyped vcf
    // Needs to be done separately for nonTR and TR regions
    chromosome_vcfs = SEP_CHROMS_PRE_GENOTYPE(jasmine_merged_tuple)

    genotyped_chroms_nonTR = GENOTYPE_NONTR( chromosome_vcfs.join(haplotag_bam_tuple), chromNum )
    genotyped_nonTR = MERGE_GENOTYPE_CHROMS_NONTR( genotyped_chroms_nonTR.groupTuple() )

    genotyped_chroms_TR = GENOTYPE_TR( chromosome_vcfs.join(haplotag_bam_tuple), chromNum )
    genotyped_TR = MERGE_GENOTYPE_CHROMS_TR( genotyped_chroms_TR.groupTuple() )
}

