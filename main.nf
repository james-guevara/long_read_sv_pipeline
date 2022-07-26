params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    long_read_sv_pipeline
    ============================
    Required:
    --sample_fastq_table    Tab-delimited file, where the first column corresponds to the sample name and the second column corresponds to the file path of the fastq reads input.
    --reference_fasta       Reference genome to which the reads are aligned.
    --joint_called_vcf      Joint-called (or merged) VCF with all the samples in the cohort.
    --pedigree_file         Pedigree file (.fam, .ped, or .psam) containing all the samples in the cohort.
    --sequencing_mode       The sequencing technology used to generate the input reads. Must be one of: ccs_hifi, ccs_subread_fallback, or ont
    """.stripIndent()
    exit 0
}

params.sample_fastq_table_tsv = file("sample_files_Ashkenazi_tests.tsv", type: "file", checkIfExists: true)
params.reference_fasta = file("/expanse/projects/sebat1/tjena/LongReadAnalysisPipeline/Homo_sapiens_assembly38.fasta", type: "file", checkIfExists: true)
params.reference_fasta_fai = file("${params.reference_fasta}.fai", type: "file", checkIfExists: true)
params.pedigree_file = file("AshkenazimTrio.ped", type: "file", checkIfExists: true)
params.sequencing_mode = "ccs_hifi"
sequencing_modes = ["ccs_hifi", "ccs_subread_fallback", "ont"]
if (sequencing_modes.contains(params.sequencing_mode) == false) {
    log.info "Sequencing mode parameter must be one of: ccs_hifi, ccs_subread_fallback, or ont"
    exit(0)
}

process MAP {
    input: 
    tuple val(sample_name), path(fastq_gz)
    output: 
    tuple val(sample_name), path("${sample_name}.bam")
 
    script:
    if (params.sequencing_mode == "ccs_hifi") 
    """ 
    minimap2 -t ${task.cpus} -o ${sample_name}.sam -a -x map-hifi --MD -Y  ${params.reference_fasta} $fastq_gz | samtools sort --reference ${params.reference_fasta} --threads ${task.cpus} --output-fmt BAM -o ${sample_name}.bam

    """
    else if (params.sequencing_mode == "ccs_subread_fallback") // I might have to change the -Q parameter (though I don't see it in the most recent version of minimap2)
    """
    minimap2 -t ${task.cpus} -a -x map-hifi --MD -Y  ${params.reference_fasta} $fastq_gz | samtools sort --reference ${params.reference_fasta} --threads ${task.cpus} --output-fmt BAM -o ${sample_name}.bam

    """
    else if (params.sequencing_mode == "ont")
    """
    minimap2 -t ${task.cpus} -a -x map-ont  --MD -Y  ${params.reference_fasta} $fastq_gz | samtools sort --reference ${params.reference_fasta} --threads ${task.cpus} --output-fmt BAM -o ${sample_name}.bam
    """
 
    stub: 
    """
    touch ${sample_name}.bam 
    """ 
}

process ADD_READGROUP {
    input:
    tuple val(sample_name), path(bam, stageAs: "input.bam")
    output:
    tuple val(sample_name), path("${sample_name}.bam")

    script:
    """
    samtools addreplacerg --threads ${task.cpus} -r "@RG\tID:$sample_name\tSM:$sample_name" --output-fmt BAM -o ${sample_name}.bam $bam
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
    input:
    tuple val(sample_name), path(bam), path(bai)

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


workflow {
    sample_fastq_tuple_channel = Channel
                                        .fromPath(params.sample_fastq_table_tsv, type: "file", checkIfExists: true)
                                        .splitCsv(sep: "\t", header: ["sample_name", "fastq_file_path"])
                                        .map { row -> tuple(row.sample_name, row.fastq_file_path) }
    // Use minimap2 for mapping, and samtools for adding read group in the header and for sorting and indexing the BAM file
    sample_bam_tuple_channel = INDEX_BAM( ADD_READGROUP( MAP(sample_fastq_tuple_channel) ) )

    // Get coverage using mosdepth
    COVERAGE(sample_bam_tuple_channel)

    // Phasing step uses pedigree structure (though it's still optional)
    pedigree_tuple_channel = Channel
                                    .fromPath(params.pedigree_file, type: "file", checkIfExists: true)
                                    .splitCsv(sep: "\t", header: ["family_name", "sample_name", "father_name", "mother_name", "sex", "phenotype"])
                                    .map { row -> tuple(row.sample_name, row.family_name) }

    // Join the sample bam tuples with their family names, and then group samples by family (the output tuples will look like this: [ family_name, [sample_name_0, sample_name_1, sample_name_2], [bam_0, bam_1, bam_2], [bai_0, bai_1, bai_2] ]
    family_sample_bam_tuple_channel = sample_bam_tuple_channel
                                                              .join(pedigree_tuple_channel)
                                                              .map { it -> tuple(it[3], it[0], it[1], it[2]) }
                                                              .groupTuple().view()
                                                              .branch {
                                                                  large_families: it[1].size() >= 5
                                                                  standard_families: true
                                                              }
                                                              .set { families_to_be_phased } 
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




}
