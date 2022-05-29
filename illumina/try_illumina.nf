nextflow.enable.dsl=2

p = "/home/ubuntu/data/belson/projects/projects_2021/napa/input_data/clean_illumina_reads"
// Illumina reads
// params.reads = "$projectDir/*/*{1,2}.fastq.gz"
params.reads = "$p/*/*_R{1,2}.fastq.gz"

// ch = Channel.fromFilePairs(params.reads)
// ch.view()

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)

    FASTQC(reads_ch)
}

process FASTQC{
    tag "Perform read QC"
    publishDir "$projectDir/${sample}",mode:"copy"
    conda "/home/belson/miniconda3/envs/fastqc"

    input:
    tuple val(sample),path(reads)

    output:
    path( "*" ), emit: fastqc_out

    script:
    """
    mkdir fastqc
    fastqc ${reads} -o fastqc
    """

}

