nextflow.enable.dsl=2

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)

    FASTQC(reads_ch)
    BBDUK(reads_ch)
    SKESA(reads_ch,BBDUK.out)
    ASSEMBLY_STATS(reads_ch,SKESA.out)
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

process BBDUK{
    tag "Perform read trimming"
    publishDir "$projectDir/${sample}/bbduk",mode:"copy"
    conda "/home/belson/miniconda3/envs/bbduk"

    input:
    tuple val(sample),path(reads)

    output:
    path "${sample}_trimmed_r1.fastq.gz"
    path "${sample}_trimmed_r2.fastq.gz"

    script:
    """
    bbduk.sh threads=8 in=${reads[0]} in2=${reads[1]} out=${sample}_trimmed_r1.fastq.gz out2=${sample}_trimmed_r2.fastq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=20 minlength=50
    """
}


process SKESA{
    tag "Perform read assembly"
    publishDir "$projectDir/${sample}/skesa",mode:"copy"
    conda "/home/belson/miniconda3/envs/skesa"

    input:
    tuple val(sample),path(reads)
    path "${sample}_trimmed_r1.fastq.gz"
    path "${sample}_trimmed_r2.fastq.gz"

    output:
    file "${sample}_skesa.fa"

    script:
    """
    skesa --fastq ${sample}_trimmed_r1.fastq.gz ${sample}_trimmed_r2.fastq.gz --cores 4 --memory 48 > ${sample}_skesa.fa
    """
}

process ASSEMBLY_STATS{
    tag "Perform assembly QC"
    publishDir "$projectDir/${sample}/assembly-stat",mode:"copy"
    conda "/home/belson/miniconda3/envs/assembly-stat"

    input:
    tuple val(sample),path(reads)
    file "${sample}_skesa.fa"

    output:
    file "${sample}_contigs.assembly_stats.tsv"

    script:
    """
    assembly-stats -t ${sample}_skesa.fa > ${sample}_contigs.assembly_stats.tsv
    """
}