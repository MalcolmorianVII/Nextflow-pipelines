nextflow.enable.dsl=2

p = "/home/ubuntu/data/belson/projects/projects_2021/napa/input_data/clean_illumina_reads"
// Illumina reads
params.reads = "$projectDir/*/*{1,2}.fastq.gz"
//bbduk ref
params.bbRef = "/home/ubuntu/external_tb/references/2019.04.22/adapters.fa"
// ref_genome
params.ref_genome = "/home/ubuntu/data/belson/reference/2021.04.01/17762-33892_1_71_contigs.fa"
workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)

    FASTQC(reads_ch)
    BBDUK(reads_ch)
    MULTIQC(BBDUK.out)
    PREPARE-PHENIX(params.ref_genome)
    PHENIX(BBDUK.out,reads_ch)
    SKESA(BBDUK.out)
    ASSEMBLY-STATS(SKESA.out)
}

FASTQC{
    tag "Perform read QC"
    publishDir "$projectDir/$sample",mode:"copy"
    conda "/home/belson/miniconda3/envs/fastqc"

    input:
    tuple val(sample),path(reads)

    output:
    path( "*_fastqc*" ), emit: fastqc_out

    script:
    """
    fastqc ${reads}
    """

}

BBDUK{
    tag "Perform read trimming"
    publishDir "$projectDir/$sample",mode:"copy"
    conda "/home/belson/miniconda3/envs/bbduk"

    input:
    tuple val(sample),path(reads)

    output:
    path trimmed

    script:
    """
    bbduk.sh threads=8 in=${reads[0]} in2=${reads[1]} out=${sample}_trimmed_r1.fastq.gz out2=${sample}_trimmed_r2.fastq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=20 minlength=50
    """
}

MULTIQC{
    tag "Perform collective read QC"
    publishDir "$projectDir/$sample",mode:"copy"
    conda "/home/belson/miniconda3/envs/multiqc"

    input:
    path trimmed

    output:
    file "{qc_results_dir}/multiqc_report.html"

    script:
    """
    multiqc -o {qc_results_dir} $trimmed
    """
}

PREPARE-PHENIX{
    tag "Prepare reference genome"
    publishDir "$projectDir/$sample",mode:"copy"
    conda "/home/belson/.conda/envs/phenix"

    input:
    path ref

    output:
    'ref.fai'

    script:
    """
    phenix.py prepare_reference --mapper bwa --variant gatk --reference $ref
    """
}

PHENIX{
    tag "Perform read QC"
    publishDir "$projectDir/$sample",mode:"copy"
    conda "/home/belson/.conda/envs/phenix"

    input:
    path trimmed
    file ref_genome
    output:
    path "$sample/Phenix"
    script:
    """
    phenix.py run_snp_pipeline.py -r1 ${trimmed[0]} -r2 ${trimmed[1]} -r $ref_genome --sample-name $sample --mapper bwa --variant gatk --filters min_depth:5,mq_score:30
    """
}

SKESA{
   tag "Perform read QC"
    publishDir "$projectDir/$sample",mode:"copy"
    conda "/home/belson/miniconda3/envs/skesa"
    input:
    path trimmed
    output:
    path "$sample/skesa/${sample}_skesa.fa"
    script:
    """
    skesa --fasta $trimmed --cores 4 --memory 48 > $sample/skesa/${sample}_skesa.fa
    """
}

ASSEMBLY-STATS{
    tag "Perform read QC"
    publishDir "$projectDir/$sample",mode:"copy"
    conda "/home/belson/miniconda3/envs/assembly-stat"
    input:
    path "$sample/skesa/${sample}_skesa.fa"
    output:
    path "$sample/assembly-stat"
    script:
    """
    assembly-stats -t $sample/skesa > $sample/assembly-stat/${sample}_contigs.assembly_stats.tsv
    """
}