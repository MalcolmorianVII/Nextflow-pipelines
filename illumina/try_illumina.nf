nextflow.enable.dsl=2

p = "/home/ubuntu/data/belson/projects/projects_2021/napa/input_data/clean_illumina_reads"
// Illumina reads
// params.reads = "$projectDir/*/*{1,2}.fastq.gz"
params.reads = "$p/*/*_R{1,2}.fastq.gz"
//bbduk ref
params.bbRef = "/home/ubuntu/external_tb/references/2019.04.22/adapters.fa"
// ref_genome
params.ref_genome = "/home/ubuntu/data/belson/reference/2021.04.01/17762-33892_1_71_contigs.fa"

ch = Channel.fromFilePairs(params.reads)
ch.subscribe { 
    sample,reads ->
    sample }.view()


