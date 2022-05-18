nextflow.enable.dsl=2

workflow {
    //Get long read
    reads_ch = Channel.fromPath(params.nanoRead, checkIfExists:true )
    //First run stats
    ASSEMBLY_STATS(reads_ch)
    //Run kraken2
    KRAKEN2(reads_ch)
    //Run flye
    FLYE(reads_ch)
    //Run racon
    RACON(reads_ch,FLYE.out)
    //Run medaka
    MEDAKA(reads_ch,RACON.out)
    //Run bakta
    BAKTA(MEDAKA.out)

}

process ASSEMBLY_STATS {
    tag "Perform QC on long read"
    publishDir "$projectDir/assemblyStat",mode:"copy",overwrite: false

    conda '/home/belson/miniconda3/envs/assembly-stat'

    input:
    path (read)

    output:
    path "assemblyStats.tsv" 

    script:
    """
    assembly-stats <(gzip -cd $read) | tee assemblyStats.tsv
    """
}

process KRAKEN2 {
    tag "Perform metagenomic sampling"
    publishDir "$projectDir/kraken2",mode:"copy",overwrite: false

    conda "/home/belson/miniconda3/envs/kraken2"

    input:
    path (read)

    output:
    path 'kraken_report.txt'
    path 'kraken'

    script:
    """
    kraken2 --use-names --threads 8 --db $params.krakenDb --report kraken_report.txt --gzip-compressed $read > kraken
    """
}

process FLYE {
    tag "Genome assembly"
    publishDir "$projectDir/Flye",mode:"copy",overwrite: false

    conda "/home/belson/miniconda3/envs/flye"

    input:
    path (read)

    output:
    path 'Flye'

    script:
    """
    flye --nano-hq $read -g 5m -o Flye -t 8  
    """
}



process RACON{
    tag "Perform initial polishing"
    publishDir "$projectDir/racon",mode:"copy",overwrite: false

    conda "/home/belson/miniconda3/envs/racon"

    input:
    path (read) 
    path 'Flye'

    output:
    file "racon.paf"
    file "racon.fasta"

    script:
    """
    minimap2 -x map-ont Flye/assembly.fasta $read > racon.paf
    racon -t 4 $read racon.paf Flye/assembly.fasta > racon.fasta
    """
}

process MEDAKA{
    tag "Medaka miracle genome polishing"
    publishDir "$projectDir/Medaka",mode:"copy",overwrite: false

    conda "/home/belson/miniconda3/envs/medaka"

    input:
    path (read)
    file "racon.paf"
    file "racon.fasta"

    output:
    path 'Medaka'

    script:
    """
    medaka_consensus -i $read -d racon.fasta -t 8  -m $params.medakaMod -o Medaka
    """
}

process BAKTA{
    tag "Perform genome annotation"
    publishDir "$projectDir/Bakta",mode:"copy",overwrite: false

    conda "/home/belson/.conda/envs/bakta"

    input:
    path 'Medaka'

    output:
    path 'Bakta'

    script:
    """
    bakta --db $params.baktaDb Medaka/consensus.fasta -o Bakta
    """
}

workflow.onComplete { 
	log.info ( workflow.success ? "\n<<<<<Pipeline run successfully @ M7>>>>>" : "Oops .. something went wrong" )
}