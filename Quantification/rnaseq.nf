#!/usr/bin/env nextflow 

nextflow.enable.dsl = 2



params.raw = "./data/*{1,2}.fastq.gz"
params.outdir="./results"
params.bbduk=file("/work/users/pz192nijo/Tools/bbmap/bbduk.sh")
params.adapters=file("/work/users/pz192nijo/Tools/bbmap/resources/adapters.fa")
params.genome="/work/users/pz192nijo/Database/GenomeDB/GRCh38"
params.annot="/work/users/pz192nijo/Database/GenomeDB/GRCh38/gencode.v43.annotation.gtf"


process FASTQC {

        publishDir "${params.outdir}/fastqc"

        input:
        tuple val(sample_id), path(reads)

        output:
        file("*.{html,zip}")

        script:

        """
        fastqc -t 10 ${reads[0]} ${reads[1]}
        """
}


process REMOVE_ADAPTER {

        publishDir "${params.outdir}/clean"

        input:
        tuple val(sample_id), path(reads)


        output:
        tuple val(sample_id), path("${sample_id}{1,2}_clean.fastq.gz")

        script:

        """
        ${params.bbduk} in1=${reads[0]} in2=${reads[1]} out1=${sample_id}1_clean.fastq.gz out2=${sample_id}2_clean.fastq.gz ref=${params.adapters}  tpe tbo ktrim=r k=21 mink=11 hdist=2

        """
}


process POSTQC {

       publishDir "${params.outdir}/postqc"

       input:
       tuple val(sample_id), path(reads)

       output:
       file("*.{html,zip}")

       script:

       """
       fastqc -t 10 ${reads[0]} ${reads[1]}

       """
}


process MULTIQC {
    
    //conda 'bioconda::multiqc=1.11'
    
    publishDir "${params.outdir}/postqc", mode:'copy'

    input:
    path('*') 
    
    output:
    path('multiqc_report.html')

    script:
    """
    
    multiqc .
    
    """
}


process ALIGN {

        publishDir "${params.outdir}/bam", mode:'copy'

        input:
        tuple val(sample_id), path(reads)


        output:
        tuple val(sample_id), path("${sample_id}Aligned.sortedByCoord.out.bam")

        script:

        """
        STAR --quantMode GeneCounts --runThreadN 6 --genomeDir "${params.genome}"  --readFilesIn ${reads[0]} ${reads[1]} --readFilesCommand zcat --outFileNamePrefix ${sample_id} --outSAMtype BAM SortedByCoordinate
        """
}


process COUNT {

        publishDir "${params.outdir}/counts", mode:'copy'

        input:
        tuple val(sample_id), path(bam)


        output:
        tuple val(sample_id), path("${sample_id}.gtf"),path("${sample_id}.tab")

        script:

        """
        stringtie ${bam} -G ${params.annot} -o ${sample_id}.gtf -e -B -A ${sample_id}.tab

        """
}




reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true )

workflow {
  REMOVE_ADAPTER(reads_ch)
  POSTQC(REMOVE_ADAPTER.out)
  MULTIQC(POSTQC.out)
  ALIGN(REMOVE_ADAPTER.out)
  COUNT(ALIGN.out)
}