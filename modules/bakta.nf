#!/usr/bin/env nextflow


process BAKTA {
    label 'bakta'
    container 'oschwengers/bakta:v1.11.3'
    publishDir "results/bakta", mode: 'copy'
    containerOptions '--entrypoint ""'// necessary for the way bakta image is built
    

    input:
    path fasta_file
    path bakta_db_dir

    output:
    path "${fasta_file.baseName}_bakta" // each bakta output directory is named after the fasta file

    script: 
    """
    bakta --db ${bakta_db_dir} --output ${fasta_file.baseName}_bakta --prefix ${fasta_file.baseName} ${fasta_file}
    """
}

