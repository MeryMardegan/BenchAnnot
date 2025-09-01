#!/usr/bin/env nextflow

process PROKKA {
    label 'prokka'
    container 'staphb/prokka:1.14.6'
    publishDir "results/prokka", mode: 'copy' 

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_prokka" // each prokka output directory is named after the fasta file

    script: 
    """
    prokka --outdir ${fasta_file.baseName}_prokka --prefix ${fasta_file.baseName} ${fasta_file}
    """
}

