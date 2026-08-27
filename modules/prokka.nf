#!usr/bin/env nextflow

process PROKKA {
    label 'prokka'
    container "https://depot.galaxyproject.org/singularity/prokka%3A1.14.6--pl5321hdfd78af_5"
    publishDir "data/reproduced/prokaryote_output_tools/prokka", mode: 'copy'

    input:


    script:
    """
    prokka \
    --outdir ${fasta_file.baseName} \
    --prefix ${fasta_file.baseName} ${fasta_file} \
    --cpus ${task.cpus} \
    --compliant
    """
}