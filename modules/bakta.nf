#!/usr/bin/env nextflow

process BAKTA {
    label 'bakta'
    container 'oschwengers/bakta:v1.11.3'
    publishDir "data/reproduced/prokaryote_output_tools/bakta", mode: 'copy'
    
    input:
    path fasta_file
    
    script:
    """
    bakta \
    --db ${bakta_db_dir} \
    --output ${fasta_file.baseName} \
    --prefix ${fasta_file.baseName} ${fasta_file} \
    -t ${task.cpus} \
    --skip-plot \
    --compliant
    """
}