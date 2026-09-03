#!/usr/bin/env nextflow

process BAKTA {
    label 'bakta'
    tag "Bakta annotation for ${sample_id}"
    publishDir "${params.outdir}/prokaryote_output_tools/bakta", mode: 'copy'

    input:
    tuple val (sample_id), path(fasta_file)
    path bakta_db

    output:
    tuple val (sample_id), path("results/${sample_id}.*"), emit: results

    script:
    """
    set -euo pipefail

    bakta \
        --db /database \
        --output results \
        --prefix ${sample_id} \
        --threads ${task.cpus} \
        --skip-plot \
        --compliant \
        ${fasta_file}
    """
}