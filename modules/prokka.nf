#!/usr/bin/env nextflow

process PROKKA {
    label 'prokka'
    tag "Prokka annotation for ${sample_id}"
    publishDir "${params.outdir}/prokaryote_output_tools/prokka", mode: 'copy', saveAs: { filename -> file(filename).name }

    input:
    tuple val (sample_id), path(fasta_file)

    output:
    tuple val (sample_id), path("results/${sample_id}.*"), emit: results

    script:
    """
    set -euo pipefail

    prokka \
        --outdir results \
        --prefix ${sample_id} \
        --cpus ${task.cpus} \
        --compliant \
        ${fasta_file}
    """
}