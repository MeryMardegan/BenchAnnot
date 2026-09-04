#!/usr/bin/env nextflow

process PREPARE_BAKTA {

    label 'bakta_prepare'
    tag 'Bakta database download'

    output:
    path "bakta_db", emit: database_dir

    script:
    """
    set -euo pipefail

    mkdir -p download

    bakta_db download \
        --output download \
        --type full

    db_dir=\$(find download \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -name 'db*' \
        | head -n 1)

    if [ -z "\${db_dir}" ]; then
        echo "Bakta database directory was not found after download." >&2
        exit 1
    fi

    mv "\${db_dir}" bakta_db
    """
}

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