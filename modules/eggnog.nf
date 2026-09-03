#!/usr/bin/env nextflow
process PREPARE_EGGNOG {
    label 'eggnog_prepare'
    tag 'Prepare eggNOG database'

    output:
    path "eggnog_v5", emit: database_dir

    script:
    """
    set -euo pipefail

    mkdir -p eggnog_db

    export MAMBA_SKIP_ACTIVATE=""
    source /usr/local/bin/_activate_current_env.sh

    download_eggnog_data.py \
    -y \
    -M \
    --data_dir eggnog_db
    """
}


process EGGNOG_PROKARYOTE {
    label "eggnog_mapper_v2"
    tag "${fasta_file.baseName}"
    publishDir "data/reproduced/prokaryote_output_tools/eggnog", mode: 'copy'

    input:
    path fasta_file
    path eggnog_db

    output:
    path("${fasta_file.baseName}_eggnog.emmaper.*"), emit: eggnog_results

    script:
    """
    set -euo pipefail
    export MAMBA_SKIP_ACTIVATE=""
    export EGGNOG_DATA_DIR=/eggnog-data
    export NXT_TASK_MONITOR=0
    source /usr/local/bin/_activate_current_env.sh

    # Use a local temp directory to avoid polluting the work directory
    mkdir -p tmp
    mkdir -p ${fasta_file.baseName}_eggnog
    emmapper.py \
    --itype genome \
    --genepred prodigal \
    --decorate_gff yes \
    -i ${fasta_file} \
    -o ${fasta_file.baseName} \
    --data_dir ${eggnog_db} \
    -m mmseqs \
    --cpu ${task.cpus} \
    --output_dir ${fasta_file.baseName}_eggnog \
    --dbmem
    """

}

process EGGNOG_EUKARYOTE {
    label "eggnog_mapper_v2"
    tag "$sample_id"
    publishDir "data/reproduced/eukaryote_output_tools/eggnog", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path eggnog_db

    output:
    path("${sample_id}_eggnog.emmaper.*"), emit: eggnog_results

    script:
    """
    set -euo pipefail
    # Activate the container's environment for eggNOG-mapper.
    export MAMBA_SKIP_ACTIVATE=""
    export NXT_TASK_MONITOR=0
    source /usr/local/bin/_activate_current_env.sh

    # Use a local temp directory to avoid polluting the work directory
    mkdir -p tmp

    emapper.py \
        -i ${proteins} \
        --itype proteins \
        -o ${sample_id}_eggnog \
        -m mmseqs \
        --cpu ${task.cpus} \
        --data_dir ${eggnog_db} \
        --temp_dir ./tmp
        """
}