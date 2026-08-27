#!/usr/bin/env nextflow

process EGGNOG_PROKARYOTE {
    label "eggnog_prokaryote"
    tag "$sample_id"
    container "brunoholiva/eggnog_v2:2.1.13-sse41"
    publishDir "data/reproduced/prokaryote_output_tools/eggnog", mode: 'copy'

    input:
    path fasta_file

    script:
    """
    set -euo pipefail
    export MAMBA_SKIP_ACTIVATE=""
    source
    /usr/local/bin/_activate_current_env.sh
    mkdir ${fasta_file.baseName}_eggnog
    emmaper.py \
    --itype genome \
    --genepred prodigal \
    --decorate_gff yes \
    -i ${fasta_file} \
    -o ${fasta_file.baseName} \
    --data_dir ${eggnog_db_dir} \
    -m mmseqs \
    --cpu ${task.cpus} \
    --output_dir ${fasta_file.baseName}_eggnog \
    --dbmen
    """

}

process EGGNOG_EUKARYOTE {
    label "eggnog_eukaryote"
    tag "$sample_id"
    publishDir "data/reproduced/eukaryote_output_tools/eggnog", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)

    output:
    path("${sample_id}.eggnog.emmaper.*")

    script:
    """
    set -euo pipefail
    # Activate the container's environment for eggNOG-mapper.
    export MAMBA_SKIP_ACTIVATE=""
    source /usr/local/bin/_activate_current_env.sh

    # Use a local temp directory to avoid polluting the work directory
    mkdir -p tmp
    emapper.py \
    -i ${proteins} \
    --itype proteins \
    -o ${sample_id}_eggnog \
    -m mmseqs \
    --cpus ${task.cpus} \
    --data_dir /eggnog-data \
    --temp_dir ./tmp
    """
}