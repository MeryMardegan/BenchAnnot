process KOFAMSCAN {
    label 'kofamscan'
    tag "KofamScan annotation for $sample_id"
    publishDir "/data/reproduced/eukaryote_output_tools/kofamscan", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)

    path kofam_db

    output:
    tuple val(sample_id), path("${sample_id}.kofam.txt"), emit: results

    script:
    """
    set -euo pipefail

    /usr/local/bin/exec_annotation \
      --cpu ${task.cpus} \
      --profile /database/profiles \
      --ko-list /database/ko_list \
      -f detail-tsv \
      --report-unannotated \
      ${proteins} \
      -o ${sample_id}.kofam.txt
    """
}