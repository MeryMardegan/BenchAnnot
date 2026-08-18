process KOFAMSCAN {
    label 'kofamscan'
    publishDir { 'results/eukaryotes/kofamscan' }, mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)
    path (profiles)
    path (ko_list)

    output:
    path "${sample_id}.kofam.txt"

    script:
    """
    /usr/local/bin/exec_annotation \
      -o ${sample_id}.kofam.txt \
      --profile ${profiles} \
      --ko-list ${ko_list} \
      -f detail-tsv --cpu ${task.cpus} \
      ${proteins}
    """
}