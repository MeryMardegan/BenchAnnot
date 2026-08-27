process KOFAMSCAN {
    label 'kofamscan'
    tag "$sample_id"
    publishDir "results/eukaryotes/kofamscan", mode: 'copy', pattern: "${sample_id}.kofam.*"

    input:
    tuple val(sample_id), path(proteins)
    path (profiles)
    path (ko_list)

    output:
    path "${sample_id}.kofam.txt"

    script:
    """
    # Run KO assignment using the provided profile and KO list databases.
    /usr/local/bin/exec_annotation \
      --cpu ${task.cpus} \
      --profile ${profiles} \
      --ko-list ${ko_list} \
      -f detail-tsv \
      -f mapper-one-line \
      --report-unannotated \
      ${proteins} \
      -o ${sample_id}.kofam.txt
    """