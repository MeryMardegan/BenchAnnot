process INTERPROSCAN {
    label 'interproscan'
    tag "InterProScan annotation for $sample_id"
    publishDir "data/reproduced/eukaryote_output_tools/interproscan", mode: 'copy'

    input:
    // Input comes from GFFREAD: tuple(val(sample_id), path("${sample_id}.faa")).
    tuple val(sample_id), path(faa)
    path ips_db

    output:
    // Standardize output names to a stable module prefix.
    tuple val(sample_id),
          path ("${sample_id}.interpro.*"),
          emit: results

    script:
    def outbase = "${sample_id}.interpro"
    def fmt = (params.ips_formats ?: 'tsv,gff3')

    """
    set -euo pipefail
    # Keep temporary files scoped to the task directory.

    mkdir -p temp

    /opt/interproscan/interproscan.sh \
      -i ${faa} \
      -f ${fmt} \
      -cpu ${task.cpus} \
      -goterms \
      --iprlookup \
      --pathways \
      -b ${outbase} \
      --tempdir temp
    """
}