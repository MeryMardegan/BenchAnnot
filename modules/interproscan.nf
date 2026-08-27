process INTERPROSCAN {
    label 'interproscan'
    tag "$sample_id"
    publishDir "results/eukaryotes/interproscan", mode: 'copy', pattern: "${sample_id}.*"

    input:
    // Input comes from GFFREAD: tuple(val(sample_id), path("${sample_id}.faa")).
    tuple val(sample_id), path(faa)

    output:
    // Standardize output names to a stable module prefix.
    tuple val(sample_id),
          path ("${sample_id}.interpro.*")

    script:
    def outbase = "${sample_id}.interpro"
    def inputFa = faa
    def fmt = (params.ips_formats ?: 'tsv,gff3')

    """
    set -euo pipefail
    # Keep temporary files scoped to the task directory.
    mkdir -p temp

    /opt/interproscan/interproscan.sh \
      -i ${inputFa} \
      -f ${fmt} \
      -cpu ${task.cpus} \
      -goterms \
      --iprlookup \
      --pathways \
      -b ${outbase} \
      --tempdir temp
    """