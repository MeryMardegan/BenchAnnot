process INTERPROSCAN {
    label 'interproscan'
    tag "$sample_id"
    publishDir { 'results/eukaryotes/interproscan' }, mode: 'copy', pattern: "${sample_id}.*"

    input:
    // O GFFREAD emite a tuple val(sample_id), path("${sample_id}.faa") e por consequencia vamos usar ele
    tuple val(sample_id), path(faa)

    output:
    // O output padrão do InterProScan é nomeado com base no input, mas aqui renomeamos para um padrão estável do módulo
    tuple val(sample_id),
          path ("${sample_id}.interpro.tsv"),
          path ("${sample_id}.interpro.gff3")

    script:
    def outbase = "${sample_id}.interpro"
    def inputFa = params.ips_test ? "${sample_id}.subset.faa" : faa
    def fmt = (params.ips_formats ?: 'tsv,gff3')
    """
    set -euo pipefail

    mkdir -p temp

    /opt/interproscan/interproscan.sh \
      -i ${inputFa} \
      -f ${fmt} \
      -cpu ${task.cpus} \
      -goterms \
      -b ${outbase} \
      --tempdir temp

    """
}