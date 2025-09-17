process GFFREAD {
  label 'gffread'
  tag "$sample_id"
  publishDir "results/eukaryotes/gffread", mode: 'copy'

  input:
  tuple val(sample_id), path(fasta), path(anno)

  output:
  tuple val(sample_id), path("${sample_id}.faa")

  script:
  """
  set -euo pipefail

  # remover transcritos/linhas com trans-splicing e strand indefinido ('?'), mantendo cabeçalhos (#), e qualquer linha com strand +/-
  awk 'BEGIN{FS=OFS="\\t"}
       /^#/ {print; next}
       \$7=="?" {next}
       \$9 ~ /exception=trans-splicing/ {next}
       {print}' "${anno}" > "${sample_id}.filtered.gff"

  # Métrica simples do filtro
  total=\$(grep -vc '^#' "${anno}" || true)
  kept=\$(grep -vc '^#' "${sample_id}.filtered.gff" || true)
  removed=\$(( total - kept ))
  echo "[GFFREAD/${sample_id}] total=\$total kept=\$kept removed=\$removed" >&2

  # 2) Se não sobrou nada útil
  if [ "\$kept" -le 0 ]; then
    echo "[GFFREAD/${sample_id}] No usable entries after filtering (all were trans-splicing or strand '?')." >&2
    exit 2
  fi

  # 3) Execução tolerante do gffread para extrair proteína
  gffread "${sample_id}.filtered.gff" \
    -g "${fasta}" \
    -y "${sample_id}.faa" \
    --force-exons --gene2exon --t-adopt -E
  """
}
