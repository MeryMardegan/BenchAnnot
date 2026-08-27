process GFFREAD {
  label 'gffread'
  tag "$sample_id"
  publishDir "data/reproduced/eukaryote_output_tools/gffread", mode: 'copy'

  input:
  tuple val(sample_id), path(fasta), path(anno)

  output:
  tuple val(sample_id), path("${sample_id}.faa")

  script:
  """
  set -euo pipefail

  # Filter out trans-splicing entries and undefined strand ('?') while preserving headers.
    awk 'BEGIN{FS=OFS="\\t"}
        /^#/ {print; next}
        \$7=="?" {next}
        \$9 ~ /exception=trans-splicing/ {next}
        {print}' "${anno}" > "${sample_id}.filtered.gff"

    # Emit a simple filtering summary for traceability.
    total=\$(grep -vc '^#' "${anno}" || true)
    kept=\$(grep -vc '^#' "${sample_id}.filtered.gff" || true)
    removed=\$(( total - kept ))
    echo "[GFFREAD/${sample_id}] total=\$total kept=\$kept removed=\$removed" >&2

    # Fail early if the filter removed all entries.
    if [ "\$kept" -le 0 ]; then
      echo "[GFFREAD/${sample_id}] No usable entries after filtering (all were trans-splicing or strand '?')." >&2
      exit 2
    fi
    
    # Run gffread to extract protein sequences from the filtered GFF.
    gffread -F -S -C -J \
      "${sample_id}.filtered.gff" \
      -g "${fasta}" \
      -y "${sample_id}.faa" \
    """
}
