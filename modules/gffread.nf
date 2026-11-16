process GFFREAD {
    label 'gffread'
    tag "GFFread extraction for ${sample_id}"
    publishDir "${params.outdir}/eukaryote_output_tools/gffread", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta), path(anno)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.gff"), emit: filtered_gff
    tuple val(sample_id), path("${sample_id}_gffread.faa"), emit: proteins

    script:
    """
    set -euo pipefail

    # Preserve headers while excluding unsupported trans-splicing and strand records.
    awk 'BEGIN{FS=OFS="\\t"}
        /^#/ {print; next}
        \$7=="?" {next}
        \$9 ~ /exception=trans-splicing/ {next}
        {print}' "${anno}" > "${sample_id}.filtered.gff"

    # Report filtering counts for traceability.
    total=\$(grep -vc '^#' "${anno}" || true)
    kept=\$(grep -vc '^#' "${sample_id}.filtered.gff" || true)
    removed=\$(( total - kept ))
    echo "[GFFREAD/${sample_id}] total=\$total kept=\$kept removed=\$removed" >&2

    # Fail before sequence extraction if no usable feature remains.
    if [ "\$kept" -le 0 ]; then
      echo "[GFFREAD/${sample_id}] No usable entries after filtering (all were trans-splicing or strand '?')." >&2
      exit 2
    fi
    
    gffread -F -S -C -J \
      "${sample_id}.filtered.gff" \
      -g "${fasta}" \
      -y "${sample_id}_gffread.faa"
    """
}
