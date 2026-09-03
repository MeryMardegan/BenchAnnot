#!/usr/bin/env nextflow

process PGAP {
label "pgap"
tag "PGAP annotation for ${sample_id}"
publishDir "data/reproduced/prokaryote_output_tools/pgap", mode: 'copy'

input:
tuple val(sample_id),
      path(fasta_file),
      val(species)

val pgap_dir
val pgap_container

output:
tuple val(sample_id),
      path("${sample_id}_pgap"),
      emit: results

script:
"""
set -euo pipefail

export PGAP_INPUT_DIR=${pgap_dir}

\${PGAP_INPUT_DIR}/pgap.py \
-n \
-g ${fasta_file} \
-s "${species}" \
--taxcheck \
--auto-correct-tax \
-o ${sample_id}_pgap \
-D singularity \
--container-path ${pgap_container} \
--no-internet

"""
}