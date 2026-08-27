#!/usr/bin/env nextflow

process PGAP {
label "pgap"
publishDir "data/reproduced/prokaryote_output_tools/pgap", mode: 'copy'

input:
path fasta_file
val pgap_dir


output:
path "${fasta_file.baseName}_pgap"

script:
"""
export PGAP_INPUT_DIR=${pgap_dir}
echo \${PGAP_INPUT_DIR}
\${PGAP_INPUT_DIR}/pgap.py \
-n \
-g ${fasta_file} \
-s "Mycoplasmoides genitalium" \
--taxcheck \
--auto-correct-tax \
-o ${fasta_file.baseName} \
-D singularity \
--container-path ${pgap_dir}/pgap_2025-05-06.build7983.sif \
--no-internet

"""
}