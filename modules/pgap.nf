#!/usr/bin/env nextflow

process PREPARE_PGAP_PYTHON {

    label 'pgap_python_prepare'
    tag "Python ${params.pgap_python_version} for PGAP"

    input:
    val runtime_dir

    output:
    val runtime_dir, emit: runtime

    script:
    """
    set -euo pipefail

    RUNTIME="${runtime_dir}"
    MAMBA_ROOT="/data/runtime/micromamba"
    BOOTSTRAP="\$PWD/micromamba-bootstrap"

    mkdir -p "/data/runtime"
    mkdir -p "\$BOOTSTRAP"

    # If a valid runtime already exists, reuse it.
    if [ -x "\$RUNTIME/bin/python" ]; then
        "\$RUNTIME/bin/python" -c \
            'import sys; assert sys.version_info[:2] == (3, 11)'
        "\$RUNTIME/bin/python" --version
        exit 0
    fi

    # Remove an incomplete installation from a previous failed attempt.
    rm -rf "\$RUNTIME"

    curl -Ls \
        https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj -C "\$BOOTSTRAP" bin/micromamba

    export MAMBA_ROOT_PREFIX="\$MAMBA_ROOT"

    "\$BOOTSTRAP/bin/micromamba" create \
        -y \
        -p "\$RUNTIME" \
        -c conda-forge \
        "python=${params.pgap_python_version}"

    "\$RUNTIME/bin/python" --version

    "\$RUNTIME/bin/python" -c \
        'import sys; assert sys.version_info[:2] == (3, 11)'
    """
}

process PREPARE_PGAP {

    label 'pgap_prepare'
    tag "PGAP ${params.pgap_version}"

    input:
    val pgap_python_dir

    output:
    path "pgap", emit: runtime

    script:
    """
    set -euo pipefail

    mkdir -p pgap

    mkdir -p pgap

    curl -L \
        https://raw.githubusercontent.com/ncbi/pgap/${params.pgap_version}/scripts/pgap.py \
        -o pgap/pgap.py

    chmod +x pgap/pgap.py

    cd pgap

    mkdir -p tmp

    export PGAP_INPUT_DIR="\$PWD"
    export TMPDIR="\$PWD/tmp"
    export SINGULARITY_TMPDIR="\$PWD/tmp"

    "${pgap_python_dir}/bin/python" pgap.py \
        --use-version ${params.pgap_version} \
        --no-self-update \
        -D singularity

    rm -rf tmp
    """
}

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
val pgap_python_dir

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