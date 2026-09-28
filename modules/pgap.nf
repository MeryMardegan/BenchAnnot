#!/usr/bin/env nextflow

process PREPARE_PGAP_PYTHON {

    label 'pgap_python_prepare'
    tag "Python ${params.pgap_python_version} for PGAP"

    publishDir "${projectDir}/data/database", mode: 'copy'

    input:
    val runtime_dir

    output:
    val runtime_dir, emit: runtime

    script:
    """
    set -euo pipefail

    RUNTIME="${runtime_dir}"
    RUNTIME_ROOT=\$(dirname "\$RUNTIME")
    MAMBA_ROOT="\$RUNTIME_ROOT/micromamba"
    BOOTSTRAP="\$PWD/micromamba-bootstrap"
    RUNTIME_TMP="\$RUNTIME_ROOT/tmp"

    mkdir -p "\$RUNTIME_ROOT"
    mkdir -p "\$BOOTSTRAP"
    mkdir -p "\$MAMBA_ROOT"
    mkdir -p "\$RUNTIME_TMP"

    export TMPDIR="\$RUNTIME_TMP"
    export MAMBA_ROOT_PREFIX="\$MAMBA_ROOT"


    # Reuse a complete runtime to keep preparation idempotent.
    if [ -x "\$RUNTIME/bin/python" ]; then
        "\$RUNTIME/bin/python" -c \
            'import sys; assert sys.version_info[:2] == (3, 11)'
        "\$RUNTIME/bin/python" --version
        exit 0
    fi

    # Remove partial state before recreating the runtime.
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

    publishDir "${projectDir}/data/database", mode: 'copy', overwrite: true

    input:
    val pgap_python_dir

    output:
    path "pgap", emit: runtime

    script:
    """
    set -euo pipefail

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
    publishDir "${params.outdir}/prokaryote_output_tools/pgap", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta_file), val(species)

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

    export PGAP_INPUT_DIR="${pgap_dir}"

    "${pgap_python_dir}/bin/python" \
        "${pgap_dir}/pgap.py" \
        -n \
        -g ${fasta_file} \
        -s "${species}" \
        --taxcheck \
        --auto-correct-tax \
        -o ${sample_id}_pgap \
        -D singularity \
        --container-path "${pgap_container}" \
        --no-internet
    """
}
