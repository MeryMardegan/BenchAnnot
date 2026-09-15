#!/usr/bin/env nextflow

process PREPARE_BAKTA {

    label 'bakta_prepare'
    tag "Bakta database v${params.bakta_db_version} download"

    input:
    val cache_dir

    output:
    path "bakta_db", emit: database_dir

    script:
    """
    set -euo pipefail

    # ============================================================
    # Paths
    # ============================================================

    CACHE_DIR="${cache_dir}"

    TARBALL="\$CACHE_DIR/db.tar.xz"
    PARTIAL="\${TARBALL}.part"

    ZENODO_URL="https://zenodo.org/record/14916843/files/db.tar.xz"
    BACKUP_URL="https://s3.computational.bio.uni-giessen.de/bakta-db/db-v${params.bakta_db_version}.tar.xz"

    mkdir -p download
    mkdir -p tmp

    export TMPDIR="\$PWD/tmp"
    export TMP="\$TMPDIR"
    export TEMP="\$TMPDIR"

    echo "========================================"
    echo "Preparing Bakta database"
    echo "Bakta DB version: ${params.bakta_db_version}"
    echo "AMRFinder DB: ${params.bakta_amrfinder_db_version}"
    echo "Cache: \$CACHE_DIR"
    echo "========================================"


    # ============================================================
    # Bakta database archive
    # ============================================================

    if [ -s "\$TARBALL" ] && \
       echo "${params.bakta_db_md5}  \$TARBALL" | md5sum -c -
    then

        echo "Using cached Bakta DB archive: \$TARBALL"

    else

        rm -f "\$TARBALL"

        echo "Downloading Bakta DB ${params.bakta_db_version} from Zenodo..."

        if ! wget \
            -c \
            --tries=3 \
            --timeout=60 \
            --waitretry=30 \
            -O "\$PARTIAL" \
            "\$ZENODO_URL"
        then

            echo "Zenodo download failed. Trying Bakta backup repository..." >&2

            wget \
                -c \
                --tries=0 \
                --timeout=60 \
                --waitretry=30 \
                -O "\$PARTIAL" \
                "\$BACKUP_URL"
        fi

        echo "Validating Bakta database archive..."

        echo "${params.bakta_db_md5}  \$PARTIAL" \
            | md5sum -c -

        mv "\$PARTIAL" "\$TARBALL"
    fi


    # ============================================================
    # Extract Bakta DB
    # ============================================================

    echo "Extracting Bakta database..."

    rm -rf download/db

    tar -xJf "\$TARBALL" \
        -C download


    # ============================================================
    # Validate Bakta DB using Bakta itself
    # ============================================================

    echo "Validating Bakta database..."

    python -c 'from pathlib import Path; from bakta.db import check; info=check(Path("download/db")); expected="${params.bakta_db_version}"; detected=str(info["major"])+"."+str(info["minor"]); assert detected == expected, "Bakta DB version mismatch: expected="+expected+", detected="+detected; assert info["type"] == "full", "Bakta DB type mismatch: expected=full, detected="+str(info["type"]); print("Bakta database validated: "+detected+" ("+str(info["type"])+")")'


    # ============================================================
    # Prepare fixed AMRFinderPlus database
    # ============================================================

    AMR_ROOT="download/db/amrfinderplus-db"
    AMR_VERSION="${params.bakta_amrfinder_db_version}"

    AMR_DIR="\$AMR_ROOT/\$AMR_VERSION"

    AMR_URL="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/${params.bakta_amrfinder_db_series}/${params.bakta_amrfinder_db_version}"

    echo "Preparing AMRFinderPlus DB \$AMR_VERSION..."

    # Remove the version bundled with the Bakta archive.
    # BenchAnnot installs the explicitly pinned version instead.
    rm -rf "\$AMR_ROOT"

    mkdir -p "\$AMR_DIR"

    echo "Downloading fixed AMRFinderPlus database:"
    echo "\$AMR_URL"

    wget \
        --recursive \
        --level=1 \
        --no-parent \
        --no-host-directories \
        --cut-dirs=6 \
        --reject="index.html*,robots.txt" \
        --execute robots=off \
        --tries=10 \
        --timeout=60 \
        --waitretry=30 \
        --directory-prefix="\$AMR_DIR" \
        "\$AMR_URL/"


    # ============================================================
    # Validate downloaded AMRFinder DB
    # ============================================================

    if [ ! -s "\$AMR_DIR/version.txt" ]; then
        echo "ERROR: AMRFinderPlus version.txt was not downloaded." >&2
        exit 1
    fi

    if ! grep -Fq "\$AMR_VERSION" "\$AMR_DIR/version.txt"; then
        echo "ERROR: unexpected AMRFinderPlus database version." >&2
        echo "Expected: \$AMR_VERSION" >&2
        echo "Detected:" >&2
        cat "\$AMR_DIR/version.txt" >&2
        exit 1
    fi


    # ============================================================
    # Index AMRFinder DB
    # ============================================================

    echo "Indexing AMRFinderPlus database..."

    amrfinder_index "\$AMR_DIR"

    # Bakta/AMRFinder expects this symlink.
    ln -sfn \
        "\$AMR_VERSION" \
        "\$AMR_ROOT/latest"


    # ============================================================
    # Validate AMRFinder installation
    # ============================================================

    if [ ! -L "\$AMR_ROOT/latest" ]; then
        echo "ERROR: AMRFinderPlus latest symlink was not created." >&2
        exit 1
    fi

    DETECTED_AMR_VERSION=\$(basename "\$(readlink -f "\$AMR_ROOT/latest")")

    if [ "\$DETECTED_AMR_VERSION" != "\$AMR_VERSION" ]; then
        echo "ERROR: AMRFinderPlus database version mismatch." >&2
        echo "Expected: \$AMR_VERSION" >&2
        echo "Detected: \$DETECTED_AMR_VERSION" >&2
        exit 1
    fi

    echo "AMRFinderPlus database validated:"
    echo "  version: \$DETECTED_AMR_VERSION"


    # ============================================================
    # Final Bakta DB validation
    # ============================================================

    echo "Final Bakta database validation..."

    python -c 'from pathlib import Path; from bakta.db import check; info=check(Path("download/db")); print("Final Bakta database validation successful: "+str(info["major"])+"."+str(info["minor"])+" ("+str(info["type"])+")")'


    # ============================================================
    # Final Nextflow output
    # ============================================================

    mv download/db bakta_db

    echo "========================================"
    echo "Bakta database preparation complete"
    echo "Bakta DB: ${params.bakta_db_version}"
    echo "AMRFinder DB: ${params.bakta_amrfinder_db_version}"
    echo "========================================"
    """
}

process BAKTA {
    label 'bakta'
    tag "Bakta annotation for ${sample_id}"
    publishDir "${params.outdir}/prokaryote_output_tools/bakta", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta_file)
    path bakta_db

    output:
    tuple val(sample_id),
          path("results/${sample_id}.*"),
          emit: results

    script:
    """
    set -euo pipefail

    mkdir -p tmp
    mkdir -p matplotlib-cache

    export TMPDIR="\$PWD/tmp"
    export TMP="\$TMPDIR"
    export TEMP="\$TMPDIR"
    export MPLCONFIGDIR="\$PWD/matplotlib-cache"

    bakta \
        --db /database \
        --output results \
        --prefix ${sample_id} \
        --threads ${task.cpus} \
        --tmp-dir "\$PWD/tmp" \
        --skip-plot \
        --compliant \
        --debug \
        ${fasta_file}
    """
}
