process PREPARE_KOFAM {

    label 'kofamscan_prepare'
    tag "Kofam database ${params.kofamscan_db_release}"

    input:
    val cache_dir

    output:
    path "kofam_db", emit: database_dir

    script:
    """
    set -euo pipefail

    RELEASE="${params.kofamscan_db_release}"

    CACHE_DIR="${cache_dir}"
    RELEASE_CACHE="\$CACHE_DIR/\$RELEASE"

    PROFILES_ARCHIVE="\$RELEASE_CACHE/profiles.tar.gz"
    KO_LIST_ARCHIVE="\$RELEASE_CACHE/ko_list.gz"

    BASE_URL="https://www.genome.jp/ftp/db/kofam/archives/\${RELEASE}"

    mkdir -p "\$CACHE_DIR"
    mkdir -p "\$RELEASE_CACHE"

    echo "========================================"
    echo "Preparing Kofam database"
    echo "Release: \$RELEASE"
    echo "Cache: \$CACHE_DIR"
    echo "========================================"


    # ============================================================
    # Profiles archive
    # ============================================================

    if [ -s "\$PROFILES_ARCHIVE" ] && gzip -t "\$PROFILES_ARCHIVE" 2>/dev/null; then

        echo "Using cached Kofam profiles:"
        echo "  \$PROFILES_ARCHIVE"

    else

        echo "Downloading/resuming Kofam profiles..."

        cd "\$RELEASE_CACHE"

        wget \
            -c \
            --tries=0 \
            --timeout=60 \
            --waitretry=30 \
            "\$BASE_URL/profiles.tar.gz"

        cd - >/dev/null

        gzip -t "\$PROFILES_ARCHIVE"
    fi


    # ============================================================
    # KO list archive
    # ============================================================

    if [ -s "\$KO_LIST_ARCHIVE" ] && gzip -t "\$KO_LIST_ARCHIVE" 2>/dev/null; then

        echo "Using cached Kofam KO list:"
        echo "  \$KO_LIST_ARCHIVE"

    else

        echo "Downloading/resuming Kofam KO list..."

        cd "\$RELEASE_CACHE"

        wget \
            -c \
            --tries=0 \
            --timeout=60 \
            --waitretry=30 \
            "\$BASE_URL/ko_list.gz"

        cd - >/dev/null

        gzip -t "\$KO_LIST_ARCHIVE"
    fi


    # ============================================================
    # Extract database
    # ============================================================

    echo "Extracting Kofam database..."

    rm -rf kofam_db
    mkdir -p kofam_db

    echo "Extracting Kofam profiles..."

    tar -xzf "\$PROFILES_ARCHIVE" \
        -C kofam_db

    echo "Extracting KO list..."

    gzip -dc "\$KO_LIST_ARCHIVE" \
        > kofam_db/ko_list


    # ============================================================
    # Validate extracted database
    # ============================================================

    if [ ! -d "kofam_db/profiles" ]; then
        echo "ERROR: Kofam profiles directory was not created." >&2
        exit 1
    fi

    if [ ! -s "kofam_db/ko_list" ]; then
        echo "ERROR: Kofam ko_list was not created." >&2
        exit 1
    fi

    PROFILE_COUNT=\$(find kofam_db/profiles -type f | wc -l)

    if [ "\$PROFILE_COUNT" -eq 0 ]; then
        echo "ERROR: Kofam profiles directory is empty." >&2
        exit 1
    fi

    echo "========================================"
    echo "Kofam database validated"
    echo "Release: \$RELEASE"
    echo "Profiles: \$PROFILE_COUNT"
    echo "========================================"
    """
}

process KOFAMSCAN {
    label 'kofamscan'
    tag "KofamScan annotation for $sample_id"
    publishDir "data/reproduced/eukaryote_output_tools/kofamscan", mode: 'copy'

    input:
    tuple val(sample_id), path(proteins)

    path kofam_db

    output:
    tuple val(sample_id), path("${sample_id}.kofam.txt"), emit: results

    script:
    """
    set -euo pipefail

    rm -rf kofam_tmp
    mkdir -p kofam_tmp

    /usr/local/bin/exec_annotation \
      --cpu ${task.cpus} \
      --profile /database/profiles/eukaryote.hal \
      --ko-list /database/ko_list \
      --tmp-dir "\$PWD/kofam_tmp" \
      -f detail-tsv \
      --report-unannotated \
      ${proteins} \
      -o ${sample_id}.kofam.txt
    """
}
