process PREPARE_INTERPROSCAN {

    label 'interproscan_prepare'
    tag "InterProScan data ${params.interproscan_version}"

    input:
    val cache_dir

    output:
    path "interproscan_data", emit: database_dir

    script:
    """
    set -euo pipefail

    VERSION="${params.interproscan_version}"

    CACHE_DIR="${cache_dir}"

    ARCHIVE="\$CACHE_DIR/interproscan-data-\${VERSION}.tar.gz"
    PARTIAL="\${ARCHIVE}.part"

    MD5_FILE="\$CACHE_DIR/interproscan-data-\${VERSION}.tar.gz.md5"

    BASE_URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/\${VERSION}/alt"

    mkdir -p "\$CACHE_DIR"
    mkdir -p extracted

    echo "========================================"
    echo "Preparing InterProScan database"
    echo "Version: \$VERSION"
    echo "Cache: \$CACHE_DIR"
    echo "========================================"


    # ============================================================
    # Download archive
    # ============================================================

    if [ -s "\$ARCHIVE" ] && [ -s "\$MD5_FILE" ]; then

        echo "Checking cached InterProScan archive..."

    else

        echo "Downloading InterProScan data..."

        if command -v wget >/dev/null 2>&1; then

            wget \
                -c \
                --tries=0 \
                --timeout=60 \
                --waitretry=30 \
                -O "\$PARTIAL" \
                "\$BASE_URL/interproscan-data-\${VERSION}.tar.gz"

            wget \
                --tries=0 \
                --timeout=60 \
                --waitretry=30 \
                -O "\$MD5_FILE" \
                "\$BASE_URL/interproscan-data-\${VERSION}.tar.gz.md5"

        else

            python -c 'import urllib.request; urllib.request.urlretrieve("'"\\\$BASE_URL"'/interproscan-data-'"\\\$VERSION"'.tar.gz", "'"\\\$PARTIAL"'")'

            python -c 'import urllib.request; urllib.request.urlretrieve("'"\\\$BASE_URL"'/interproscan-data-'"\\\$VERSION"'.tar.gz.md5", "'"\\\$MD5_FILE"'")'
        fi

        mv "\$PARTIAL" "\$ARCHIVE"
    fi


    # ============================================================
    # Verify archive
    # ============================================================

    echo "Validating InterProScan archive..."

    EXPECTED_MD5=\$(awk '{print \$1}' "\$MD5_FILE")

    DETECTED_MD5=\$(md5sum "\$ARCHIVE" | awk '{print \$1}')

    if [ "\$EXPECTED_MD5" != "\$DETECTED_MD5" ]; then
        echo "ERROR: InterProScan archive checksum mismatch." >&2
        echo "Expected: \$EXPECTED_MD5" >&2
        echo "Detected: \$DETECTED_MD5" >&2
        exit 1
    fi

    echo "InterProScan archive checksum OK"


    # ============================================================
    # Extract data
    # ============================================================

    echo "Extracting InterProScan database..."

    rm -rf extracted
    mkdir -p extracted

    tar -pxzf "\$ARCHIVE" \
        -C extracted


    DATA_DIR=\$(find extracted \
        -type d \
        -path "*/interproscan-\${VERSION}/data" \
        -print \
        -quit)


    if [ -z "\$DATA_DIR" ]; then
        echo "ERROR: InterProScan data directory was not found after extraction." >&2
        exit 1
    fi


    if ! find "\$DATA_DIR" -type f -print -quit | grep -q .; then
        echo "ERROR: InterProScan data directory is empty." >&2
        exit 1
    fi


    mv "\$DATA_DIR" interproscan_data


    echo "========================================"
    echo "InterProScan database prepared"
    echo "Version: \$VERSION"
    echo "========================================"
    """
}

process INTERPROSCAN {
    label 'interproscan'
    tag "InterProScan annotation for $sample_id"
    publishDir "${params.outdir}/eukaryote_output_tools/interproscan", mode: 'copy'

    input:
    // Input comes from GFFREAD: tuple(val(sample_id), path("${sample_id}.faa")).
    tuple val(sample_id), path(faa)
    path ips_db

    output:
    // Standardize output names to a stable module prefix.
    tuple val(sample_id),
          path ("${sample_id}.interpro.*"),
          emit: results

    script:
    def outbase = "${sample_id}.interpro"
    def fmt = (params.ips_formats ?: 'tsv,gff3')

    """
    set -euo pipefail
    # Keep temporary files scoped to the task directory.

    mkdir -p temp

    /opt/interproscan/interproscan.sh \
      -i ${faa} \
      -f ${fmt} \
      -cpu ${task.cpus} \
      -goterms \
      --iprlookup \
      --pathways \
      -b ${outbase} \
      --tempdir temp
    """
}
