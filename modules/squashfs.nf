process PACK_SQUASHFS {

    label 'squashfs_tools'
    tag "${database_name}"

    publishDir { "data/database/${database_name}" }, mode: 'copy', overwrite: true

    input:
    path(database_dir)
    val (database_name) 
    val (image_name)

    output:
    tuple val(database_name), path("${image_name}.sqsh"), emit: database

    path "${image_name}.sqsh.size", emit: size_manifest

    script:
    """
    set -euo pipefail

    REAL_DB=\$(readlink -f "${database_dir}")
    
    if [ ! -d "\$REAL_DB" ]; then
        echo "ERROR: database directory not found: \$REAL_DB" >&2
        exit 1
    fi

    mksquashfs \
        "\$REAL_DB" \
        "${image_name}.sqsh" \
        -noappend \
        -processors ${task.cpus}

   stat -c '%s' \
        "${image_name}.sqsh" \
        > "${image_name}.sqsh.size"
 
    """
}
