process PACK_SQUASHFS {

    label 'squashfs_tools'
    tag "${database_name}"

    publishDir { "data/database/${database_name}" }, mode: 'copy'

    input:
    path(database_dir)
    val (database_name) 
    val(image_name)

    output:
    tuple val(database_name), path("${image_name}.sqsh"), emit: database

    script:
    """
    set -euo pipefail

    mksquashfs \
        ${database_dir} \
        ${image_name}.sqsh \
        -noappend
    """
}