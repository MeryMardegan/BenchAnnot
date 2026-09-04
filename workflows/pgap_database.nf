#!/usr/bin/env nextflow

include { PREPARE_PGAP_PYTHON; PREPARE_PGAP } from '../modules/pgap'


workflow RESOLVE_PGAP {

    main:

    /*
     * Resolve the Python runtime required by the PGAP wrapper.
     */

    def expectedPython =
        file("${params.pgap_python_dir}/bin/python")

    if (
        expectedPython.exists() &&
        expectedPython.isFile()
    ) {

        log.info "Using existing PGAP Python runtime: ${params.pgap_python_dir}"

        pgap_python_ch =
            channel.value(params.pgap_python_dir)

    } else {

        log.info "PGAP Python runtime not found. Preparing Python ${params.pgap_python_version}..."

        PREPARE_PGAP_PYTHON(
            params.pgap_python_dir
        )

        pgap_python_ch =
            PREPARE_PGAP_PYTHON.out.runtime
    }


    /*
     * Resolve PGAP itself.
     */

    def expectedDir =
        file(params.pgap_dir)

    def expectedContainer =
        file(params.pgap_container)

    if (
        expectedDir.exists() &&
        expectedDir.isDirectory() &&
        expectedContainer.exists() &&
        expectedContainer.isFile() &&
        expectedContainer.size() > 0
    ) {

        log.info "Using existing PGAP installation: ${params.pgap_dir}"

        pgap_dir_ch =
            channel.value(params.pgap_dir)

        pgap_container_ch =
            channel.value(params.pgap_container)

    } else {

        log.info "PGAP installation not found. Preparing PGAP ${params.pgap_version}..."

        prepared = PREPARE_PGAP(
            pgap_python_ch
        )

        pgap_dir_ch =
            prepared.runtime.map { runtime ->
                runtime.toString()
            }

        pgap_container_ch =
            prepared.runtime.map { runtime ->
                file("${runtime}/pgap_${params.pgap_version}.sif")
            }
    }


    emit:

    pgap_dir =
        pgap_dir_ch

    pgap_container =
        pgap_container_ch

    pgap_python =
        pgap_python_ch
}