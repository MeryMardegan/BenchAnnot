#!/usr/bin/env nextflow

include { PREPARE_BAKTA } from '../modules/bakta'
include { PACK_SQUASHFS } from '../modules/squashfs'
include { databaseImageIsValid } from '../lib/validation'


workflow RESOLVE_BAKTA_DB {

    main:

    if (databaseImageIsValid(params.bakta_db)) {

        log.info "Using existing Bakta database: ${params.bakta_db}"

        bakta_db_ch = channel.value(
            file(params.bakta_db)
        )

    } else {

        log.info "Bakta database image not found. Preparing database..."

        prepared = PREPARE_BAKTA()

        packed = PACK_SQUASHFS(
            prepared.database_dir,
            'bakta',
            'bakta_2026-09'
        )

        bakta_db_ch = packed.database
    }

    emit:
    database = bakta_db_ch
}