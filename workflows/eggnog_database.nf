#!/usr/bin/env nextflow

include { PREPARE_EGGNOG } from '../modules/eggnog'
include { PACK_SQUASHFS }  from '../modules/squashfs'
include { databaseImageIsValid } from '../lib/validation'

workflow RESOLVE_EGGNOG_DB {

    main:

    if (databaseImageIsValid(params.eggnog_db)) {

        log.info "Using existing eggNOG database: ${params.eggnog_db}"
        eggnog_db_ch = channel.value(file(params.eggnog_db))

    } else {

        log.info "eggNOG database image not found. Preparing database..."

        prepared = PREPARE_EGGNOG()
        packed = PACK_SQUASHFS(prepared.database_dir, 'eggnog', 'eggnog_2026-09')
        eggnog_db_ch = packed.database.map { _database_name, database_image -> database_image }
    }

    emit:
    database = eggnog_db_ch
}
