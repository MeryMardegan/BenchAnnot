#!/usr/bin/env nextflow

include { PREPARE_KOFAM } from '../modules/KofamScan'
include { PACK_SQUASHFS } from '../modules/squashfs'

include {
    databaseImageIsValid;
    managedDatabaseImageIsValid
} from '../lib/validation'


workflow RESOLVE_KOFAM_DB {

    main:

    def expectedDb = params.kofamscan_db \
        ? params.kofamscan_db \
        : "${projectDir}/data/database/kofamscan/kofamscan_${params.kofamscan_db_release}.sqsh"

    def dbIsValid = params.kofamscan_db \
        ? databaseImageIsValid(expectedDb) \
        : managedDatabaseImageIsValid(expectedDb)


    if (dbIsValid) {

        log.info "Using existing Kofam database: ${expectedDb}"

        kofam_db_ch = channel.value(
            file(expectedDb)
        )

    } else {

        log.info "Kofam database image not found or incomplete. Preparing release ${params.kofamscan_db_release}..."

        def cacheDir = file(
            "${projectDir}/data/database/kofamscan/download"
        )

        cacheDir.mkdirs()

        prepared = PREPARE_KOFAM(
            cacheDir.toString()
        )

        packed = PACK_SQUASHFS(
            prepared.database_dir,
            'kofamscan',
            "kofamscan_${params.kofamscan_db_release}"
        )

        kofam_db_ch = packed.database.map {
            _database_name, database_image ->
                database_image
        }
    }


    emit:

    database = kofam_db_ch
}
