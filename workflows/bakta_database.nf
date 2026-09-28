#!/usr/bin/env nextflow

include { PREPARE_BAKTA } from '../modules/bakta'
include { PACK_SQUASHFS } from '../modules/squashfs'
include {
    databaseImageIsValid
    managedDatabaseImageIsValid
} from '../lib/validation'


workflow RESOLVE_BAKTA_DB {

    main:

    // Prefer an external image, then fall back to the BenchAnnot-managed image.
    def expectedDb = params.bakta_db \
        ? params.bakta_db \
        : "${projectDir}/data/database/bakta/bakta_${params.bakta_db_version}.sqsh"

    // Managed images must also match the recorded size manifest.
    def dbIsValid = params.bakta_db \
        ? databaseImageIsValid(expectedDb) \
        : managedDatabaseImageIsValid(expectedDb)


    if (dbIsValid) {

        log.info "Using existing Bakta database: ${expectedDb}"

        bakta_db_ch = channel.value(file(expectedDb))

    } else {

        log.info "Bakta database image not found or incomplete. Preparing Bakta DB ${params.bakta_db_version}..."

        def cacheDir = file("${projectDir}/data/database/bakta/download")

        cacheDir.mkdirs()

        prepared = PREPARE_BAKTA(cacheDir.toString())

        packed = PACK_SQUASHFS(
            prepared.database_dir,
            'bakta',
            "bakta_${params.bakta_db_version}"
        )

        // BAKTA consumes only the image path from the packaging tuple.
        bakta_db_ch = packed.database.map { _database_name, database_image -> database_image }
    }


    emit:

    database = bakta_db_ch
}
