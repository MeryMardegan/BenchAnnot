#!/usr/bin/env nextflow

include { PREPARE_INTERPROSCAN } from '../modules/interproscan'
include { PACK_SQUASHFS } from '../modules/squashfs'

include {
    databaseImageIsValid;
    managedDatabaseImageIsValid
} from '../lib/validation'


workflow RESOLVE_INTERPROSCAN_DB {

    main:

    def expectedDb = params.ips_db \
        ? params.ips_db \
        : "${projectDir}/data/database/interproscan/interproscan_${params.interproscan_version}.sqsh"


    def dbIsValid = params.ips_db \
        ? databaseImageIsValid(expectedDb) \
        : managedDatabaseImageIsValid(expectedDb)


    if (dbIsValid) {

        log.info "Using existing InterProScan database: ${expectedDb}"

        ips_db_ch = channel.value(
            file(expectedDb)
        )

    } else {

        log.info "InterProScan database image not found or incomplete. Preparing ${params.interproscan_version}..."

        def cacheDir = file(
            "${projectDir}/data/database/interproscan/download"
        )

        cacheDir.mkdirs()

        prepared = PREPARE_INTERPROSCAN(
            cacheDir.toString()
        )

        packed = PACK_SQUASHFS(
            prepared.database_dir,
            'interproscan',
            "interproscan_${params.interproscan_version}"
        )

        ips_db_ch = packed.database.map {
            database_name, database_image ->
                database_image
        }
    }


    emit:

    database = ips_db_ch
}
