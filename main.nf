#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PROKARYOTE_INPUTS; EUKARYOTE_INPUTS } from './workflows/input.nf'
include { RESOLVE_EGGNOG_DB } from './workflows/eggnog_database.nf'
include { RESOLVE_BAKTA_DB } from './workflows/bakta_database.nf'
include { RESOLVE_PGAP } from './workflows/pgap_database.nf'
include { RESOLVE_KOFAM_DB } from './workflows/kofam_database.nf'
include { PROKARYOTE_ANNOTATION } from './workflows/prokaryote.nf'
include { EUKARYOTE_ANNOTATION } from './workflows/eukaryote.nf'
include { RESOLVE_INTERPROSCAN_DB } from './workflows/interproscan_database.nf'

workflow {

    main:

    // ============================================================================
    // Pipeline orchestration
    // ============================================================================

    // ---------------------------------------------------------------------------
    // Workflow selection
    // ---------------------------------------------------------------------------

    def validTypes = ['prokaryote', 'eukaryote', 'both']

    if (!(params.annotation_type in validTypes)) {
        error """
        Invalid --annotation_type: ${params.annotation_type}

        Valid options:
          --annotation_type prokaryote
          --annotation_type eukaryote
          --annotation_type both
        """
    }
    // Resolve this shared database once for both annotation branches.
    RESOLVE_EGGNOG_DB()

    // ---------------------------------------------------------------------------
    // Prokaryotic annotation
    // ---------------------------------------------------------------------------

    if (params.annotation_type in ['prokaryote', 'both']) {

        PROKARYOTE_INPUTS()
        RESOLVE_BAKTA_DB()
        RESOLVE_PGAP()

        PROKARYOTE_ANNOTATION(
            PROKARYOTE_INPUTS.out.samples,
            RESOLVE_BAKTA_DB.out.database,
            RESOLVE_EGGNOG_DB.out.database,
            RESOLVE_PGAP.out.pgap_dir,
            RESOLVE_PGAP.out.pgap_container,
            RESOLVE_PGAP.out.pgap_python
        )
    }
    // ---------------------------------------------------------------------------
    // Eukaryotic annotation
    // ---------------------------------------------------------------------------

    if (params.annotation_type in ['eukaryote', 'both']) {

        EUKARYOTE_INPUTS()
        RESOLVE_KOFAM_DB()
        RESOLVE_INTERPROSCAN_DB()

        EUKARYOTE_ANNOTATION(
            EUKARYOTE_INPUTS.out.samples,
            RESOLVE_KOFAM_DB.out.database,
            RESOLVE_INTERPROSCAN_DB.out.database,
            RESOLVE_EGGNOG_DB.out.database
        )
    }

    // ---------------------------------------------------------------------------
    // Completion status
    // ---------------------------------------------------------------------------

    onComplete:
    if (workflow.success) {
        println "Workflow completed successfully!"
    } else {
        println "Workflow failed."
    }
}
