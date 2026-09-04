#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PROKARYOTE_INPUTS; EUKARYOTE_INPUTS } from './workflows/input.nf'
include { RESOLVE_EGGNOG_DB } from './workflows/eggnog_database.nf'
include { RESOLVE_BAKTA_DB } from './workflows/bakta_database.nf'
include { RESOLVE_PGAP } from './workflows/pgap_database.nf'
include { PROKARYOTE_ANNOTATION } from './workflows/prokaryote.nf'
include { EUKARYOTE_ANNOTATION }  from './workflows/eukaryote.nf'

workflow {

    main:

    // =====================================================================
    // Validate workflow selection
    // =====================================================================

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


    // =====================================================================
    // Shared eggNOG database
    //
    // Both annotation workflows use the same eggNOG database.
    // Therefore, it is resolved/prepared only once per execution.
    // =====================================================================
    
    RESOLVE_EGGNOG_DB()

    // =====================================================================
    // Prokaryotic annotation workflow
    // =====================================================================

    if (params.annotation_type in ['prokaryote', 'both']) {

        PROKARYOTE_INPUTS()
        RESOLVE_BAKTA_DB()
        RESOLVE_PGAP()

        PROKARYOTE_ANNOTATION(
            PROKARYOTE_INPUTS.out.samples,
            RESOLVE_BAKTA_DB.out.database,
            RESOLVE_EGGNOG_DB.out.database,
            RESOLVE_PGAP.out.pgap_python,
            RESOLVE_PGAP.out.pgap_dir,
            RESOLVE_PGAP.out.pgap_container
        )
    }


    // =====================================================================
    // Eukaryotic annotation workflow
    // =====================================================================

    if (params.annotation_type in ['eukaryote', 'both']) {

        EUKARYOTE_INPUTS()

        EUKARYOTE_ANNOTATION(
            EUKARYOTE_INPUTS.out.samples,
            EUKARYOTE_INPUTS.out.kofam_db,
            EUKARYOTE_INPUTS.out.ips_db,
            RESOLVE_EGGNOG_DB.out.database
        )
    }

    onComplete:
    if (workflow.success) {
        println "Workflow completed successfully!"
    } else {
        println "Workflow failed."
    }
}
