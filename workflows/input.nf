#!/usr/bin/env nextflow

include {
    validateProkaryoteRow;
    validateEukaryoteRow;
    validateBaktaInputs;
    validatePgapInputs;
    validateKofamInputs;
    validateInterproscanInputs
} from '../lib/validation.nf'


workflow PROKARYOTE_INPUTS {

    main:

    /*
     * Validate external resources required by the prokaryotic workflow.
     * eggNOG is resolved separately by RESOLVE_EGGNOG_DB.
     */
    validateBaktaInputs()
    validatePgapInputs()

    /*
     * validateProkaryoteRow() returns:
     * tuple(sample_id, genome_fasta, species, taxid, genetic_code)
     */
    samples_ch = channel
        .fromPath(
            params.prokaryote_samplesheet,
            checkIfExists: true
        )
        .splitCsv(header: true)
        .map { row ->
            validateProkaryoteRow(row)
        }

    /*
     * Bakta DB is now a SquashFS image.
     */
    bakta_db_ch = channel.value(
        file(params.bakta_db)
    )

    /*
     * PGAP remains an external installation directory
     * plus its container image.
     */
    pgap_dir_ch = channel.value(
        file(params.pgap_dir)
    )

    pgap_container_ch = channel.value(
        params.pgap_container
    )

    emit:
    samples        = samples_ch
    bakta_db       = bakta_db_ch
    pgap_dir       = pgap_dir_ch
    pgap_container = pgap_container_ch
}


workflow EUKARYOTE_INPUTS {

    main:

    /*
     * Validate external resources required by the eukaryotic workflow.
     * eggNOG is resolved separately by RESOLVE_EGGNOG_DB.
     */
    validateKofamInputs()
    validateInterproscanInputs()

    /*
     * validateEukaryoteRow() returns:
     * tuple(
     *     sample_id,
     *     genome_fasta,
     *     reference_gff,
     *     reference_faa,
     *     organism_id
     * )
     */
    samples_ch = channel
        .fromPath(
            params.eukaryote_samplesheet,
            checkIfExists: true
        )
        .splitCsv(header: true)
        .map { row ->
            validateEukaryoteRow(row)
        }

    /*
     * Databases are provided as SquashFS images.
     */
    kofam_db_ch = channel.value(
        file(params.kofamscan_db)
    )

    ips_db_ch = channel.value(
        file(params.ips_db)
    )

    emit:
    samples  = samples_ch
    kofam_db = kofam_db_ch
    ips_db   = ips_db_ch
}