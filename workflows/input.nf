#!/usr/bin/env nextflow

include {
    validateProkaryoteRow;
    validateEukaryoteRow;
    validateBaktaInputs;
    validateKofamInputs;
    validateInterproscanInputs
} from '../lib/validation.nf'


workflow PROKARYOTE_INPUTS {

    main:
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


    emit:
    samples        = samples_ch
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