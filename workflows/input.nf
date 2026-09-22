#!/usr/bin/env nextflow

include {
    validateProkaryoteRow;
    validateEukaryoteRow;
    validateBaktaInputs;
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

    emit:
    samples  = samples_ch
}
