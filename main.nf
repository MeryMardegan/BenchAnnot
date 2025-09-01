#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PROKKA } from './modules/prokka.nf'
include { BAKTA } from './modules/bakta.nf'

workflow {
    fasta_ch = Channel.fromPath('data/*.fna')
    bakta_db_ch = Channel.value(params.bakta_db_dir)

    if (params.bakta_db_dir == null) {
        error "Error: Missing required parameter. Use --bakta_db_dir absolute/path/to/bakta/db on the command line."
    }
    
    PROKKA(fasta_ch)
    BAKTA(fasta_ch, bakta_db_ch)
}

