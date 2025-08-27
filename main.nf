#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PROKKA } from './modules/prokka.nf'

workflow {
    fasta_ch = Channel.fromPath('data/*.fna')

    PROKKA(fasta_ch)
}

