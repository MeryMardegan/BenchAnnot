#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Module imports

include { GFFREAD     } from './modules/gffread.nf'
include { KOFAMSCAN   } from './modules/kofamScan.nf'
include { INTERPROSCAN } from './modules/interproscan.nf'
include { PROKKA } from './modules/prokka.nf'
include { BAKTA } from './modules/bakta.nf'
include { EGGNOG } from './modules/eggnog.nf'
include { PGAP } from './modules/pgap.nf'

//===========
//PROKARYOTES
//===========

// Validate required parameters for prokaryotic annotation
def validateProkaryoticInputs() {
    def required = ['bakta_db_dir', 'eggnog_db_dir', 'pgap_dir']
    def missing = required.findAll { !params[it] }

    if (missing) {
        error "Missing required parameters: ${missing.collect { "--${it}" }.join(', ')}"
    }

    required.each { name ->
        def location = file(params[name])
        if (!location.exists() || !location.isDirectory()) {
            error "Parameter --${name} must be an existing directory: ${params[name]}"
        }
    }
}

// Build all channels needed for prokaryotic annotation
def createProkaryoticChannels() {
    validateProkaryoticInputs()

    def fasta_ch = Channel.fromPath("$projectDir/data/genome_prokaryote/*.fna", checkIfExists: true)
    def bakta_db_ch = Channel.value(params.bakta_db_dir)
    def eggnog_db_ch = Channel.value(params.eggnog_db_dir)
    def pgap_dir_ch = Channel.value(params.pgap_dir)

    return [fasta_ch, bakta_db_ch, eggnog_db_ch, pgap_dir_ch]

}

//===========
//EUKARYOTES
//===========

// Pair each genome FASTA with its corresponding GFF annotation
def createEukaryoticGenomePairs() {
    Channel
    .fromPath("$projectDir/data/genome_eukaryote/*.fna", checkIfExists: true)
    .map { fa ->
        def id  = fa.baseName
        def ann = file("data/eukaryotes/${id}.gff")   // use GFF
        tuple(id, fa, ann)
      }
    .filter { id, fa, ann -> 
        if (!ann.exists()) {
            log.warn "Missing GFF annotation for genome: '${id}'. (${ann}) excluded."
            return false
        }
        return true
      } 
}

// KofamScan database eukaryote inputs
def createEukaryoticKofamChannels() {
    def profiles = Channel.fromPath("$projectDir/data/eukaryotes/db/kofamscam/profiles/*", checkIfExists: true)
    def ko_list = Channel.fromPath("$projectDir/data/eukaryotes/db/kofamscam/ko_list", checkIfExists: true)
    return [profiles, ko_list]
}

workflow prokaryote_annotation {
  // Create channels for prokaryotic annotation
  (fasta_ch, bakta_db_ch, eggnog_db_ch, pgap_dir_ch) = createProkaryoticChannels()
  // 1) Run Prokka on the genome FASTA files
  PROKKA(fasta_ch)
  // 2) Run BAKTA on the genome FASTA files
  BAKTA(fasta_ch, bakta_db_ch)
  // 3) Run eggNOG-mapper on the genome FASTA files
  EGGNOG(fasta_ch, eggnog_db_ch)
  // 4) Run PGAP on the genome FASTA files
  PGAP(fasta_ch, pgap_dir_ch)
}

workflow eukaryote_annotation {
  // Create channels for eukaryotic annotation
  genome_pairs = createEukaryoticGenomePairs()
  (profiles, ko_list) = createEukaryoticKofamChannels()
  // 1) Build protein sequences from genome + GFF using gffread
  proteins = GFFREAD(genome_pairs)
  // 2) Run InterProScan on the protein sequences
  INTERPROSCAN(proteins)
  // 3) Run KofamScan on the protein sequences
  KOFAMSCAN(proteins, profiles, ko_list)
  // 4) Run eggNOG-mapper on the protein sequences
  EGGNOG(proteins, eggnog_db_ch)
}

workflow.onComplete { println "Workflow completed successfully!" }