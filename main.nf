#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PROKKA } from './modules/prokka.nf'
include { BAKTA } from './modules/bakta.nf'
include { GFFREAD     } from './modules/gffread.nf'
include { KOFAMSCAN   } from './modules/kofamscan.nf'


Channel
  .fromPath('data/eukaryotes/*.fna')
  .map { fa ->
    def id  = fa.baseName
    def ann = file("data/eukaryotes/${id}.gff")   // use GFF
    tuple(id, fa, ann)
  }
  .filter { id, fa, ann -> ann.exists() }
  .set { genome_pairs }

// DB já baixado manualmente
profiles = Channel.fromPath('data/eukaryotes/db/profiles', checkIfExists: true)
ko_list = Channel.fromPath('data/eukaryotes/db/ko_list', checkIfExists: true)

workflow {
    fasta_ch = Channel.fromPath('data/*.fna')
    bakta_db_ch = Channel.value(params.bakta_db_dir)

    if (params.bakta_db_dir == null) {
        error "Error: Missing required parameter. Use --bakta_db_dir absolute/path/to/bakta/db on the command line."
    }
    
    PROKKA(fasta_ch)
    BAKTA(fasta_ch, bakta_db_ch)

    proteins = GFFREAD(genome_pairs)
    KOFAMSCAN(proteins, profiles, ko_list)
}

workflow.onComplete { println "Workflow completed successfully!" }

