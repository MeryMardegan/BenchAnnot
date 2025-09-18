#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { GFFREAD     } from './modules/gffread.nf'
include { KOFAMSCAN   } from './modules/kofamscan.nf'
include { INTERPROSCAN } from './modules/interproscan.nf'


Channel
  .fromPath('$projectDir/data/eukaryotes/*.fna')
  .map { fa ->
    def id  = fa.baseName
    def ann = file("data/eukaryotes/${id}.gff")   // use GFF
    tuple(id, fa, ann)
  }
  .filter { id, fa, ann -> ann.exists() }
  .set { genome_pairs }

// DBs específicos do KOFAM
profiles = Channel.fromPath('$projectDir/data/eukaryotes/db/kofamscam/profiles/*', checkIfExists: true)
ko_list = Channel.fromPath('$projectDir/data/eukaryotes/db/kofamscam/ko_list', checkIfExists: true)


workflow gff_interpro {

    proteins = GFFREAD(genome_pairs)
    INTERPROSCAN(proteins)
}

workflow gff_kofamscan {

    proteins = GFFREAD(genome_pairs)
    KOFAMSCAN(proteins, profiles, ko_list)
}

workflow.onComplete { println "Workflow completed successfully!" }

