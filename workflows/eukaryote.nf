include { GFFREAD } from '../modules/gffread'
include { KOFAMSCAN } from '../modules/KofamScan'
include { INTERPROSCAN } from '../modules/interproscan'
include { EGGNOG_EUKARYOTE } from '../modules/eggnog'

workflow EUKARYOTE_ANNOTATION {
    
    take:

    samples_ch
    kofam_db
    ips_db_ch
    eggnog_db_ch

    main:    

    samples_ch
        .map { sample_id, genome_fasta, reference_gff, _reference_faa, _organism_id ->
            tuple(sample_id, genome_fasta, reference_gff)
        }
        .set { gffread_ch }

    GFFREAD(gffread_ch)
    proteins_ch = GFFREAD.out.proteins
    KOFAMSCAN(proteins_ch, kofam_db)
    INTERPROSCAN(proteins_ch, ips_db_ch)
    EGGNOG_EUKARYOTE(proteins_ch, eggnog_db_ch)
}