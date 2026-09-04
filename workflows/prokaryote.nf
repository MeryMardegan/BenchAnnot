include { PROKKA }            from '../modules/prokka'
include { BAKTA }             from '../modules/bakta'
include { PGAP }              from '../modules/pgap'
include { EGGNOG_PROKARYOTE } from '../modules/eggnog'


workflow PROKARYOTE_ANNOTATION {

    take:
    samples_ch
    bakta_db_ch
    eggnog_db_ch
    pgap_dir_ch
    pgap_container_ch
    pgap_python_ch

    main:

    samples_ch.multiMap { sample_id, genome_fasta, species, taxid, genetic_code ->

        prokka: tuple(
            sample_id,
            genome_fasta
        )

        bakta: tuple(
            sample_id,
            genome_fasta,
            species,
            taxid,
            genetic_code
        )

        eggnog: tuple(
            sample_id,
            genome_fasta
        )

        pgap: tuple(
            sample_id,
            genome_fasta,
            species
        )
    }
    .set { inputs }


    PROKKA(
        inputs.prokka
    )

    BAKTA(
        inputs.bakta,
        bakta_db_ch
    )

    EGGNOG_PROKARYOTE(
        inputs.eggnog,
        eggnog_db_ch
    )

    PGAP(
        inputs.pgap,
        pgap_dir_ch,
        pgap_container_ch,
        pgap_python_ch
    )
}