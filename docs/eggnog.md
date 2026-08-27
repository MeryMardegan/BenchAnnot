# eggNOG-mapper

## Overview
eggNOG-mapper assigns functional annotations to protein sequences using precomputed orthology and HMM profiles. In this pipeline, it consumes the protein FASTA produced by GFFREAD.

## Inputs
- Protein FASTA (`.faa`) produced by GFFREAD.

## Outputs
- `${sample_id}_eggnog.emapper.*` stored under:
	- `results/eukaryotes/eggnog/`

## Database Setup
Download the recommended genome assembly databases:
- `mmseqs.tar.gz`
- `eggnog.db.gz`

Extract them into a directory accessible by the pipeline, for example:
`data/eukaryotes/db/eggnog/`

## Configuration
Set the database path in `nextflow.config`:
- `params.emapper_data_dir = "${projectDir}/data/eukaryotes/db/eggnog"`

The container binds this directory to `/eggnog-data` and exports `EGGNOG_DATA_DIR`.

## Example Run
```bash
nextflow run main.nf -entry eukaryotes_annot
```

## Notes
- The module uses `emapper.py -m mmseqs` for speed on large proteomes.
- Ensure the database path is readable by the container runtime (Docker/Apptainer).