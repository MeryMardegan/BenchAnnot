# eggNOG-mapper

## Overview
eggNOG-mapper assigns functional annotations to protein sequences using precomputed orthology and HMM profiles. In this pipeline, it consumes the protein FASTA produced by GFFREAD.

## Inputs
- Protein FASTA (`.faa`) produced by GFFREAD.

## Outputs
- `${sample_id}_eggnog.emapper.*` stored under:
	- `data/reproduced/eukaryote_output_tools/eggnog/` by default

## Database Setup
Download the recommended genome assembly databases:
- `mmseqs.tar.gz`
- `eggnog.db.gz`

The pipeline packages the prepared files as:
`data/database/eggnog/eggnog_2026-09.sqsh`

## Configuration
Set the database image path in `nextflow.config`:
- `params.eggnog_db = "${projectDir}/data/database/eggnog/eggnog_2026-09.sqsh"`

The container mounts the image at `/database` and exports that location through
`EGGNOG_DATA_DIR`.

## Example Run
```bash
nextflow run main.nf --annotation_type eukaryote
```

## Notes
- The module uses `emapper.py -m mmseqs` for speed on large proteomes.
- Ensure the database path is readable by the container runtime (Docker/Apptainer).
