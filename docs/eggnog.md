# eggNOG-mapper

## Purpose

eggNOG-mapper assigns orthology-based functional annotations. BenchAnnot uses
eggNOG-mapper 2.1.13 in both annotation branches through `modules/eggnog.nf`.

## Processes and container

| Process | Label | Container |
| --- | --- | --- |
| `PREPARE_EGGNOG` | `eggnog_prepare` | `brunoholiva/eggnog_v2:2.1.13-sse41` |
| `EGGNOG_PROKARYOTE` | `eggnog_mapper_v2` | `brunoholiva/eggnog_v2:2.1.13-sse41` |
| `EGGNOG_EUKARYOTE` | `eggnog_mapper_v2` | `brunoholiva/eggnog_v2:2.1.13-sse41` |

The preparation process patches only the obsolete download host embedded in
the pinned eggNOG-mapper release. It then verifies the SQLite, DIAMOND, and
MMseqs2 database files before packaging.

## Database

`params.eggnog_db` defaults to:

```text
data/database/eggnog/eggnog_2026-09.sqsh
```

If that image is missing or fails basic file validation, `RESOLVE_EGGNOG_DB`
runs database preparation with the downloader's `-M` option and then runs
`PACK_SQUASHFS`. Existing images are checked only for the `.sqsh` extension,
existence, and non-zero size; the generated `.size` manifest is not consulted.

The `eggnog_2026-09` name is hardcoded, but the remote database content is not
pinned by release or checksum. The container version is fixed; the downloaded
database snapshot is not. The image is mounted read-only at `/database`; the
annotation process exports `EGGNOG_DATA_DIR=/database` and passes
`--data_dir /database`.

## Inputs

The prokaryotic process consumes:

```text
tuple(sample_id, genome_fasta), eggnog_database_image
```

It runs genome mode with Prodigal gene prediction, requests decorated GFF
output, and enables `--dbmem`. The eukaryotic process consumes the protein tuple
emitted by GFFread:

```text
tuple(sample_id, proteins_faa), eggnog_database_image
```

It runs protein mode and explicitly assigns `./tmp` with `--temp_dir`. Both
branches use the MMseqs2 search mode.

## Outputs

Outputs are flattened when published so the generated files are directly below
the tool directory:

```text
<outdir>/prokaryote_output_tools/eggnog/<sample_id>.emapper.*
<outdir>/eukaryote_output_tools/eggnog/<sample_id>_eggnog.emapper.*
```

The eukaryotic annotations filename matches the downstream contract:

```text
<sample_id>_eggnog.emapper.annotations
```

## Resources and retry policy

Database preparation requests 2 CPUs, 128 GB, and 72 hours. Annotation inherits
the default CPU count and requests 128 GB for 24 hours. Annotation retries once
to tolerate transient container or cache failures.

## Example

```bash
nextflow run main.nf --annotation_type eukaryote
```
