# Bakta

## Purpose

Bakta annotates bacterial genomes using a versioned reference database.
BenchAnnot uses Bakta 1.11.3 with database schema version 6.0 and AMRFinderPlus
series 4.0 release `2025-07-16.1`.

## Processes and container

| Process | Label | Container |
| --- | --- | --- |
| `PREPARE_BAKTA` | `bakta_prepare` | `oschwengers/bakta:v1.11.3` |
| `BAKTA` | `bakta` | `oschwengers/bakta:v1.11.3` |

## Database

An external database can be supplied with `--bakta_db /path/to/database.sqsh`.
When no image is supplied, `RESOLVE_BAKTA_DB` prepares the configured Bakta and
AMRFinderPlus releases and packages:

```text
data/database/bakta/bakta_6.0.sqsh
```

Managed images have a `.size` manifest and are reused only when their current
size matches that manifest. External images receive basic file and extension
validation. The selected image is mounted read-only at `/database`.

The Bakta version parameter is coupled to the fixed Zenodo record, configured
MD5, AMRFinder release, and Bakta 1.11.3 container. Changing
`bakta_db_version` alone is not a safe database upgrade.

## Input

```text
tuple(sample_id, genome_fasta), bakta_database_image
```

## Outputs

Bakta runs in compliant mode, skips plot generation, and publishes its standard
files directly below:

```text
<outdir>/prokaryote_output_tools/bakta/<sample_id>.*
```

The standard set includes GFF3, GenBank/GBFF, protein and nucleotide FASTA,
TSV, JSON, EMBL, text, and other files matching `results/<sample_id>.*`; plots
are disabled.

## Resources

Database preparation requests 1 CPU, 128 GB, and 72 hours. Annotation inherits
the default CPU count and requests 128 GB for 24 hours. Temporary and Matplotlib
cache files remain inside the task directory.

## Example

```bash
nextflow run main.nf --annotation_type prokaryote
```
