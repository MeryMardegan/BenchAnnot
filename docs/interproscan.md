# InterProScan

## Purpose

InterProScan combines protein family, domain, site, GO, and pathway signatures.
BenchAnnot uses InterProScan 5.75-106.0 for eukaryotic protein annotation.

## Processes and containers

| Process | Label | Container |
| --- | --- | --- |
| `PREPARE_INTERPROSCAN` | `interproscan_prepare` | `docker://python:3.11-bookworm` |
| `INTERPROSCAN` | `interproscan` | `interpro/interproscan:5.75-106.0` |

The Python preparation container downloads, verifies, and extracts the data
archive for the configured InterProScan release. The annotation container
provides the InterProScan executable.

## Database

`params.ips_db` may point to an existing non-empty `.sqsh` image. When it is
unset or invalid, `RESOLVE_INTERPROSCAN_DB` prepares the pinned release and
packages it as:

```text
data/database/interproscan/interproscan_5.75-106.0.sqsh
```

The SquashFS root is mounted read-only at `/opt/interproscan/data`, where the
container expects its data directory.

`interproscan_version` controls data download and image naming, but it does not
change the hardcoded `interpro/interproscan:5.75-106.0` container. Software and
data releases must be updated together in code.

## Input

`INTERPROSCAN` receives the protein tuple emitted by GFFread plus the database
image:

```text
tuple(sample_id, <sample_id>_gffread.faa), interproscan_database_image
```

## Outputs

`params.ips_formats` controls the requested formats and defaults to `tsv,gff3`.
The stable output prefix is `<sample_id>.interpro`:

```text
<outdir>/eukaryote_output_tools/interproscan/<sample_id>.interpro.tsv
<outdir>/eukaryote_output_tools/interproscan/<sample_id>.interpro.gff3
```

The command also requests GO terms, InterPro lookups, and pathway mappings.

## Resources

Preparation requests 1 CPU, 8 GB, and 48 hours. Annotation requests 10 CPUs,
128 GB, and 24 hours. Temporary files remain inside the Nextflow task directory.

## Example

```bash
nextflow run main.nf --annotation_type eukaryote --ips_formats tsv,gff3
```
