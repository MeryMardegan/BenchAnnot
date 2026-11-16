# KofamScan

## Purpose

KofamScan assigns KEGG Ortholog identifiers using profile HMMs and
profile-specific score thresholds. BenchAnnot uses KofamScan 1.3.0 for the
eukaryotic branch.

## Processes and container

| Process | Label | Container |
| --- | --- | --- |
| `PREPARE_KOFAM` | `kofamscan_prepare` | host tools configured by the executor |
| `KOFAMSCAN` | `kofamscan` | `docker://quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2` |

Preparation downloads the profiles and KO list for the pinned archive release,
extracts them, and verifies that `profiles/` is non-empty and `ko_list` exists.
Preparation does not explicitly verify `profiles/eukaryote.hal`, which is the
profile used during annotation.

## Database

`params.kofamscan_db` may point to an existing non-empty `.sqsh` image. If no
valid image is supplied, the pipeline uses `params.kofamscan_db_release`
(`2026-07-02` by default), prepares the database, and packages:

```text
data/database/kofamscan/kofamscan_2026-07-02.sqsh
```

The image root contains `profiles/` and `ko_list`. It is mounted read-only at
`/database`; annotation uses `/database/profiles/eukaryote.hal` and
`/database/ko_list`.

## Input

`KOFAMSCAN` consumes the protein tuple emitted by GFFread plus the database
image:

```text
tuple(sample_id, <sample_id>_gffread.faa), kofamscan_database_image
```

## Output

KofamScan runs `exec_annotation` in `detail-tsv` mode and reports unannotated
queries:

```text
<outdir>/eukaryote_output_tools/kofamscan/<sample_id>.kofam.txt
```

The downstream analysis retains significant rows marked with `*`; that filter
is not applied by the Nextflow process.

## Resources

Preparation requests 1 CPU, 4 GB, and 24 hours. Annotation inherits the default
CPU count and requests 128 GB for 24 hours.

## Example

```bash
nextflow run main.nf --annotation_type eukaryote
```
