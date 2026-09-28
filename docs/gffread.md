# GFFread

## Purpose

GFFread extracts eukaryotic protein sequences from the reference genome and GFF
before functional annotation. BenchAnnot uses GFFread 0.12.7.

## Process and container

| Process | Label | Container |
| --- | --- | --- |
| `GFFREAD` | `gffread` | `quay.io/biocontainers/gffread:0.12.7--h9a82719_0` |

## Input

```text
tuple(sample_id, genome_fasta, reference_gff)
```

The workflow validates `reference_faa` and `organism_id` from the samplesheet,
but the current annotation branch does not pass them to GFFread.

## Filtering

Before sequence extraction, the process preserves GFF headers and removes:

- records with undefined strand `?`;
- records whose attributes contain `exception=trans-splicing`.

The process reports the total, retained, and removed record counts and fails if
no feature remains. This filtered GFF is a pipeline artifact; the downstream
audit reconstructs its own filtered reference.

Sequence extraction runs `gffread -F -S -C -J`. These options are part of the
current biological extraction contract and must not be changed as formatting
cleanup.

## Outputs

```text
<outdir>/eukaryote_output_tools/gffread/<sample_id>.filtered.gff
<outdir>/eukaryote_output_tools/gffread/<sample_id>_gffread.faa
```

The protein tuple is shared by KofamScan, InterProScan, and eggNOG-mapper.

## Resources

GFFread inherits the default 10 CPUs, 128 GB, and 24-hour allocation, although
the command itself does not receive an explicit thread argument.

## Example

```bash
nextflow run main.nf --annotation_type eukaryote
```
