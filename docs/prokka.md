# Prokka

## Purpose

Prokka performs structural and functional annotation of prokaryotic genomes.
BenchAnnot runs Prokka 1.14.6 through `modules/prokka.nf`.

## Process and container

| Process | Label | Container |
| --- | --- | --- |
| `PROKKA` | `prokka` | `https://depot.galaxyproject.org/singularity/prokka%3A1.14.6--pl5321hdfd78af_5` |

## Input

```text
tuple(sample_id, genome_fasta)
```

The process uses `sample_id` as the output prefix, enables compliant output,
and passes the allocated CPU count to Prokka. Species, taxid, and genetic code
metadata are not passed to this process.

## Outputs

Prokka writes its standard output set and BenchAnnot flattens those files into:

```text
<outdir>/prokaryote_output_tools/prokka/<sample_id>.*
```

This includes the GFF and GenBank files consumed by the prokaryotic analysis.
The standard Prokka set also includes protein, nucleotide, feature-table,
submission, TSV, text, and log files matching `results/<sample_id>.*`.

## Resources

The process inherits the default 10 CPUs and requests 64 GB for 24 hours. A
task-local temporary directory is assigned through `TMPDIR`.

## Example

```bash
nextflow run main.nf --annotation_type prokaryote
```
