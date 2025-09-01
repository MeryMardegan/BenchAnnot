# BenchAnnot

A Nextflow pipeline to run genome annotation on FASTA assemblies.

## Purpose

- Orchestrate annotation tools and produce reproducible results for benchmarking and downstream analysis.

## Inputs

- Put FASTA files (`.fna`) in the `data/` directory.

## Outputs

- Per-sample annotation directories under `results/module/sample>_module/`, where `module` is the tool used and `sample` is the genome used.

## Requirements

- Nextflow
- Docker
- Specific databases

## Modules

- `modules/prokka.nf` — Prokka 1.14.6 ([GitHub](https://github.com/tseemann/prokka), [Docker Image](https://hub.docker.com/r/staphb/prokka))
- `modules/bakta.nf` — Bakta 1.11.3 ([GitHub](https://github.com/oschwengers/bakta), [Docker Image](https://hub.docker.com/r/oschwengers/bakta))

## Parameters

Bakta requires a specific database (in this case, version 6), that can be downloaded [here](https://zenodo.org/records/14916843). After downloading and extracting, provide the path in the run command, as follows `--bakta_db_dir absolute/path/to/bakta/db`.

## Run

- `nextflow run main.nf --bakta_db_dir path/to/bakta/db`