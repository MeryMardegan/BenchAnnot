# Functional Annotation Pipeline

This pipeline is designed to test and compare multiple functional annotation tools for genome annotation in a reproducible and scalable way.


## Purpose

- Orchestrate annotation tools and produce reproducible results for benchmarking and downstream analysis.

## Inputs
Place your genome files in the `data/` folder.
Run the pipeline with:
    ```bash
    nextflow run main.nf
The pipeline will automatically execute all modules defined in the workflow.

- Put FASTA files (`.fna`) in the `data/` directory.

## Outputs

- Annotated results for each genome are stored in the results/ directory.
- Intermediate files are stored in the work/ directory.

## Requirements

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)
- Input genome files in FASTA format (`.fna`) placed in the `data/` directory

## Modules

- `modules/prokka.nf` — Prokka 1.14.6 ([GitHub](https://github.com/tseemann/prokka), [Docker Image](https://hub.docker.com/r/staphb/prokka))
- `modules/bakta.nf` — Bakta 1.11.3 ([GitHub](https://github.com/oschwengers/bakta), [Docker Image](https://hub.docker.com/r/oschwengers/bakta))
- Prokka (modules/prokka.nf)
    Performs rapid genome annotation and outputs annotated GFF, GBK, and protein FASTA files.
    Results: results/prokka/<sample>_prokka/

- GFFREAD (modules/gffread.nf)
    Processes and manipulates GFF/GTF files (e.g., extracting transcript sequences).
    Results: results/gffread/<sample>.fa

- KofamScan (modules/kofamscan.nf)
    Assigns KEGG Orthologs (KOs) to protein sequences using HMM profiles. Requires KEGG database download.
    Results: results/kofamscan/<sample>_ko.txt


## Parameters

Bakta requires a specific database (in this case, version 6), that can be downloaded [here](https://zenodo.org/records/14916843). After downloading and extracting, provide the path in the run command, as follows `--bakta_db_dir absolute/path/to/bakta/db`.

## Run

- `nextflow run main.nf --bakta_db_dir path/to/bakta/db`

# Notes
- Pipeline ensures reproducibility and scalability using Docker containers.
- Additional tools for functional annotation (e.g., InterProScan, eggNOG-mapper, Funannotate) will be integrated in future updates.
- All intermediate files remain in work/, while final structured results are under results/.