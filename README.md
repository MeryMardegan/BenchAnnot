# BenchAnnot

## Requirements

- Nextflow
- Docker
- Input genome files in FASTA format (`.fna`) placed in the `data/` directory


## Output

- Annotated results for each genome will be found in `results/prokka/<sample>_prokka/`, where `<sample>` is the base name of your input file.

## Module

- The Prokka process is defined in `modules/prokka.nf` and automatically handles each input genome.