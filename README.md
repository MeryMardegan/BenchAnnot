# BenchAnnot

BenchAnnot is a reproducible framework for benchmarking genome annotation
tools against curated reference annotations. This repository has two parts:

- **Analysis workflow** — notebooks and Python modules that audit, clean, and
  compare annotation tool outputs (eukaryotic and prokaryotic).
- **Nextflow pipeline** — `main.nf` that orchestrates functional annotation
  tools (Prokka, Bakta, GFFread, KofamScan) in a reproducible and scalable way.

## Analysis workflow

The current refactor focuses on the eukaryotic post-processing workflow for
*Saccharomyces cerevisiae* and *Drosophila melanogaster*. Prokaryotic notebooks
remain available, but their planned improvements are documented separately.

Run the eukaryotic notebooks from `2_run/notebooks/eukaryotic/` in this order:

1. `1_audit_prepare_reference.ipynb` audits GFF/GFFread processing and builds
   an independent NCBI reference universe.
2. `2_prepare_output_tools.ipynb` cleans Kofam, Pannzer, EggNOG, and
   InterProScan outputs using exact `RNA_ID` values.
3. `3_functional_analysis.ipynb` compares all reference RNA records, performs
   the deterministic first-transcript sensitivity analysis, and evaluates the
   reviewed Swiss-Prot canonical subset.

The NCBI FAA is the functional reference universe. GFFread is used to audit
sequence-extraction losses; it does not define the reference denominator.

### Quick start

Create the locked environment and install the package:

```bash
conda-lock install --name benchannot 1_setup/conda-lock.yml
conda activate benchannot
python -m pip install -e .
python -m ipykernel install --user --name benchannot --display-name "Python (BenchAnnot)"
```

Then open the project in Jupyter and execute the three eukaryotic notebooks
top to bottom. Detailed setup and workflow notes are in [`docs/`](docs/).

Run the Python test suite with:

```bash
python -m unittest discover -s tests -v
```

### Inputs

Reference files are stored under `data/genome_eukaryote/`. Tool outputs are
stored under `data/origin/eukaryote_output_tools/`. File names and required
schemas are documented in [`docs/data-contracts.md`](docs/data-contracts.md).

### Outputs

The notebooks write organized outputs under `2_run/output/eukaryotic/`:

- `audit/gff_gffread/<organism>/`: GFF/GFFread processing audit;
- `audit/reference/<organism>/`: NCBI transcript and locus reference audit;
- `tool_preparation/tables/` and `tool_preparation/plots/`: cleaned tool
  tables, exact RNA membership tables, summaries, and UpSet plots;
- `functional_analysis/tables/` and `functional_analysis/plots/`: all-RNA,
  first-transcript, and canonical analyses;
- `functional_analysis/uniprot_cache/`: verified local UniProt snapshots used
  for canonical mapping.

Generated outputs are analysis artifacts, not source inputs. The scientific
definitions and denominator rules are described in
[`docs/eukaryotic-workflow.md`](docs/eukaryotic-workflow.md).

### Project layout

```text
1_setup/                 locked environment and setup notes
2_run/notebooks/         eukaryotic and prokaryotic notebooks
2_run/output/            generated analysis artifacts
data/                    reference genomes and tool outputs
docs/                    workflow, schema, and roadmap documentation
src/benchannot/          tested parsing and analysis functions
tests/                   unit tests for the extracted logic
```

### Scope and status

The eukaryotic workflow is the validated active path. The prokaryotic
notebooks have not been changed in this refactor; correctness and
reproducibility priorities for a future pass are listed in
[`docs/prokaryotic-roadmap.md`](docs/prokaryotic-roadmap.md).

## Nextflow pipeline

Orchestrates annotation tools and produces reproducible results for
benchmarking and downstream analysis.

### Inputs

Place FASTA files (`.fna`) in the `data/` directory.

### Run

```bash
nextflow run main.nf --bakta_db_dir path/to/bakta/db
```

Bakta requires a specific database (version 6), downloadable
[here](https://zenodo.org/records/14916843). After downloading and extracting,
provide the path with `--bakta_db_dir`.

### Outputs

- Annotated results for each genome are stored in the `results/` directory.
- Intermediate files are stored in the `work/` directory.

### Modules

- `modules/prokka.nf` — Prokka 1.14.6 ([GitHub](https://github.com/tseemann/prokka), [Docker Image](https://hub.docker.com/r/staphb/prokka))
- `modules/bakta.nf` — Bakta 1.11.3 ([GitHub](https://github.com/oschwengers/bakta), [Docker Image](https://hub.docker.com/r/oschwengers/bakta))
- `modules/gffread.nf` — processes GFF/GTF files (e.g., extracting transcript sequences)
- `modules/kofamscan.nf` — assigns KEGG Orthologs using HMM profiles

### Requirements

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)
- Input genome files in FASTA format (`.fna`) placed in the `data/` directory

### Notes

- Pipeline ensures reproducibility and scalability using Docker containers.
- Additional tools for functional annotation (e.g., InterProScan, eggNOG-mapper, Funannotate) will be integrated in future updates.
- All intermediate files remain in `work/`, while final structured results are under `results/`.