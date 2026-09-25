# BenchAnnot

BenchAnnot is a reproducible framework for benchmarking genome annotation
tools against curated reference annotations. This repository has two parts:

- **Nextflow pipeline** — `main.nf` multiple functional annotation tools for
 genomes using Nextflow DSL2. It is designed for reproducible, scalable benchmarking 
 and downstream analysis.
- **Analysis workflow** — notebooks and Python modules that audit, clean, and
  compare annotation tool outputs (eukaryotic and prokaryotic), under
  [`analysis_benchannot/`](analysis_benchannot/).


## Nextflow pipeline

Orchestrates annotation tools and produces reproducible results for
benchmarking and downstream analysis.

### Inputs

Place FASTA files (`.fna`) in `data/genome_eukaryote/` or
`data/genome_prokaryote/` and describe them in the corresponding samplesheet
under `assets/`.
Required files for eukaryotes:
- Genome FASTA files: `data/genome_eukaryote/<sample>.fna`
- Matching GFF files: `data/genome_eukaryote/<sample>.gff`
- Reference protein files: `data/genome_eukaryote/<sample>.faa`

### Run

#### How you can starts
mkdir -p "$PWD/apptainer_cache/tmp"

export SINGULARITY_TMPDIR="$PWD/apptainer_cache/tmp"
export TMPDIR="$PWD/apptainer_cache/tmp"

```bash
nextflow run main.nf --annotation_type both
```

Bakta requires a specific database (version 6), downloadable
[here](https://zenodo.org/records/14916843). After downloading and extracting,
provide a SquashFS image with `--bakta_db`, or leave the parameter unset so
BenchAnnot prepares the pinned database release.

### Outputs

- Annotated results are stored under `data/reproduced/eukaryote_output_tools/`
  and `data/reproduced/prokaryote_output_tools/` by default.
- Intermediate files are stored in the `work/` directory.

### Modules

**Prokaryote**
- `modules/prokka.nf` — Prokka 1.14.6 ([GitHub](https://github.com/tseemann/prokka), [Docker Image](https://hub.docker.com/r/staphb/prokka))
- `modules/bakta.nf` — Bakta 1.11.3 ([GitHub](https://github.com/oschwengers/bakta), [Docker Image](https://hub.docker.com/r/oschwengers/bakta))
- `modules/eggnog.nf` — eggnog-mapper-v2 2.1.13 ([GitHub](https://github.com/eggnogdb/eggnog-mapper))
- `modules/pgap.nf`— Prokaryotic Genome Annotation Pipeline 2026-06-18.build8602 ([GitHub](https://github.com/ncbi/pgap))
**Eukaryote**
- `modules/gffread.nf` — processes GFF/GTF files (e.g., extracting transcript sequences)
- `modules/kofamscan.nf` — assigns KEGG Orthologs using HMM profiles
- `modules/interproscan.nf` — functional domain and GO/pathway annotation
- `modules/eggnog.nf` — Orthology-based functional annotation.

### Parameters

Some modules require specific database directories to be provided as parameters.
Set database locations in [nextflow.config](nextflow.config):

#### Prokaryote

**Bakta**
Requires Bakta database version 6. 
Download from [Zenodo](https://zenodo.org/records/14916843) and package it as a
SquashFS image, then use `--bakta_db /absolute/path/to/bakta.sqsh`. When the
parameter is unset, the pipeline prepares its pinned release under
`data/database/bakta/`.

**EggNog Mapper**
Requires Eggnog database files.
The configured image is `data/database/eggnog/eggnog_2026-09.sqsh`. If it is
absent, the pipeline prepares and packages the database.

**PGAP**
Requires a local PGAP installation (including its databases and container images). Follow the official Quick Start to install PGAP and download required data: https://github.com/ncbi/pgap/wiki/Quick-Start. Once installed, provide the absolute installation path with:
`--pgap_dir /absolute/path/to/pgap`

#### Eukaryote

- `params.ips_db` optionally supplies an InterProScan SquashFS image.
- `params.eggnog_db` identifies the eggNOG SquashFS image.
- `params.kofamscan_db` optionally supplies a KofamScan SquashFS image.

#### Run command
Choose eukaryotic, prokaryotic, or combined annotation:
**Only Eukaryote**
```bash
nextflow run main.nf --annotation_type eukaryote
```
**Only Procaryote**
```
nextflow run main.nf --annotation_type prokaryote
```
**Run all**
nextflow run main.nf

### Requirements

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)
- Databases downloaded locally for InterProScan, KofamScan, and eggNOG-mapper
- Input genome files in FASTA format (`.fna`) placed in the `data/` directory

### Outputs
For prokaryotes:
- Per-tool outputs under `data/reproduced/prokaryote_output_tools/`.

### Documentation
- eggNOG-mapper: [docs/eggnog.md](docs/eggnog.md)
- InterProScan: [docs/interproscan.md](docs/interproscan.md)
- KofamScan: [docs/KofamScan.md](docs/KofamScan.md)

### Notes

- Pipeline ensures reproducibility and scalability using Docker containers.
- Additional tools for functional annotation (e.g., InterProScan, eggNOG-mapper, Funannotate) will be integrated in future updates.
- All intermediate files remain in `work/`, while final structured results are
  published under `data/reproduced/` by default.

## Analysis workflow

The current refactor focuses on the eukaryotic post-processing workflow for
*Saccharomyces cerevisiae* and *Drosophila melanogaster*. Prokaryotic notebooks
remain available, but their planned improvements are documented separately.

Run the eukaryotic notebooks from `analysis_benchannot/2_run/notebooks/eukaryotic/` in this order:

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
conda-lock install --name benchannot analysis_benchannot/1_setup/conda-lock.yml
conda activate benchannot
python -m pip install -e analysis_benchannot/
python -m ipykernel install --user --name benchannot --display-name "Python (BenchAnnot)"
```

Then open the project in Jupyter and execute the three eukaryotic notebooks
top to bottom. Detailed setup and workflow notes are in
[`analysis_benchannot/docs/`](analysis_benchannot/docs/).

Run the Python test suite with:

```bash
python -m unittest discover -s analysis_benchannot/tests -v
```

### Inputs

Reference files are stored under `data/genome_eukaryote/`. Each notebook
explicitly selects `origin` or `reproduced`; tool inputs are then read from
`data/<dataset>/eukaryote_output_tools/`. File names and required schemas are documented in
[`analysis_benchannot/docs/data-contracts.md`](analysis_benchannot/docs/data-contracts.md).

### Outputs

The notebooks write organized outputs under
`analysis_benchannot/2_run/output/<dataset>/eukaryotic/`:

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
[`analysis_benchannot/docs/eukaryotic-workflow.md`](analysis_benchannot/docs/eukaryotic-workflow.md).

### Project layout

```text
data/                        shared references, datasets, databases, runtimes
analysis_benchannot/
    1_setup/                 locked environment and setup notes
    2_run/notebooks/         eukaryotic and prokaryotic notebooks
    2_run/output/<dataset>/  generated analysis artifacts
    docs/                    workflow, schema, and roadmap documentation
    src/benchannot/          tested parsing and analysis functions
    tests/                   unit tests for the extracted logic
```

### Scope and status

The eukaryotic workflow is the validated active path. The prokaryotic notebooks
now use the same explicit dataset and output namespaces, while scientific
correctness and reproducibility priorities for a future pass remain listed in
[`analysis_benchannot/docs/prokaryotic-roadmap.md`](analysis_benchannot/docs/prokaryotic-roadmap.md).
