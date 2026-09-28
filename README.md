# BenchAnnot

BenchAnnot is a reproducible framework for benchmarking genome annotation tools
against curated reference annotations. The repository contains two connected
components:

- **Nextflow pipeline:** runs structural and functional annotation tools through
  a DSL2 workflow.
- **Analysis workflow:** audits, cleans, and compares annotation outputs through
  notebooks and Python modules under [`analysis_benchannot/`](analysis_benchannot/).

## Nextflow pipeline

The pipeline supports prokaryotic annotation, eukaryotic functional annotation,
or both branches in one execution. The current configuration targets Slurm and
Apptainer and can prepare the required local databases when valid images are not
already available.

### Requirements

- [Nextflow](https://www.nextflow.io/) with DSL2 support;
- Slurm access using the configured `SP2` queue;
- Apptainer or Singularity with SquashFS image-mount support;
- network access while containers, databases, or the PGAP runtime are prepared;
- input files referenced by the selected samplesheet.

The executor, queue, resources, containers, and runtime configuration are
described in [`docs/execution.md`](docs/execution.md).

### Inputs

Input files are declared in CSV samplesheets rather than discovered by scanning
directories. Relative paths are resolved from the repository root.

The prokaryotic samplesheet defaults to `assets/prokaryote_samplesheet.csv`:

```text
sample_id,genome_fasta,species,taxid,genetic_code
```

`sample_id`, `genome_fasta`, `species`, and `taxid` are required. `genetic_code`
is optional and defaults to `11`. `taxid` and `genetic_code` must be integers,
but are not currently passed to annotation processes.

The eukaryotic samplesheet defaults to `assets/eukaryote_samplesheet.csv`:

```text
sample_id,genome_fasta,reference_gff,reference_faa,organism_id
```

All five columns are required and all three paths must identify existing files.
The current annotation branch passes the genome FASTA and reference GFF to
GFFread; `reference_faa` and `organism_id` are retained as samplesheet metadata
but are not process inputs.

Use unique `sample_id` values containing only filename-safe characters. The
pipeline interpolates these values directly into output names and does not
currently validate uniqueness or sanitize them.

### Run

Prepare a local temporary directory for Apptainer:

```bash
mkdir -p "$PWD/apptainer_cache/tmp"
export SINGULARITY_TMPDIR="$PWD/apptainer_cache/tmp"
export TMPDIR="$PWD/apptainer_cache/tmp"
```

Run one branch or both:

```bash
nextflow run main.nf --annotation_type eukaryote
nextflow run main.nf --annotation_type prokaryote
nextflow run main.nf --annotation_type both
```

The default is `--annotation_type both`. Override a samplesheet or output root
when needed:

```bash
nextflow run main.nf \
  --annotation_type eukaryote \
  --eukaryote_samplesheet /path/to/eukaryotes.csv \
  --outdir /path/to/reproduced
```

### Outputs

`--outdir` defaults to `data/reproduced`. All annotation modules publish below
that root:

```text
<outdir>/eukaryote_output_tools/
<outdir>/prokaryote_output_tools/
```

Database images and managed runtimes remain under `data/database/` and
`data/runtime/`; they are not analysis results. Nextflow task outputs and logs
remain under `work/`, which is also required by `-resume`.

The analysis distinguishes archived source outputs in `data/origin/` from
pipeline reruns in `data/reproduced/`. See
[`analysis_benchannot/docs/data-contracts.md`](analysis_benchannot/docs/data-contracts.md)
for the complete boundary between pipeline and notebooks.

### Tools

| Branch | Tool | Documentation |
| --- | --- | --- |
| Prokaryotic | Prokka | [`docs/prokka.md`](docs/prokka.md) |
| Prokaryotic | Bakta | [`docs/bakta.md`](docs/bakta.md) |
| Prokaryotic | PGAP | [`docs/pgap.md`](docs/pgap.md) |
| Both | eggNOG-mapper | [`docs/eggnog.md`](docs/eggnog.md) |
| Eukaryotic | GFFread | [`docs/gffread.md`](docs/gffread.md) |
| Eukaryotic | KofamScan | [`docs/KofamScan.md`](docs/KofamScan.md) |
| Eukaryotic | InterProScan | [`docs/interproscan.md`](docs/interproscan.md) |

Database preparation, SquashFS packaging, container mounts, and configured
resources are documented in [`docs/execution.md`](docs/execution.md).

## Analysis workflow

The current refactor focuses on the eukaryotic post-processing workflow for
*Saccharomyces cerevisiae* and *Drosophila melanogaster*. Prokaryotic notebooks
remain available, but their planned improvements are documented separately.

Run the eukaryotic notebooks from `analysis_benchannot/2_run/notebooks/eukaryotic/` in this order:

1. `1_audit_prepare_reference.ipynb` audits GFF/GFFread processing and builds
   an independent NCBI reference universe.
2. `2_prepare_output_tools.ipynb` cleans Kofam, Pannzer, EggNOG, and
   InterProScan outputs using exact `RNA_ID` values.
3. `3_functional_analysis.ipynb` compares the first submitted transcript per
   locus among records returned by all four tools, then evaluates the reviewed
   Swiss-Prot canonical subset separately.

The NCBI FAA and GFF define the complete functional reference. The primary
tool comparison restricts that reference to RNA_ID values present in the
GFFread FASTA and selects one submitted transcript per locus.

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
- `functional_analysis/tables/` and `functional_analysis/plots/`: submitted
  first-transcript common-tool comparison and canonical analysis;
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
