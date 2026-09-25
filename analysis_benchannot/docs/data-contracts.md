# Data contracts

BenchAnnot keeps archived tool outputs, newly reproduced tool outputs, shared
references, and downstream analysis artifacts in separate namespaces. This
document defines the paths and file-level assumptions currently enforced by the
pipeline, notebooks, and Python modules. Where producers and consumers do not
yet agree, the mismatch is documented rather than hidden by a fallback.

## Path authority and dataset selection

`src/benchannot/paths.py` is the authority for Python analysis paths.
`project_paths(dataset)` accepts exactly `"origin"` or `"reproduced"`; any
other value raises `ValueError`.

Each analysis notebook has one explicit selector:

```python
DATASET = "origin"  # "origin" or "reproduced"
PATHS = project_paths(DATASET)
```

Changing `DATASET` changes both the tool-input root and analysis-output root.
It does not change shared reference paths. There is no fallback from a missing
`reproduced` input to `origin`, or in the opposite direction.

The selected layout is:

```text
<repository-root>/
    data/
        genome_eukaryote/                  shared references
        genome_prokaryote/                 shared references
        origin/                            archived tool outputs
        reproduced/                        newly generated tool outputs
        database/                          pipeline databases
        runtime/                           pipeline runtimes
    analysis_benchannot/
        2_run/output/<dataset>/            downstream analysis artifacts
```

`project_paths()` resolves this layout without using the process working
directory and without creating files or directories. Analysis writers create
their output directories only after required-input validation succeeds.

Nextflow's `params.annotation_type` is independent of `DATASET`: it selects
`prokaryote`, `eukaryote`, or `both`, not `origin` or `reproduced`.

## Directory ownership

| Location | Owner and role |
| --- | --- |
| `data/genome_eukaryote/` | shared eukaryotic reference inputs |
| `data/genome_prokaryote/` | shared prokaryotic reference inputs |
| `data/origin/` | archived outputs used by the original analysis |
| `data/reproduced/` | outputs from current pipeline or external reruns |
| `data/database/` | database images and PGAP installation data |
| `data/runtime/` | managed runtime files, including PGAP Python |
| `analysis_benchannot/2_run/output/<dataset>/` | notebook-generated analysis artifacts |

`origin` is intended to be immutable. This is a project rule, not a filesystem
permission: users must not target it with `--outdir` or write analysis results
there. Nextflow defaults to `${projectDir}/data/reproduced`.

## Input preflight

Python notebooks call `validate_required_inputs()` with the files consumed by
the next stage. The validator:

- resolves relative paths from the repository root, not `cwd`;
- aggregates all missing paths and non-file paths in one `FileNotFoundError`;
- reports repository-relative paths;
- does not create directories;
- does not inspect size, checksum, schema, or parseability.

Parsers remain responsible for schema and identifier validation. A directory's
existence alone does not satisfy a notebook input contract.

Nextflow separately validates samplesheets, required row fields, referenced
files, annotation type, and prokaryotic numeric fields. Samplesheet paths are
repository-relative and are resolved against `projectDir`.

## Shared references

### Eukaryotes

For each `<organism>`:

```text
data/genome_eukaryote/<organism>.fna
data/genome_eukaryote/<organism>.gff
data/genome_eukaryote/<organism>.faa
```

The FAA defines the functional reference universe. The GFF supplies transcript
and locus mappings. The FNA is a pipeline input and provenance record. GFFread
outputs audit sequence extraction; they do not define the reference
denominator.

The eukaryotic samplesheet columns are:

```text
sample_id,genome_fasta,reference_gff,reference_faa,organism_id
```

### Prokaryotes

For each `<organism>`:

```text
data/genome_prokaryote/<organism>.fna
data/genome_prokaryote/<organism>.gff
data/genome_prokaryote/<organism>.gbff
```

The prokaryotic samplesheet columns are:

```text
sample_id,genome_fasta,species,taxid,genetic_code
```

`genetic_code` is optional in validation and defaults to `11`. The notebooks
currently derive GFF and GBFF names from their fixed organism IDs; those files
are not supplied through the samplesheet.

## Dataset-specific tool inputs

The roots are:

```text
data/<dataset>/eukaryote_output_tools/
data/<dataset>/prokaryote_output_tools/
```

All files required by one notebook run must belong to the selected dataset
except shared references. Mixed-dataset analysis is unsupported.

### Eukaryotic filenames

```text
eukaryote_output_tools/gffread/<organism>_gffread.faa
eukaryote_output_tools/gffread/<organism>.filtered.gff
eukaryote_output_tools/kofamscan/<organism>.kofam.txt
eukaryote_output_tools/pannzer/<organism>_pannzer.txt
eukaryote_output_tools/eggnog/<organism>_eggnog.emapper.annotations
eukaryote_output_tools/interproscan/<organism>.interpro.tsv
eukaryote_output_tools/interproscan/<organism>.interpro.gff3
```

The analysis audit reconstructs its own filtered GFF and consumes the GFFread
FAA. The pipeline-published filtered GFF is not a downstream notebook input.

PANNZER is an external web-service result. No PANNZER Nextflow process exists.
Both datasets therefore require a separately supplied PANNZER file at the path
above. Missing reproduced PANNZER output is an error, not a reason to use the
origin file.

### Prokaryotic filenames

```text
prokaryote_output_tools/prokka/<organism>.gff
prokaryote_output_tools/prokka/<organism>.gbk
prokaryote_output_tools/bakta/<organism>.gff3
prokaryote_output_tools/bakta/<organism>.gbff
prokaryote_output_tools/eggnog/<organism>.emapper.decorated.gff
prokaryote_output_tools/pgap/<organism>_annot.gff
prokaryote_output_tools/pgap/<organism>_annot.gbk
```

The PGAP paths are the current notebook consumer contract and match the archived
origin layout. The Nextflow PGAP process currently publishes a
`<sample_id>_pgap/` directory instead of flattening and renaming individual
files. Reproduced PGAP compatibility is therefore unresolved. Until the real
published contents are verified and a publication rule is implemented, the
notebooks intentionally fail preflight rather than search alternative paths.

## Analysis outputs

Analysis artifacts are namespaced by the selected dataset:

```text
analysis_benchannot/2_run/output/<dataset>/eukaryotic/
analysis_benchannot/2_run/output/<dataset>/prokaryotic/
```

Legacy artifacts or stored notebook output that mention
`2_run/output/eukaryotic/` or `2_run/output/prokaryotic/` without `<dataset>`
do not define the current output contract.

### Eukaryotic output namespaces

```text
eukaryotic/audit/gff_gffread/<organism>/
eukaryotic/audit/reference/<organism>/
eukaryotic/tool_preparation/tables/
eukaryotic/tool_preparation/plots/
eukaryotic/functional_analysis/tables/
eukaryotic/functional_analysis/plots/
eukaryotic/functional_analysis/uniprot_cache/
```

Prepared references include:

```text
audit/reference/<organism>/prepared_reference/transcript_reference.tsv
audit/reference/<organism>/prepared_reference/locus_reference.tsv
```

### Prokaryotic output namespace

The gene-prediction notebook writes:

```text
prokaryotic/<organism>_prokka_cleaned.gff
prokaryotic/<organism>_gene_prediction.png
prokaryotic/<organism>_concordant_CDS.csv
```

The functional-analysis and identifier-matching notebooks currently display
results but do not persist their summary tables or figures.

## Eukaryotic table contracts

The benchmark uses exact identifiers and explicit schemas to prevent silent row
expansion, identifier normalization, or denominator changes.

The prepared transcript table must have unique, non-empty `RNA_ID` values. The
NCBI FAA identifiers are parsed into versioned `protein_id`, exact `RNA_ID`, and
locus metadata from the associated reference annotation.

| Tool | Source file | Join key | Functional rule |
| --- | --- | --- | --- |
| KofamScan | `kofamscan/<organism>.kofam.txt` | `RNA_ID` | retain significant rows where `marker == "*"` |
| PANNZER | `pannzer/<organism>_pannzer.txt` | `RNA_ID` | retain one cleaned description row per exact ID |
| eggNOG | `eggnog/<organism>_eggnog.emapper.annotations` | `RNA_ID` | retain one cleaned description row per exact ID |
| InterProScan | `interproscan/<organism>.interpro.tsv` | `RNA_ID` | select the numerically smallest eligible e-value with stable tie ordering |

Cleaned tool tables contain exactly `RNA_ID` and the tool description column,
with unique non-empty `RNA_ID` values. Tool IDs outside the reference universe
are rejected rather than silently dropped.

The cleaning audit distinguishes:

- `gffread_input_rna_ids`: unique non-mitochondrial IDs submitted to tools;
- `tool_identified_rna_ids`: submitted IDs identified by raw tool output;
- `identified_outside_submission`: raw-output IDs absent from the submission;
- `retained_rna_ids`: cleaned IDs retained in the NCBI reference universe;
- `represented_loci`: unique reference loci represented by retained IDs.

`<organism>_functional_annotations.tsv` contains:

```text
protein_id  RNA_ID  locus_tag  Reference
Kofam  Kofam_returned
Pannzer  Pannzer_returned
EggNOG  EggNOG_returned
InterProScan  InterProScan_returned
```

`<tool>_returned == True` means at least one cleaned row exists for the exact
`RNA_ID`; it does not mean that the description is informative or biologically
correct. Coverage and UpSet membership use returned flags. Functional
classification uses description content.

Canonical selection preserves one row per selected locus and retains
`resolution_status` and `submission_status`. It does not use a tool description
to choose a canonical protein. Canonical coverage columns are:

```text
Source  RNA_ID count  Represented loci  Recovered canonical RNA_ID
```

## Provenance

Canonical UniProt cache files live at:

```text
eukaryotic/functional_analysis/uniprot_cache/<organism>_<proteome_id>.json
```

The cache records and verifies `proteome_id`, request and response URLs,
retrieval time, response SHA-256, and record count.

Tool/database versions are pinned in Nextflow configuration, but the pipeline
does not yet write a complete per-run manifest containing input checksums,
effective parameters, image checksums, container digests, and Nextflow version.
Origin is intended to represent archived study outputs, but the setup notebook
still contains placeholder archive URLs and no complete checksum manifest.

## Current non-contracts

- Stored notebook outputs and historical absolute paths are not path contracts.
- No dataset fallback or automatic mixed-dataset analysis exists.
- PANNZER reproduction is external to Nextflow.
- Reproduced PGAP layout compatibility is unresolved.
- Python preflight does not validate file size, checksum, or schema.
- Prokaryotic ORF keys currently omit `seqid`; multi-contig correctness is not
  guaranteed by the current coordinate key.
- Prokaryotic concordant CSVs contain only the current concordant subset; they
  are not a full-reference denominator.
- Live prokaryotic UniProt mappings are not persistently cached and may collapse
  one-to-many mappings.
- Description informativeness is not biological correctness or semantic
  equivalence.
- `origin` immutability is a documented rule, not a technical write barrier.

Before any merge or export, code should validate required columns, key
non-emptiness, expected uniqueness/cardinality, reference containment, row
counts, and denominator definitions while preserving the distinction between
missing, empty, placeholder, vague, and informative descriptions.
