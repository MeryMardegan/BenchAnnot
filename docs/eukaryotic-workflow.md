# Eukaryotic workflow

This document defines the current eukaryotic benchmark workflow. The three
notebooks are intentionally separated so that reference construction,
tool-output preparation, and functional comparison can be audited independently.

## Reference universe

The reference is built from the NCBI protein FASTA in
`data/genome_eukaryote/` after mitochondrial records are excluded. The GFF and
GFFread FASTA are used to audit processing losses, not to replace the NCBI
reference.

The transcript reference has one row per exact `RNA_ID` and contains:

- `protein_id`: versioned NCBI protein accession;
- `RNA_ID`: exact transcript key used to join tool outputs;
- `gene_id`: gene identifier, when available;
- `locus_tag`: locus-level aggregation key;
- `product`: NCBI reference description.

The locus reference has one row per `locus_tag`, with deterministic counts and
lists of its transcript and protein records.

The current reference sizes are 6,002 transcript records and 6,002 loci for
*S. cerevisiae*, and 30,789 transcript records and 13,973 loci for
*D. melanogaster*. The exact values are also written to each reference audit
workbook.

## Notebook sequence

### 1. Audit and prepare the reference

`2_run/notebooks/eukaryotic/1_audit_prepare_reference.ipynb` writes two audit
workbooks per organism:

- `audit/gff_gffread/<organism>/gff_gffread_audit.xlsx`;
- `audit/reference/<organism>/ncbi_reference_audit.xlsx`.

It records pre-GFFread filtering, GFFread extraction losses, mitochondrial
removals, transcript-level reference records, and locus-level reference records.

### 2. Prepare tool outputs

`2_run/notebooks/eukaryotic/2_prepare_output_tools.ipynb` reads the exact
reference and the raw outputs from KofamScan, Pannzer, EggNOG, and InterProScan.
It writes cleaned per-tool tables and the combined
`<organism>_functional_annotations.tsv` table.

Kofam produces two retained variants for overlap analysis:

- all Kofam hits;
- significant hits where `marker == "*"`.

Only significant Kofam hits enter the downstream functional table. Presence
plots count any returned row, even when its description is empty, `-`, or
otherwise non-informative. The combined table therefore stores an explicit
`<tool>_returned` Boolean flag for every tool in addition to the description
column.

The notebook also writes exact RNA-level UpSet tables and plots. For
*D. melanogaster*, it additionally selects the first reference transcript in
table order for each `locus_tag` and writes a second, noncanonical
representative-per-locus UpSet analysis. The corresponding files are named
`drosophila_melanogaster_*_first_reference_transcript_{all_kofam,significant_kofam}`.
UpSet categories are displayed by cardinality from largest to smallest; the
exported membership tables keep the fixed source-column order.

The cleaning summary is a 4-by-2 panel grid: tools are rows, *S. cerevisiae*
is the first column, and *D. melanogaster* is the second. Every panel uses its
own count scale and has three stages: GFFread input RNA_ID, tool-identified
RNA_ID after mitochondrial removal, and retained RNA_ID after the tool-specific
selection and exact reference match. The exact selection rule is shown in the
panel text. `represented_loci` remains available in the summary table for the
Drosophila locus-level audit but is not mixed into the RNA_ID plot.

### 3. Analyze functional annotations

`2_run/notebooks/eukaryotic/3_functional_analysis.ipynb` has three analyses:

1. **All RNA records:** every reference `RNA_ID` is evaluated independently.
2. **First transcript per locus:** the first transcript in reference table
   order is selected as a deterministic, explicitly noncanonical sensitivity
   analysis.
3. **Reviewed Swiss-Prot canonical subset:** canonical selection is external
   to tool performance and is resolved conservatively.

Functional descriptions are non-informative when they are missing, empty,
`-`, or contain the configured vague terms such as `hypothetical`,
`uncharacterized`, `unknown`, `putative`, `predicted`, or `DUF`. A category
such as `Both Informative` describes description informativeness only; it does
not establish biological equivalence.

## Canonical mapping rules

The canonical subset uses the local verified UniProt proteome cache. A locus is
resolved only when a reviewed Swiss-Prot record provides:

1. an exact, versioned RefSeq protein cross-reference;
2. an exact amino-acid sequence match to the NCBI reference FAA; and
3. a unique candidate, or an exact equivalent-transcript tie resolved by the
   deterministic NCBI reference order.

Ambiguous mappings, sequence conflicts, absent reviewed mappings, and reviewed
canonical proteins not present in submitted GFFread input do not fall back to
another transcript. A canonical RNA absent from submitted input remains in the
canonical denominator with missing tool annotations.

The canonical coverage table reports:

- `RNA_ID count`: reference RNA records returned by the source;
- `Represented loci`: unique loci represented by those records;
- `Recovered canonical RNA_ID`: selected canonical RNA records present in the
  source's returned set.

These are separate from informative-description counts in the functional shift
tables.

## Validated run

The current validated outputs report:

| Organism | Reference RNA_ID | Reference loci | Canonical denominator |
| --- | ---: | ---: | ---: |
| *Saccharomyces cerevisiae* | 6,002 | 6,002 | 5,823 |
| *Drosophila melanogaster* | 30,789 | 13,973 | 3,847 |

For *D. melanogaster*, the canonical denominator consists of 2,386 uniquely
resolved canonical loci and 1,461 equivalent-transcript ties. The mapping table
retains the other statuses rather than silently converting them into canonical
records.

## Reproducibility

Run notebooks from the project environment and execute each notebook from its
first cell. The notebooks discover the project root from `pyproject.toml`
where supported. UniProt data are stored in a local cache with request metadata,
record counts, and a SHA-256 response hash. Unit tests do not access the network.
