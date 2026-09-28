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

Each notebook explicitly selects `origin` or `reproduced`. Paths below are
relative to `analysis_benchannot/2_run/output/<dataset>/eukaryotic/`; raw tool
inputs come from `data/<dataset>/eukaryote_output_tools/`.

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

Only significant Kofam hits (`marker == "*"`) enter the downstream functional
table and figures. Presence plots count a returned row even when its
description is empty, `-`, or otherwise non-informative. The combined table
therefore stores an explicit `<tool>_returned` Boolean flag for every tool in
addition to the description column.

The UpSet reference universe is first restricted to exact RNA_ID values present
in the GFFread FASTA, then reduced to the first submitted transcript per
`locus_tag`. Tool memberships are restricted to those representatives, so
RNA_ID values not returned by a tool remain explicit no-hits. The combined
two-organism figure uses significant Kofam hits and the fixed source-column
order Reference, Kofam, Pannzer, EggNOG, InterProScan.

The processing summary is a 4-by-2 panel grid: tools are rows,
*S. cerevisiae* is the first column, and *D. melanogaster* is the second. Every
panel shows four cumulative stages: GFFread input RNA_ID, tool-identified
RNA_ID, retained RNA_ID after cleaning and exact reference matching, and
retained representative RNA_ID after selecting the first submitted transcript
per locus.

### 3. Analyze functional annotations

`2_run/notebooks/eukaryotic/3_functional_analysis.ipynb` has two analyses:

1. **Submitted common-tool comparison:** restrict the reference to RNA_ID
   values present in the GFFread FASTA, select the first submitted transcript
   per locus, and retain only proteins returned by Kofam, Pannzer, EggNOG, and
   InterProScan. This adds a common-hit filter to the same representative
   universe used by the UpSet figure.
2. **Reviewed Swiss-Prot canonical subset:** canonical selection is external
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
first cell. Paths are resolved by the installed `benchannot.paths` module and
do not depend on the Jupyter working directory. UniProt data are stored in a
local cache with request metadata, record counts, and a SHA-256 response hash.
Unit tests do not access the network.
