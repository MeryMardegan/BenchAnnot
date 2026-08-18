# Data contracts

The benchmark uses exact identifiers and explicit schemas at each boundary.
These contracts prevent silent row expansion, identifier normalization, or
denominator changes.

## Reference inputs

| Path | Required role |
| --- | --- |
| `data/genome_eukaryote/<organism>.gff` | source annotation audited against the NCBI reference |
| `data/genome_eukaryote/<organism>.faa` | authoritative protein and description reference |
| `data/genome_eukaryote/<organism>.fna` | genome sequence input/provenance |
| `data/origin/eukaryote_output_tools/gffread/<organism>_gffread.faa` | GFFread extraction inventory used for loss auditing |

The NCBI FAA identifiers are parsed into versioned `protein_id`, exact
`RNA_ID`, and locus metadata from the associated reference annotation. The
prepared transcript table must have unique, non-empty `RNA_ID` values.

## Tool inputs

| Tool | Source file | Join key | Functional rule |
| --- | --- | --- | --- |
| KofamScan | `kofamscan/<organism>.kofam.txt` | `RNA_ID` | retain significant rows where `marker == "*"` for functional analysis |
| Pannzer | `pannzer/<organism>_pannzer.txt` | `RNA_ID` | retain one cleaned row per exact ID |
| EggNOG | `eggnog/<organism>_eggnog.emapper.annotations` | `RNA_ID` | retain one cleaned row per exact ID |
| InterProScan | `interproscan/<organism>.interpro.tsv` | `RNA_ID` | select the numerically smallest e-value with stable tie ordering |

All cleaned tool tables must contain exactly `RNA_ID` and the tool's
description column, with unique non-empty `RNA_ID` values. Tool IDs not present
in the reference are rejected rather than silently dropped.

The cleaning summary distinguishes the three plotted processing stages:

- `gffread_input_rna_ids`: unique non-mitochondrial IDs in the GFFread FASTA
  supplied to tools;
- `tool_identified_rna_ids`: unique non-mitochondrial IDs identified by the raw
  tool output and belonging to that GFFread input;
- `identified_outside_submission`: raw-output IDs absent from the submission;
- `retained_rna_ids`: cleaned IDs retained in the NCBI reference universe;
- `represented_loci`: unique reference `locus_tag` values represented by the
  retained IDs.

The figure plots only `gffread_input_rna_ids`, `tool_identified_rna_ids`, and
`retained_rna_ids`. Raw hit-row counts remain available in the TSV audit but are
not plotted because tools emit different numbers of evidence rows per RNA_ID.

## Combined functional table

`<organism>_functional_annotations.tsv` contains the reference provenance and
one description plus one presence flag for each tool:

```text
protein_id  RNA_ID  locus_tag  Reference
Kofam  Kofam_returned
Pannzer  Pannzer_returned
EggNOG  EggNOG_returned
InterProScan  InterProScan_returned
```

The description and presence columns have different meanings:

- `<tool>_returned == True` means at least one cleaned line was returned for
  that exact `RNA_ID`;
- the description may still be empty, `-`, or non-informative;
- functional classification uses description content;
- coverage and UpSet membership use the returned flag.

## Canonical tables

The canonical analysis writes mapping, selection, coverage, and functional input
tables. The canonical selection must preserve one row per selected locus and
retain `resolution_status` and `submission_status`. It must not use a tool
description to select a canonical protein.

The coverage table has these columns:

```text
Source  RNA_ID count  Represented loci  Recovered canonical RNA_ID
```

All counts are exact integer counts. `RNA_ID count` is not a count of
informative descriptions.

## Validation requirements

Before a merge or export, validate:

- required columns are present;
- key columns are non-empty and unique where the contract is one-to-one;
- tool IDs are contained in the reference universe;
- row counts and denominator definitions are printed and exported;
- missing, empty, placeholder, and vague descriptions remain distinguishable
  from returned-row presence.
