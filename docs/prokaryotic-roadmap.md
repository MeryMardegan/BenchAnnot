# Recommended prokaryotic notebook changes

## Correctness priorities

1. Include `seqid` in every CDS key. Replace `Start_End_Strand` identifiers
   with an explicit `(seqid, start, end, strand)` key so equal coordinates on
   different contigs or replicons cannot be merged.
2. Preserve one-to-many UniProt and UniRef mappings. Do not collapse mapping
   responses with `dict(zip(...))`; classify each query as unique, ambiguous,
   or unmapped and resolve candidates only with explicit sequence and organism
   evidence.
3. Parse PGAP inference values by structured prefix. Only actual
   `similar to AA sequence:RefSeq:<accession>` evidence should enter the RefSeq
   protein mapping route; motif identifiers must remain separate evidence.
4. Evaluate against the full NCBI RefSeq CDS universe. The current
   `concordant_only=True` export selects CDSs recovered by every tool and can
   inflate downstream functional and identifier performance. Retain every
   reference CDS, left-join each tool, add explicit `<tool>_returned` flags,
   and report concordant-only results only as a sensitivity analysis.
5. Validate merge cardinality. Require unique compound keys per source and use
   `validate="one_to_one"` where appropriate so duplicate records cannot
   silently expand denominators.
6. Use strain-specific taxonomy and proteomes for identifier mapping. Record
   all candidates and the evidence used to select a target instead of taking
   the first returned candidate.
7. Rename functional categories to describe informativeness rather than
   biological correctness. Prefer `both informative`, `both non-informative`,
   `tool-only informative`, and `reference-only informative`, or retain the
   current labels only with an explicit caveat.
8. Document eggNOG structural provenance. A decorated eggNOG GFF generally
   inherits an upstream submitted CDS set and should not be presented as an
   independent gene predictor without naming that source.

## Reproducibility and maintenance

1. Define and validate named schemas for every GFF, TSV, and mapping table;
   replace positional parsing and broad substring-based column deletion with
   explicit contracts.
2. Cache external mapping responses with sorted requests, bounded retries,
   timeouts, release/access dates, taxon, namespace, checksums, and per-query
   status. Fail clearly when a mapping job is incomplete.
3. Discover the project root from `pyproject.toml` rather than assuming the
   notebook working directory.
4. Organize outputs under per-analysis `tables/`, `plots/`, `intermediate/`,
   and `cache/` directories. Save the exact plotted tables alongside static
   figures.
5. Record input checksums, reference assembly accessions, tool and database
   versions, package versions, parameters, and execution date in a manifest.
6. Use stable tool ordering and shared color constants. State denominators in
   figure captions, label diverging axes explicitly, replace `Ground Truth`
   with `NCBI RefSeq benchmark reference`, and avoid calling accession
   concordance `accuracy` unless exact sequence identity validates it.

## Notebook-specific scope

- `1_gene_prediction.ipynb`: compound CDS keys, eggNOG provenance, full-reference
  export, merge validation, and structural metric exports.
- `2_functional_analysis.ipynb`: full-reference denominators, explicit
  informativeness terminology, saved plot tables, and reproducible figures.
- `3_id_matching.ipynb`: structured PGAP inference parsing, ambiguity-preserving
  UniProt/UniRef mappings, strain-specific filtering, cache provenance, and
  evidence-based mapping statuses.
