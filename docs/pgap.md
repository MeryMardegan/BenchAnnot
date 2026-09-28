# PGAP

## Purpose

The NCBI Prokaryotic Genome Annotation Pipeline provides a second prokaryotic
annotation path. BenchAnnot pins PGAP release `2026-06-18.build8602`.

## Processes and runtime

| Process | Label | Execution environment |
| --- | --- | --- |
| `PREPARE_PGAP_PYTHON` | `pgap_python_prepare` | managed Python 3.11 runtime |
| `PREPARE_PGAP` | `pgap_prepare` | `pgap.py` launches Singularity |
| `PGAP` | `pgap` | `pgap.py` launches the pinned PGAP SIF |

PGAP differs from the other modules: Nextflow does not assign a container to
the annotation process. The `pgap.py` wrapper invokes the PGAP container itself.

## Managed installation

The resolver checks for existing files at:

```text
data/runtime/pgap-python-3.11/
data/database/pgap/
data/database/pgap/pgap_2026-06-18.build8602.sif
```

The reuse check confirms that the Python path is a file and that the PGAP
directory and non-empty SIF exist; it does not verify the installed release or
SIF checksum. If they are absent, BenchAnnot bootstraps micromamba from its
mutable `linux-64/latest` endpoint, installs the Python 3.11 series without a
lockfile, downloads the pinned `pgap.py`, and asks it to prepare the matching
Singularity distribution.

Although `pgap_python_version` is configurable, runtime assertions currently
require Python 3.11. Changing that parameter alone is unsupported.

## Input

```text
tuple(sample_id, genome_fasta, species)
```

`taxid` and `genetic_code` are validated in the prokaryotic samplesheet but are
not currently passed to PGAP. The run enables taxonomic checking and automatic
taxon correction and disables network access during annotation.

## Output

The complete PGAP output directory is published as:

```text
<outdir>/prokaryote_output_tools/pgap/<sample_id>_pgap/
```

## Downstream limitation

The prokaryotic notebooks currently expect flat files named
`<sample_id>_annot.gff` and `<sample_id>_annot.gbk`. BenchAnnot does not yet
select and rename those files from the published PGAP directory. Reproduced
PGAP output therefore does not directly satisfy the notebook input contract.

## Resources

Python preparation requests 1 CPU, 8 GB, and 1 hour. PGAP distribution
preparation inherits 10 CPUs, 128 GB, and 24 hours. Annotation uses the same
default CPU count and requests 128 GB for 24 hours.

## Example

```bash
nextflow run main.nf --annotation_type prokaryote
```
