# Pipeline execution

## Runtime model

BenchAnnot is a Nextflow DSL2 pipeline configured for Slurm and Apptainer. The
repository configuration currently sets:

```text
executor = slurm
queue    = SP2
cpus     = 10
memory   = 128 GB
time     = 24 hours
```

Tool labels override these defaults where preparation or annotation has a
different requirement. The site-specific executor and queue are currently in
the global configuration rather than a named profile.

## Apptainer cache

`nextflow.config` enables Singularity-compatible execution and stores pulled
images under `apptainer_cache/`. Prepare a writable temporary directory before
running:

```bash
mkdir -p "$PWD/apptainer_cache/tmp"
export SINGULARITY_TMPDIR="$PWD/apptainer_cache/tmp"
export TMPDIR="$PWD/apptainer_cache/tmp"
```

Container references are versioned tags or release URLs, not immutable image
digests. Reproducibility therefore also depends on the remote image references
remaining unchanged.

Configuration and lint checks for this revision used Nextflow 26.04.6. The
repository does not yet enforce minimum Nextflow, Slurm, or Apptainer versions.

## Containers

| Label | Container or execution mechanism |
| --- | --- |
| `prokka` | Galaxy Singularity Prokka 1.14.6 image |
| `bakta_prepare`, `bakta` | `oschwengers/bakta:v1.11.3` |
| `eggnog_prepare`, `eggnog_mapper_v2` | `brunoholiva/eggnog_v2:2.1.13-sse41` |
| `gffread` | Biocontainers GFFread 0.12.7 image |
| `kofamscan` | Biocontainers KofamScan 1.3.0 image |
| `interproscan_prepare` | `docker://python:3.11-bookworm` |
| `interproscan` | `interpro/interproscan:5.75-106.0` |
| `squashfs_tools` | `docker://merymardegan/squashfs-tools:4.6.1-ubuntu24.04` |
| `pgap` | `pgap.py` launches the configured PGAP SIF |

`kofamscan_prepare`, `pgap_python_prepare`, and `pgap_prepare` do not currently
declare a Nextflow container. Kofam preparation requires host `wget`, `gzip`,
`tar`, `find`, and `wc`; PGAP Python preparation requires `curl`, `tar`, and
bzip2 support; the PGAP wrapper requires a Singularity-compatible executable.

## Database lifecycle

Bakta, KofamScan, InterProScan, and eggNOG use local SquashFS database images.
Each resolver follows the same high-level contract:

1. validate an existing image;
2. prepare the configured database if validation fails;
3. package the prepared directory with `PACK_SQUASHFS`;
4. pass only the image path to the annotation process.

Managed Bakta, KofamScan, and InterProScan images include a `.size` manifest.
The resolver compares the image size with this manifest before reuse. The
database image is mounted read-only at the path expected by each tool:

| Database | Container mount |
| --- | --- |
| Bakta | `/database` |
| KofamScan | `/database` |
| eggNOG | `/database` |
| InterProScan | `/opt/interproscan/data` |

PGAP manages a directory, Python runtime, and SIF rather than a SquashFS data
image supplied as a process path input.

The eggNOG image is different: its `eggnog_2026-09` name is hardcoded, but the
downloader retrieves the current remote database without a release checksum.
An existing eggNOG image receives only basic extension, existence, and size
validation; its generated `.size` manifest is not consulted by the resolver.

Every `PACK_SQUASHFS` task requests 10 CPUs, 32 GB, and 24 hours.

## Parameter reference

| Parameter | Default or role |
| --- | --- |
| `annotation_type` | `both`; accepts `prokaryote`, `eukaryote`, or `both` |
| `prokaryote_samplesheet` | `assets/prokaryote_samplesheet.csv` |
| `eukaryote_samplesheet` | `assets/eukaryote_samplesheet.csv` |
| `outdir` | `data/reproduced` |
| `bakta_db` | optional external Bakta `.sqsh` image |
| `bakta_db_version` | `6.0` |
| `bakta_db_md5` | expected checksum for the configured Bakta archive |
| `bakta_amrfinder_db_series` | `4.0` |
| `bakta_amrfinder_db_version` | `2025-07-16.1` |
| `kofamscan_db` | optional external KofamScan `.sqsh` image |
| `kofamscan_db_release` | `2026-07-02` |
| `ips_db` | optional external InterProScan `.sqsh` image |
| `interproscan_version` | `5.75-106.0` data release and image name |
| `ips_formats` | `tsv,gff3` |
| `eggnog_db` | expected eggNOG image path |
| `pgap_version` | `2026-06-18.build8602` |
| `pgap_python_version` | `3.11` |
| `pgap_python_dir` | managed Python runtime directory |
| `pgap_dir` | managed PGAP installation directory |
| `pgap_container` | expected PGAP SIF path |

Version parameters are not independent compatibility switches. The Bakta URL,
checksum, AMRFinder release, and container must remain compatible;
`interproscan_version` does not change the hardcoded InterProScan container;
and PGAP runtime checks currently require Python 3.11. Override these values as
a coordinated code and provenance change, not as isolated command-line options.

## Directory ownership

| Path | Role |
| --- | --- |
| `work/` | task execution, logs, cache, and `-resume` source |
| `data/database/` | database images and PGAP installation |
| `data/runtime/` | managed runtimes such as PGAP Python |
| `data/reproduced/` | default annotation output root |
| `data/origin/` | archived source outputs used by the original analysis |

Do not remove `work/` while a resumed execution may still be needed. Published
files are copies for users and downstream analysis; Nextflow resumes from the
task outputs in `work/`.

## Output root

All annotation processes use `params.outdir`. Its default is:

```text
data/reproduced
```

Database and runtime preparation deliberately ignore `--outdir` because those
artifacts are reusable pipeline resources rather than annotation results.

## Known integration boundary

PANNZER is an external web-service result and has no Nextflow process. PGAP
publishes its native output directory, while the current prokaryotic notebooks
expect selected flat files. These boundaries are documented rather than hidden
by implicit file searches or fallback to archived outputs.
