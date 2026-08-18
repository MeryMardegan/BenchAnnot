# KofamScan Module

KofamScan is a tool for functional annotation that assigns KEGG Orthologs (KOs) to protein sequences using HMM profiles from the KEGG database.

## Input

- **Protein sequences** in FASTA format (`.faa`)  
  Typically generated from genome annotation (e.g., with Prokka or GFFREAD).
- **HMM profile database** (downloaded separately):  
  - `profiles/` → HMM profiles for KO assignments  
  - `ko_list` → KO ID definitions

In this pipeline, the database is expected to be located at:  
data/eukaryotes/db/

## Output

- Tab-delimited results with KO assignments for each input proteome.  
- Example output file: 
results/kofamscan/<sample>_ko.txt
The output contains:
- Query sequence ID
- Assigned KO number
- Threshold information
- E-value and score

## Running Standalone

To run KofamScan outside Nextflow (standalone example):

```bash
exec_annotation \
-o drosophila_ko.txt \
-p data/eukaryotes/db/profiles/ \
-k data/eukaryotes/db/ko_list \
drosophila_melanogaster.faa
```

## Integration in Pipeline

Module definition: modules/KofamScan.nf
Results are automatically stored under:
results/kofamscan/

## Database Setup

Before running, download the KEGG HMM profiles:
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz

tar -xzf profiles.tar.gz -C data/eukaryotes/db/
gunzip ko_list.gz -c > data/eukaryotes/db/ko_list

Ensure that:

profiles/ directory contains all .hmm files
ko_list file is in the same db/ folder

## Notes
- KofamScan uses adaptive score thresholds for higher accuracy compared to fixed E-values.
- The module is configured to run in a Docker container for reproducibility.
- Input proteins can come from Prokka or GFFREAD outputs.