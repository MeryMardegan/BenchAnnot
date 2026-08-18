"""Module containing Biopython GBK parsing and UniProt REST API mapping logic."""

import io
import time
import requests
import pandas as pd
import numpy as np
from Bio import SeqIO
from typing import List, Dict

from benchannot.io.gff import parse_gff3


def clean_id(val: str) -> str:
    """Strip generic prefixes from IDs to ensure clean API mapping."""
    return str(val).replace("cds-", "").strip() if pd.notna(val) else np.nan


def get_best_mapping(candidates, mapping_dict):
    """Return the first mapped value from a list of candidates, or NaN."""
    if not isinstance(candidates, list):
        return np.nan
    for cand in candidates:
        if cand in mapping_dict:
            return mapping_dict[cand]
    return np.nan


def parse_eggnog_candidates(eggnog_gff_path: str) -> Dict[str, List[str]]:
    """Parse eggNOG GFF and return {orf_id: [candidate_ids]} mapping."""
    eggnog_df = parse_gff3(eggnog_gff_path, parse_attributes=True)
    eggnog_df["orf_id"] = eggnog_df.apply(
        lambda row: f"{row['Start']}_{row['End']}_{row['Strand']}", axis=1
    )

    eggnog_ids_dict = {}
    for _, row in eggnog_df.iterrows():
        target = row.get("em_target", np.nan)
        pref = row.get("em_Preferred_name", np.nan)
        eggnog_ids_dict[row["orf_id"]] = [
            x for x in [target, pref] if pd.notna(x) and str(x).strip() != ""
        ]

    return eggnog_ids_dict


def parse_gbk_identifiers(filepath: str, tool: str) -> pd.DataFrame:
    """Extract ORF IDs and tool-specific database cross-references from GenBank files."""
    extracted_ids = []
    for record in SeqIO.parse(filepath, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"
                orf_id = f"{start}_{end}_{strand}"
                extracted_val = None

                if tool == "Reference":
                    for xref in feature.qualifiers.get("db_xref", []):
                        if xref.startswith("UniProt"):
                            extracted_val = xref.split(":")[1]
                elif tool == "Prokka":
                    inferences = feature.qualifiers.get("inference", [None])
                    if inferences[0]:
                        for inf in inferences:
                            if "UniProtKB" in inf:
                                extracted_val = inf.split(":")[2]
                elif tool == "Bakta":
                    for note in feature.qualifiers.get("note", []):
                        if "UniRef:UniRef" in note:
                            extracted_val = note.split("_")[1]
                elif tool == "PGAP":
                    for inf in feature.qualifiers.get("inference", []):
                        if "RefSeq" in inf or "motif" in inf:
                            extracted_val = inf.split(":")[3]

                extracted_ids.append({"orf_id": orf_id, f"{tool}_ID": extracted_val})
    return pd.DataFrame(extracted_ids)


def safe_uniprot_mapping(
    id_list: List[str], from_db: str, to_db: str = "UniProtKB", tax_id: str = None
) -> pd.DataFrame:
    """API caller with chunking (1000 IDs) and strict error catching."""
    api_map = {
        "GI number": "GI_number",
        "EMBL/GenBank/DDBJ": "EMBL-GenBank-DDBJ",
        "Gene Name": "Gene_Name",
        "STRING": "STRING",
        "UniProtKB_AC-ID": "UniProtKB_AC-ID",
        "RefSeq_Protein": "RefSeq_Protein",
    }
    real_from = api_map.get(from_db, from_db)
    real_to = api_map.get(to_db, to_db)

    clean_ids = list(
        set([str(i).strip() for i in id_list if pd.notna(i) and str(i).strip() != ""])
    )
    if not clean_ids:
        return pd.DataFrame(columns=["From", "To"])

    print(f"  > Mapping: {from_db} -> {to_db} ({len(clean_ids)} IDs)...")
    submit_url = "https://rest.uniprot.org/idmapping/run"

    chunk_size = 1000
    all_mapped = []

    for i in range(0, len(clean_ids), chunk_size):
        chunk = clean_ids[i : i + chunk_size]
        payload = {"from": real_from, "to": real_to, "ids": ",".join(chunk)}

        if tax_id is not None and real_from in ["Gene_Name", "EMBL-GenBank-DDBJ"]:
            payload["taxId"] = str(tax_id)

        try:
            response = requests.post(submit_url, data=payload)
            response.raise_for_status()
            res_json = response.json()
            if "jobId" not in res_json:
                continue
            job_id = res_json["jobId"]
        except Exception as e:
            print(f"    [-] API Payload Exception ({real_from}): {str(e)}")
            continue

        status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"

        while True:
            try:
                status_res = requests.get(status_url)
                status_res.raise_for_status()
                status_data = status_res.json()
            except Exception:
                break

            if "jobStatus" in status_data:
                if status_data["jobStatus"] in ["RUNNING", "NEW"]:
                    time.sleep(2)
                else:
                    break
            else:
                break

        result_url = f"https://rest.uniprot.org/idmapping/stream/{job_id}?format=tsv"
        try:
            result_res = requests.get(result_url)
            result_res.raise_for_status()
            if result_res.text.strip():
                mapped = pd.read_csv(io.StringIO(result_res.text), sep="\t")
                all_mapped.append(mapped)
        except Exception:
            continue

    if all_mapped:
        final_df = pd.concat(all_mapped)
        print(f"    -> Success: Mapped {len(final_df)} IDs via {from_db}.")
        return final_df
    else:
        print(f"    -> Success: 0 IDs mapped via {from_db}.")
        return pd.DataFrame(columns=["From", "To"])


def smart_eggnog_mapping(id_list: List[str], tax_id: str = None) -> pd.DataFrame:
    """4-Tier Cascade Strategy for eggNOG output (STRING -> GI -> Gene Name -> EMBL)."""
    mapped_df = pd.DataFrame()
    unmapped_ids = id_list.copy()

    string_df = safe_uniprot_mapping(
        unmapped_ids, from_db="STRING", to_db="UniProtKB", tax_id=tax_id
    )
    if not string_df.empty and "From" in string_df.columns:
        mapped_df = pd.concat([mapped_df, string_df])
        unmapped_ids = [i for i in unmapped_ids if i not in string_df["From"].values]

    if not unmapped_ids:
        return mapped_df.drop_duplicates(subset=["From"])

    numeric_ids, alpha_ids, orig_to_clean = [], [], {}

    for full_id in unmapped_ids:
        clean_val = full_id.split(".", 1)[1] if "." in str(full_id) else full_id
        orig_to_clean[clean_val] = full_id
        if str(clean_val).isdigit():
            numeric_ids.append(clean_val)
        else:
            alpha_ids.append(clean_val)

    if numeric_ids:
        gi_df = safe_uniprot_mapping(
            numeric_ids, from_db="GI number", to_db="UniProtKB", tax_id=tax_id
        )
        if not gi_df.empty and "From" in gi_df.columns:
            gi_df["From"] = gi_df["From"].map(orig_to_clean)
            mapped_df = pd.concat([mapped_df, gi_df])
            numeric_mapped = [
                orig_to_clean.get(i, i) for i in gi_df["From"].dropna().values
            ]
            unmapped_ids = [i for i in unmapped_ids if i not in numeric_mapped]

    if alpha_ids:
        gene_df = safe_uniprot_mapping(
            alpha_ids, from_db="Gene Name", to_db="UniProtKB", tax_id=tax_id
        )
        if not gene_df.empty and "From" in gene_df.columns:
            gene_df["From"] = gene_df["From"].map(orig_to_clean)
            mapped_df = pd.concat([mapped_df, gene_df])
            alpha_mapped = [
                orig_to_clean.get(i, i) for i in gene_df["From"].dropna().values
            ]
            unmapped_ids = [i for i in unmapped_ids if i not in alpha_mapped]

    unmapped_alpha = [i for i in unmapped_ids if not str(i.split(".", 1)[-1]).isdigit()]
    if unmapped_alpha:
        clean_unmapped = [
            i.split(".", 1)[1] if "." in str(i) else i for i in unmapped_alpha
        ]
        embl_df = safe_uniprot_mapping(
            clean_unmapped,
            from_db="EMBL/GenBank/DDBJ",
            to_db="UniProtKB",
            tax_id=tax_id,
        )
        if not embl_df.empty and "From" in embl_df.columns:
            embl_df["From"] = embl_df["From"].map(orig_to_clean)
            mapped_df = pd.concat([mapped_df, embl_df])

    return mapped_df.drop_duplicates(subset=["From"])


def calculate_id_metrics(
    df_func: pd.DataFrame, tools: List[str], ref_source: str = "Reference"
) -> pd.DataFrame:
    """Calculate the Exact Match and Cluster Match rates against the reference."""
    results = []
    total_uniprot_ref = df_func[f"UniProt_{ref_source}"].notna().sum()
    total_cluster_ref = df_func[f"Cluster_{ref_source}"].notna().sum()

    for tool in tools:
        uniprot_matches = (
            df_func[f"UniProt_{tool}"] == df_func[f"UniProt_{ref_source}"]
        ).sum()
        uniprot_rate = (
            (uniprot_matches / total_uniprot_ref) * 100 if total_uniprot_ref > 0 else 0
        )

        cluster_matches = (
            df_func[f"Cluster_{tool}"] == df_func[f"Cluster_{ref_source}"]
        ).sum()
        cluster_rate = (
            (cluster_matches / total_cluster_ref) * 100 if total_cluster_ref > 0 else 0
        )

        results.append(
            {
                "Annotation Tool": tool,
                "Exact UniProt Rate (%)": round(uniprot_rate, 2),
                "UniRef100 Cluster Rate (%)": round(cluster_rate, 2),
                "Ref. IDs Available (100%)": total_uniprot_ref,
                "UniProt Matches (n)": uniprot_matches,
                "Cluster Matches (n)": cluster_matches,
            }
        )
    return pd.DataFrame(results)