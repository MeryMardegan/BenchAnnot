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
                            extracted_val = note.rsplit("_", 1)[-1]
                elif tool == "PGAP":
                    for inf in feature.qualifiers.get("inference", []):
                        if "RefSeq" in inf:
                            extracted_val = inf.split(":")[3]

                extracted_ids.append({"orf_id": orf_id, f"{tool}_ID": extracted_val})
    return pd.DataFrame(extracted_ids)


def safe_uniprot_mapping(
    id_list: List[str], from_db: str, to_db: str = "UniProtKB", tax_id: str = None
) -> pd.DataFrame:
    """API caller with chunking (1000 IDs) and strict error catching.

    Chunk composition is deterministic (input order preserved, no set()) and
    failed chunks are retried with rescue passes instead of being silently
    dropped, so results do not vary between runs.
    """
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
        dict.fromkeys(
            str(i).strip() for i in id_list if pd.notna(i) and str(i).strip() != ""
        )
    )
    if not clean_ids:
        return pd.DataFrame(columns=["From", "To"])

    print(f"  > Mapping: {from_db} -> {to_db} ({len(clean_ids)} IDs)...")
    submit_url = "https://rest.uniprot.org/idmapping/run"

    def map_chunk(chunk_ids: List[str]) -> pd.DataFrame:
        payload = {"from": real_from, "to": real_to, "ids": ",".join(chunk_ids)}
        if tax_id is not None and real_from in ["Gene_Name", "EMBL-GenBank-DDBJ"]:
            payload["taxId"] = str(tax_id)

        job_id = None
        for attempt in range(6):
            try:
                response = requests.post(submit_url, data=payload, timeout=60)
                response.raise_for_status()
                job_id = response.json().get("jobId")
                if job_id:
                    break
            except Exception as e:
                print(f"    [-] API Payload Exception ({real_from}): {str(e)}")
            time.sleep(2 * (attempt + 1))
        if not job_id:
            return None

        status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
        final_status = None
        for _ in range(60):
            try:
                status_res = requests.get(status_url, timeout=60)
                status_res.raise_for_status()
                status_data = status_res.json()
            except Exception:
                time.sleep(2)
                continue
            final_status = status_data.get("jobStatus")
            if final_status in ("RUNNING", "NEW"):
                time.sleep(2)
                continue
            break

        if final_status in ("ERROR", "FAILED", "SERVER_ERROR", "RUNNING", "NEW"):
            return None

        result_url = f"https://rest.uniprot.org/idmapping/stream/{job_id}?format=tsv"
        for attempt in range(4):
            try:
                result_res = requests.get(result_url, timeout=60)
                result_res.raise_for_status()
                if result_res.text.strip():
                    return pd.read_csv(io.StringIO(result_res.text), sep="\t")
                return pd.DataFrame(columns=["From", "To"])
            except Exception as e:
                print(f"    [-] Could not retrieve results ({real_from}): {str(e)}")
                time.sleep(2 * (attempt + 1))
        return None

    chunk_size = 1000
    all_mapped = []
    pending = []
    for i in range(0, len(clean_ids), chunk_size):
        chunk = clean_ids[i : i + chunk_size]
        mapped = map_chunk(chunk)
        if mapped is None:
            pending.append(chunk)
        else:
            all_mapped.append(mapped)
        time.sleep(1)

    for round_no in range(2):
        if not pending:
            break
        print(f"    [!] {len(pending)} chunk(s) failed, rescue pass {round_no + 1}/2...")
        remaining = []
        for chunk in pending:
            mapped = map_chunk(chunk)
            if mapped is None:
                remaining.append(chunk)
            else:
                all_mapped.append(mapped)
            time.sleep(2)
        pending = remaining

    if pending:
        lost = sum(len(c) for c in pending)
        sample = ", ".join(pending[0][:5])
        print(f"    [!!] Gave up on {len(pending)} chunk(s) ({lost} IDs): {sample} ...")

    if all_mapped:
        final_df = pd.concat(all_mapped)
        final_df = final_df.sort_values(["From", "To"]).reset_index(drop=True)
        print(f"    -> Success: Mapped {len(final_df)} IDs via {from_db}.")
        return final_df
    else:
        print(f"    -> Success: 0 IDs mapped via {from_db}.")
        return pd.DataFrame(columns=["From", "To"])


def build_candidate_map(mapping_df: pd.DataFrame) -> Dict[str, List[str]]:
    """Group a (From, To) mapping table into {source_id: [target_ids...]}.

    UniProt idmapping may return several target accessions for one source ID
    (e.g. strain-level RefSeq entries). Retaining them all avoids an arbitrary
    single-accession pick and makes downstream matching deterministic.
    """
    if mapping_df is None or mapping_df.empty or not {"From", "To"}.issubset(
        mapping_df.columns
    ):
        return {}
    valid = mapping_df[["From", "To"]].dropna()
    candidates = {}
    for source, group in valid.groupby("From", sort=False):
        targets = sorted({str(v).strip() for v in group["To"] if str(v).strip() != ""})
        if targets:
            candidates[str(source)] = targets
    return candidates


def smart_eggnog_mapping(id_list: List[str], tax_id: str = None) -> pd.DataFrame:
    """2-Tier Cascade Strategy for eggNOG output (STRING -> GI number).

    eggNOG em_target identifiers are locus tags (e.g. 224308.BSU00010), not
    EMBL/GenBank/DDBJ accessions. Gene Name mapping is excluded because it
    returns ambiguous cross-species hits. UniRef100 clustering downstream
    resolves strain-level differences.
    """
    mapped_df = pd.DataFrame()
    unmapped_ids = id_list.copy()

    string_df = safe_uniprot_mapping(
        unmapped_ids, from_db="STRING", to_db="UniProtKB", tax_id=tax_id
    )
    if not string_df.empty and "From" in string_df.columns:
        mapped_df = pd.concat([mapped_df, string_df])
        unmapped_ids = [i for i in unmapped_ids if i not in string_df["From"].values]

    if not unmapped_ids:
        return mapped_df.reset_index(drop=True)

    numeric_ids, orig_to_clean = [], {}

    for full_id in unmapped_ids:
        clean_val = full_id.split(".", 1)[1] if "." in str(full_id) else full_id
        orig_to_clean[clean_val] = full_id
        if str(clean_val).isdigit():
            numeric_ids.append(clean_val)

    if numeric_ids:
        gi_df = safe_uniprot_mapping(
            numeric_ids, from_db="GI number", to_db="UniProtKB", tax_id=tax_id
        )
        if not gi_df.empty and "From" in gi_df.columns:
            gi_df["From"] = gi_df["From"].map(orig_to_clean)
            mapped_df = pd.concat([mapped_df, gi_df])

    return mapped_df.reset_index(drop=True)


def _series_match(tool_vals: pd.Series, ref_vals: pd.Series) -> pd.Series:
    """Element-wise match: True where the reference value is in the tool's candidates.

    Supports both scalar tool columns (backward compatible) and list-valued
    candidate columns produced by build_candidate_map.
    """
    hits = []
    for t, r in zip(tool_vals, ref_vals):
        if pd.isna(r):
            hits.append(False)
        elif isinstance(t, list):
            hits.append(r in t)
        else:
            hits.append(False if pd.isna(t) else t == r)
    return pd.Series(hits, index=tool_vals.index)


def calculate_id_metrics(
    df_func: pd.DataFrame, tools: List[str], ref_source: str = "Reference"
) -> pd.DataFrame:
    """Calculate the Exact Match and Cluster Match rates against the reference.

    Tool UniProt/UniRef columns may hold a single accession or a list of
    candidate accessions; a gene matches when the reference value is among them.
    """
    results = []
    total_uniprot_ref = df_func[f"UniProt_{ref_source}"].notna().sum()
    total_cluster_ref = df_func[f"Cluster_{ref_source}"].notna().sum()

    for tool in tools:
        uniprot_matches = _series_match(
            df_func[f"UniProt_{tool}"], df_func[f"UniProt_{ref_source}"]
        ).sum()
        uniprot_rate = (
            (uniprot_matches / total_uniprot_ref) * 100 if total_uniprot_ref > 0 else 0
        )

        cluster_matches = _series_match(
            df_func[f"Cluster_{tool}"], df_func[f"Cluster_{ref_source}"]
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