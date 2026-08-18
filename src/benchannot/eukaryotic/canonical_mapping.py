"""Conservative Swiss-Prot canonical mapping for eukaryotic references.

Canonical proteins are assigned externally from reviewed UniProtKB records.
An assignment requires an exact versioned RefSeq protein cross-reference and
an exact amino-acid sequence match to the NCBI reference FAA.
"""

from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Callable
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import matplotlib.pyplot as plt
import pandas as pd


_PROTEOME_PATTERN = re.compile(r"UP\d{9}")
_REFSEQ_PROTEIN_PATTERN = re.compile(r"^(?:AP|NP|TP|WP|XP|YP|ZP)_[A-Za-z0-9]+\.\d+$")
_CACHE_FIELDS = "accession,id,reviewed,sequence,cc_alternative_products,xref_refseq"


def uniprot_request_url(proteome_id: str) -> str:
    """Return the pinned-proteome REST URL used for immutable snapshots."""
    if not _PROTEOME_PATTERN.fullmatch(proteome_id):
        raise ValueError("proteome_id must be a versioned UniProt proteome ID")
    return "https://rest.uniprot.org/uniprotkb/stream?" + urlencode(
        {
            "format": "json",
            "fields": _CACHE_FIELDS,
            "includeIsoform": "true",
            "query": f"proteome:{proteome_id}",
        }
    )


def load_uniprot_cache(path: Path, expected_proteome_id: str | None = None) -> dict:
    """Load and cryptographically verify an immutable UniProt JSON cache."""
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not {"metadata", "response"} <= set(payload):
        raise ValueError(f"Invalid UniProt cache envelope: {path}")
    metadata = payload["metadata"]
    response = payload["response"]
    required = {
        "proteome_id",
        "request_url",
        "response_url",
        "retrieved_at",
        "response_sha256",
        "record_count",
    }
    if not isinstance(metadata, dict) or required - set(metadata):
        raise ValueError(f"Incomplete UniProt cache metadata: {path}")
    if expected_proteome_id and metadata["proteome_id"] != expected_proteome_id:
        raise ValueError(
            f"Cache proteome {metadata['proteome_id']} does not match {expected_proteome_id}"
        )
    records = _response_records(response)
    if metadata["record_count"] != len(records):
        raise ValueError("UniProt cache record_count does not match response")
    if metadata["response_sha256"] != _response_hash(response):
        raise ValueError("UniProt cache response SHA-256 does not match response")
    return payload


def ensure_uniprot_cache(
    proteome_id: str,
    destination: Path,
    *,
    import_snapshot: Path | None = None,
    opener: Callable = urlopen,
    timeout_seconds: float = 120.0,
) -> Path:
    """Use an existing cache, or immutably import/fetch it when absent.

    ``import_snapshot`` supports one-time migration of a prior verified cache.
    The resulting analysis depends only on ``destination``.
    """
    destination = Path(destination)
    if destination.exists():
        load_uniprot_cache(destination, proteome_id)
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    if import_snapshot is not None:
        payload = load_uniprot_cache(import_snapshot, proteome_id)
    else:
        request_url = uniprot_request_url(proteome_id)
        response = opener(
            Request(
                request_url,
                headers={"Accept": "application/json", "User-Agent": "BenchAnnot/1"},
            ),
            timeout=timeout_seconds,
        )
        try:
            body = response.read()
            response_url = response.geturl() if hasattr(response, "geturl") else request_url
        finally:
            if hasattr(response, "close"):
                response.close()
        decoded = json.loads(body)
        records = _response_records(decoded)
        payload = {
            "metadata": {
                "schema_version": "1",
                "proteome_id": proteome_id,
                "request_url": request_url,
                "response_url": response_url,
                "retrieved_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
                "response_sha256": _response_hash(decoded),
                "record_count": len(records),
            },
            "response": decoded,
        }
    with destination.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=True, indent=2, sort_keys=True)
        handle.write("\n")
    load_uniprot_cache(destination, proteome_id)
    return destination


def read_reference_faa(path: Path) -> dict[str, str]:
    """Read a reference FAA keyed by exact versioned header token."""
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    with Path(path).open(encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split(None, 1)[0]
                if current in sequences:
                    raise ValueError(f"{path}:{line_number} duplicates FAA ID {current}")
                sequences[current] = []
            elif current is None:
                raise ValueError(f"{path}:{line_number} has sequence before a header")
            else:
                sequences[current].append(line.upper())
    return {identifier: "".join(parts) for identifier, parts in sequences.items()}


def read_fasta_ids(path: Path) -> set[str]:
    """Return exact first-token identifiers from a FASTA file."""
    identifiers: set[str] = set()
    with Path(path).open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.add(line[1:].split(None, 1)[0].strip())
    return identifiers


def resolve_canonical_loci(
    reference: pd.DataFrame,
    faa_sequences: dict[str, str],
    cache_payload: dict,
    *,
    submitted_rna_ids: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resolve one unique reviewed canonical protein sequence per locus.

    Returns a transcript-level evidence mapping and one selection row for every
    uniquely resolved locus. Reference row order breaks equivalent-transcript
    ties and no unresolved locus receives a fallback representative.
    """
    required = {"RNA_ID", "protein_id", "locus_tag"}
    missing = required - set(reference)
    if missing:
        raise ValueError(f"Reference is missing columns: {sorted(missing)}")
    if reference["RNA_ID"].isna().any() or reference["RNA_ID"].duplicated().any():
        raise ValueError("Reference RNA_ID values must be present and unique")

    candidates_by_refseq = _reviewed_candidates_by_refseq(cache_payload)
    mapping_rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    submitted = submitted_rna_ids if submitted_rna_ids is not None else set(reference["RNA_ID"])

    for locus_tag, locus in reference.groupby("locus_tag", sort=False):
        evidence: list[dict[str, object]] = []
        exact_sequences: dict[str, list[dict[str, object]]] = {}
        has_reviewed_mapping = False
        for reference_order, (_, row) in enumerate(locus.iterrows()):
            protein_id = str(row["protein_id"])
            refseq_sequence = faa_sequences.get(protein_id)
            candidates = candidates_by_refseq.get(protein_id, ())
            has_reviewed_mapping |= bool(candidates)
            exact = [
                candidate
                for candidate in candidates
                if refseq_sequence is not None and candidate["sequence"] == refseq_sequence
            ]
            item = {
                "RNA_ID": row["RNA_ID"],
                "protein_id": protein_id,
                "locus_tag": locus_tag,
                "reference_order": reference_order,
                "reviewed_candidate_accessions": ";".join(
                    candidate["accession"] for candidate in candidates
                ),
                "exact_candidate_accessions": ";".join(
                    candidate["accession"] for candidate in exact
                ),
                "exact_sequence_match": bool(exact),
            }
            evidence.append(item)
            for candidate in exact:
                exact_sequences.setdefault(candidate["sequence"], []).append(
                    {**item, **candidate}
                )

        if len(exact_sequences) > 1:
            status = "ambiguous_reviewed_isoforms"
        elif len(exact_sequences) == 1:
            matching = next(iter(exact_sequences.values()))
            matching_rna = list(dict.fromkeys(item["RNA_ID"] for item in matching))
            status = (
                "equivalent_transcript_tie" if len(matching_rna) > 1 else "resolved_canonical"
            )
            selected_rna = matching_rna[0]
            selected = next(item for item in matching if item["RNA_ID"] == selected_rna)
            selection_rows.append(
                {
                    "locus_tag": locus_tag,
                    "RNA_ID": selected_rna,
                    "protein_id": selected["protein_id"],
                    "canonical_uniprot_accession": selected["accession"],
                    "canonical_sequence_sha256": sha256(
                        selected["sequence"].encode("ascii")
                    ).hexdigest(),
                    "resolution_status": status,
                    "equivalent_RNA_IDs": ";".join(matching_rna),
                    "submission_status": (
                        "submitted" if selected_rna in submitted else "canonical_not_submitted"
                    ),
                }
            )
        elif has_reviewed_mapping:
            status = "sequence_conflict"
        else:
            status = "no_reviewed_mapping"
        for item in evidence:
            mapping_rows.append({**item, "resolution_status": status})

    mapping = pd.DataFrame(mapping_rows)
    selection = pd.DataFrame(selection_rows)
    return mapping, selection


def build_canonical_functional_table(
    functional_table: pd.DataFrame, selection: pd.DataFrame
) -> pd.DataFrame:
    """Select exact canonical RNA rows; tool absence remains missing."""
    required = {"RNA_ID", "protein_id", "locus_tag"}
    if required - set(functional_table):
        raise ValueError("Functional table lacks RNA_ID, protein_id, or locus_tag")
    selected = selection.merge(
        functional_table,
        on=["RNA_ID", "protein_id", "locus_tag"],
        how="left",
        validate="one_to_one",
    )
    if len(selected) != len(selection):
        raise ValueError("Canonical selection did not preserve its locus denominator")
    return selected


def build_coverage_table(
    functional_table: pd.DataFrame,
    selection: pd.DataFrame,
    tools: list[str],
) -> pd.DataFrame:
    """Count RNA records, represented loci, and recovered canonical RNA per source."""
    canonical_ids = set(selection["RNA_ID"])
    rows = [
        {
            "Source": "NCBI RefSeq reference",
            "RNA_ID count": int(functional_table["RNA_ID"].nunique()),
            "Represented loci": int(functional_table["locus_tag"].nunique()),
            "Recovered canonical RNA_ID": len(canonical_ids),
        }
    ]
    for tool in tools:
        returned_column = f"{tool}_returned"
        if returned_column in functional_table:
            present = functional_table[returned_column].fillna(False).astype(bool)
        else:
            present = functional_table[tool].notna()
        represented = functional_table.loc[present]
        rows.append(
            {
                "Source": tool,
                "RNA_ID count": int(represented["RNA_ID"].nunique()),
                "Represented loci": int(represented["locus_tag"].nunique()),
                "Recovered canonical RNA_ID": int(represented["RNA_ID"].isin(canonical_ids).sum()),
            }
        )
    return pd.DataFrame(rows)


def plot_canonical_coverage(
    coverage: pd.DataFrame, title: str, output_path: Path
) -> None:
    """Save a 300-dpi grouped canonical-coverage bar chart."""
    metrics = ["RNA_ID count", "Represented loci", "Recovered canonical RNA_ID"]
    colors = ["#264653", "#2a9d8f", "#e76f51"]
    figure, axis = plt.subplots(figsize=(10, 6))
    x = list(range(len(coverage)))
    width = 0.24
    for offset, (metric, color) in enumerate(zip(metrics, colors)):
        bars = axis.bar(
            [value + (offset - 1) * width for value in x],
            coverage[metric],
            width,
            label=metric,
            color=color,
        )
        axis.bar_label(bars, labels=[f"{int(value):,}" for value in coverage[metric]], fontsize=7)
    axis.set_xticks(x, coverage["Source"])
    axis.set_ylabel("Exact count")
    axis.set_title(title, fontsize=14)
    axis.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    for spine in axis.spines.values():
        spine.set_visible(False)
    figure.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()
    plt.close(figure)


def _reviewed_candidates_by_refseq(payload: dict) -> dict[str, tuple[dict, ...]]:
    grouped: dict[str, list[dict]] = {}
    for record in _response_records(payload["response"]):
        entry_type = str(record.get("entryType", "")).lower()
        if "reviewed" not in entry_type or "unreviewed" in entry_type:
            continue
        accession = record.get("primaryAccession")
        sequence_value = record.get("sequence")
        sequence = sequence_value.get("value") if isinstance(sequence_value, dict) else sequence_value
        if not isinstance(accession, str) or not isinstance(sequence, str):
            continue
        is_isoform = bool(re.search(r"-\d+$", accession))
        if is_isoform and record.get("isoformSequenceStatus") != "Displayed":
            continue
        refseq_ids: list[str] = []
        for cross_reference in record.get("uniProtKBCrossReferences", []):
            if not isinstance(cross_reference, dict) or cross_reference.get("database") != "RefSeq":
                continue
            refseq_id = cross_reference.get("id")
            if isinstance(refseq_id, str) and _REFSEQ_PROTEIN_PATTERN.fullmatch(refseq_id):
                refseq_ids.append(refseq_id)
        candidate = {
            "accession": accession,
            "sequence": "".join(sequence.split()).upper(),
            "sequence_kind": "displayed_isoform" if is_isoform else "primary",
        }
        for refseq_id in dict.fromkeys(refseq_ids):
            grouped.setdefault(refseq_id, []).append(candidate)
    return {key: tuple(value) for key, value in grouped.items()}


def _response_records(response: object) -> list[dict]:
    if not isinstance(response, dict) or not isinstance(response.get("results"), list):
        raise ValueError("UniProt response must contain a results array")
    records = response["results"]
    if any(not isinstance(record, dict) for record in records):
        raise ValueError("UniProt response records must be objects")
    return records


def _response_hash(response: object) -> str:
    encoded = json.dumps(
        response, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("ascii")
    return sha256(encoded).hexdigest()
