import subprocess
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

GFF_COLUMNS = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
FEATURE_COLUMNS = [
    *GFF_COLUMNS,
    "feature_id",
    "RNA_ID",
    "gene_id",
    "locus_tag",
    "product",
]
REMOVED_COLUMNS = [*FEATURE_COLUMNS, "removal_reason"]
MRNA_COLUMNS = FEATURE_COLUMNS
REFERENCE_COLUMNS = ["protein_id", "RNA_ID", "gene_id", "locus_tag", "product"]
LOCUS_COLUMNS = [
    "gene_id",
    "locus_tag",
    "transcript_count",
    "protein_count",
    "RNA_IDs",
    "protein_ids",
    "products",
]
PROJECT_ROOT = Path(__file__).resolve().parents[3]


def display_path(path):
    path = Path(path).resolve()
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)

def run_pre_processing(anno_sc, sample_id_sc, filtered_gff_sc):
    """
    Run the original Bash preprocessing step used for GFF filtering.
    """

    anno_sc = Path(anno_sc)
    filtered_gff_sc = Path(filtered_gff_sc)

    filtered_gff_sc.parent.mkdir(parents=True, exist_ok=True)

    bash_script = r'''
    set -euo pipefail

    anno_sc="$1"
    sample_id_sc="$2"
    filtered_gff_sc="$3"

    awk 'BEGIN{FS=OFS="\t"}
         /^#/ {print; next}
         $7=="?" {next}
         $9 ~ /exception=trans-splicing/ {next}
         {print}' "${anno_sc}" > "${filtered_gff_sc}"

    total=$(grep -vc '^#' "${anno_sc}" || true)
    kept=$(grep -vc '^#' "${filtered_gff_sc}" || true)
    removed=$(( total - kept ))

    echo "[GFFREAD - ${sample_id_sc}] total=$total kept=$kept removed=$removed"
    '''

    subprocess.run(
        [
            "bash",
            "-s",
            "--",
            str(anno_sc),
            str(sample_id_sc),
            str(filtered_gff_sc),
        ],
        input=bash_script,
        text=True,
        check=True,
    )

def locate_project_root():
    cwd = Path.cwd().resolve()
    for candidate in (cwd, *cwd.parents):
        if (candidate / "input.ipynb").is_file() and (candidate / "results").is_dir():
            return candidate
    raise FileNotFoundError("Could not locate the BenchAnnot project root")

def parse_attributes(raw_attributes):
    attributes = {}
    for item in raw_attributes.split(";"):
        key, separator, value = item.partition("=")
        if separator and key:
            attributes[unquote(key)] = unquote(value)
    return attributes

def make_gff_row(fields, attributes):
    row = dict(zip(GFF_COLUMNS, fields))
    feature_type = fields[2]
    parents = attributes.get("Parent", "").split(",")
    row.update(
        {
            "feature_id": attributes.get("ID", ""),
            "RNA_ID": attributes.get("ID", "") if feature_type == "mRNA" else "",
            "gene_id": (
                attributes.get("ID", "")
                if feature_type == "gene"
                else next((parent for parent in parents if parent.startswith("gene-")), "")
            ),
            "locus_tag": attributes.get("locus_tag", ""),
            "product": attributes.get("product", ""),
        }
    )
    return row

def scan_original_gff(path):
    counts = {"total_features": 0, "kept_features": 0, "total_mrna": 0, "kept_mrna": 0}
    removed_rows, locus_rows, protein_rows, kept_mrna_ids = [], [], [], []
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number} expected 9 GFF fields, found {len(fields)}")
            counts["total_features"] += 1
            attributes = parse_attributes(fields[8])
            reasons = []
            if fields[6] == "?":
                reasons.append("unknown_strand")
            if "exception=trans-splicing" in fields[8]:
                reasons.append("trans_splicing")
            row = make_gff_row(fields, attributes)
            if fields[2] == "mRNA":
                counts["total_mrna"] += 1
                if not row["RNA_ID"]:
                    raise ValueError(f"{path}:{line_number} mRNA has no ID attribute")
                locus_rows.append(
                    {
                        "RNA_ID": row["RNA_ID"],
                        "gene_id": row["gene_id"],
                        "locus_tag": row["locus_tag"],
                    }
                )
                if not reasons:
                    counts["kept_mrna"] += 1
                    kept_mrna_ids.append(row["RNA_ID"])
            elif fields[2] == "CDS" and attributes.get("protein_id"):
                parents = attributes.get("Parent", "").split(",")
                for parent in filter(None, parents):
                    protein_rows.append(
                        {
                            "protein_id": attributes["protein_id"],
                            "RNA_ID": parent,
                            "locus_tag": attributes.get("locus_tag", ""),
                        }
                    )
            if reasons:
                row["removal_reason"] = ",".join(reasons)
                removed_rows.append(row)
            else:
                counts["kept_features"] += 1
    removed = pd.DataFrame(removed_rows, columns=REMOVED_COLUMNS)
    locus_map = pd.DataFrame(locus_rows, columns=["RNA_ID", "gene_id", "locus_tag"])
    if locus_map["RNA_ID"].duplicated().any():
        raise ValueError(f"{path} contains duplicate mRNA IDs")
    protein_map = pd.DataFrame(
        protein_rows, columns=["protein_id", "RNA_ID", "locus_tag"]
    ).drop_duplicates()
    conflicting = protein_map[protein_map["protein_id"].duplicated(keep=False)]
    if not conflicting.empty:
        raise ValueError(
            f"{path} maps protein IDs to multiple transcripts: "
            f"{conflicting['protein_id'].drop_duplicates().head().tolist()}"
        )
    protein_map = protein_map.merge(
        locus_map.rename(columns={"locus_tag": "mRNA_locus_tag"}),
        on="RNA_ID",
        how="left",
        validate="many_to_one",
    )
    protein_map["locus_tag"] = protein_map["locus_tag"].where(
        protein_map["locus_tag"].ne(""), protein_map["mRNA_locus_tag"]
    )
    protein_map = protein_map.drop(columns="mRNA_locus_tag")
    if not removed.empty:
        transcript_genes = locus_map.set_index("RNA_ID")["gene_id"].to_dict()
        removed["gene_id"] = removed.apply(
            lambda row: row["gene_id"]
            or next(
                (
                    transcript_genes[parent]
                    for parent in parse_attributes(row["attributes"])
                    .get("Parent", "")
                    .split(",")
                    if parent in transcript_genes
                ),
                "",
            ),
            axis=1,
        )
    return {
        **counts,
        "kept_mrna_ids": tuple(kept_mrna_ids),
        "removed": removed,
        "locus_map": locus_map,
        "protein_map": protein_map,
    }

def read_filtered_mrna(path):
    rows = []
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number} expected 9 GFF fields, found {len(fields)}")
            if fields[2] == "mRNA":
                rows.append(make_gff_row(fields, parse_attributes(fields[8])))
    frame = pd.DataFrame(rows, columns=MRNA_COLUMNS)
    if frame["RNA_ID"].duplicated().any():
        raise ValueError(f"{path} contains duplicate mRNA IDs")
    return frame

def read_reference_fasta(path):
    records = []
    current_id = None
    sequence_length = 0
    with path.open(encoding="ascii", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                if current_id is not None and sequence_length == 0:
                    raise ValueError(f"{path} contains an empty sequence for {current_id}")
                header = line[1:].strip()
                identifier, separator, description = header.partition(" ")
                if not identifier:
                    raise ValueError(f"{path}:{line_number} has an empty FASTA ID")
                if description.endswith("]") and " [" in description:
                    description = description.rsplit(" [", 1)[0]
                records.append({"protein_id": identifier, "product": description})
                current_id = identifier
                sequence_length = 0
            elif line.strip():
                if current_id is None:
                    raise ValueError(f"{path}:{line_number} has sequence before a header")
                sequence_length += len(line.strip())
    if current_id is not None and sequence_length == 0:
        raise ValueError(f"{path} contains an empty sequence for {current_id}")
    reference = pd.DataFrame(records, columns=["protein_id", "product"])
    if reference["protein_id"].duplicated().any():
        raise ValueError(f"{path} contains duplicate FASTA IDs")
    if reference.empty:
        raise ValueError(f"{path} contains no FASTA records")
    return reference


def read_fasta_ids(path):
    return tuple(read_reference_fasta(path)["protein_id"])

def write_workbook(path, tables):
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, table in tables.items():
            table.to_excel(writer, sheet_name=sheet_name, index=False)
            worksheet = writer.sheets[sheet_name]
            worksheet.freeze_panes = "A2"
            worksheet.auto_filter.ref = worksheet.dimensions
            for cell in worksheet[1]:
                cell.font = Font(bold=True)
            for column_index, column_name in enumerate(table.columns, 1):
                values = table[column_name].astype(str).head(1000)
                width = max([len(str(column_name)), *(len(value) for value in values)], default=10)
                worksheet.column_dimensions[get_column_letter(column_index)].width = min(width + 2, 60)


def _summary_row(stage, metric, value, definition):
    return {"stage": stage, "metric": metric, "value": value, "definition": definition}


def _join_unique(values):
    return " | ".join(dict.fromkeys(value for value in values if value))


def build_locus_reference(transcript_reference):
    """Collapse the transcript reference to exactly one row per annotated locus."""

    rows = []
    for (gene_id, locus_tag), group in transcript_reference.groupby(
        ["gene_id", "locus_tag"], sort=False, dropna=False
    ):
        rows.append(
            {
                "gene_id": gene_id,
                "locus_tag": locus_tag,
                "transcript_count": group["RNA_ID"].nunique(),
                "protein_count": group["protein_id"].nunique(),
                "RNA_IDs": _join_unique(group["RNA_ID"]),
                "protein_ids": _join_unique(group["protein_id"]),
                "products": _join_unique(group["product"]),
            }
        )
    return pd.DataFrame(rows, columns=LOCUS_COLUMNS)

def build_reference_audit(
    organism,
    original_gff,
    reference_fasta,
    filtered_gff,
    gffread_fasta,
    mitochondrial_prefix,
    output_directory,
):
    """Build separate GFF/GFFread and NCBI-reference audit reports."""

    original_gff = Path(original_gff)
    reference_fasta = Path(reference_fasta)
    filtered_gff = Path(filtered_gff)
    gffread_fasta = Path(gffread_fasta)
    output_directory = Path(output_directory)
    scan = scan_original_gff(original_gff)
    filtered_mrna = read_filtered_mrna(filtered_gff)
    filtered_ids = tuple(filtered_mrna["RNA_ID"])
    if filtered_ids != scan["kept_mrna_ids"]:
        raise ValueError(f"{organism}: filtered GFF differs from the reconstructed pre-GFFread filter")
    gffread_ids = read_fasta_ids(gffread_fasta)
    gffread_set, filtered_set = set(gffread_ids), set(filtered_ids)
    extra_gffread = sorted(gffread_set - filtered_set)
    if extra_gffread:
        raise ValueError(f"{organism}: GFFread IDs absent from filtered GFF: {extra_gffread[:5]}")
    post_missing = filtered_mrna[~filtered_mrna["RNA_ID"].isin(gffread_set)].reset_index(drop=True)
    unexpected_gffread = pd.DataFrame({"RNA_ID": extra_gffread}, columns=["RNA_ID"])
    reference_proteins = read_reference_fasta(reference_fasta)
    reference = reference_proteins.merge(
        scan["protein_map"], on="protein_id", how="left", validate="one_to_one"
    )[REFERENCE_COLUMNS]
    if reference[["RNA_ID", "locus_tag"]].isna().any().any() or reference[
        ["RNA_ID", "locus_tag"]
    ].eq("").any().any():
        raise ValueError(f"{organism}: NCBI proteins could not all be mapped through the GFF")
    if reference["RNA_ID"].duplicated().any():
        raise ValueError(f"{organism}: reference contains duplicate RNA_ID values")
    mitochondrial_mask = reference["RNA_ID"].str.startswith(mitochondrial_prefix)
    mitochondrial_removed = reference[mitochondrial_mask].reset_index(drop=True)
    reference_without_mitochondria = reference[~mitochondrial_mask].reset_index(drop=True)
    transcript_reference = reference_without_mitochondria.copy()
    locus_reference = build_locus_reference(transcript_reference)
    if locus_reference["locus_tag"].duplicated().any():
        raise ValueError(f"{organism}: locus_tag maps to more than one gene_id")
    reference_rna_set = set(transcript_reference["RNA_ID"])
    gffread_not_in_reference = sorted(gffread_set - reference_rna_set)
    gffread_not_in_reference_table = pd.DataFrame(
        {"RNA_ID": gffread_not_in_reference}, columns=["RNA_ID"]
    )
    reference_not_in_gffread = transcript_reference[
        ~transcript_reference["RNA_ID"].isin(gffread_set)
    ].reset_index(drop=True)
    gff_summary = pd.DataFrame(
        [
            _summary_row("Input", "original_gff", display_path(original_gff), "Original annotation GFF"),
            _summary_row("Input", "filtered_gff", display_path(filtered_gff), "GFF supplied to GFFread"),
            _summary_row("Input", "gffread_fasta", display_path(gffread_fasta), "Protein FASTA emitted by GFFread"),
            _summary_row("Pre-GFFread", "features_total", scan["total_features"], "Non-comment GFF records before filtering"),
            _summary_row("Pre-GFFread", "features_kept", scan["kept_features"], "Records retained for GFFread"),
            _summary_row("Pre-GFFread", "features_removed", len(scan["removed"]), "Unknown-strand or trans-splicing records"),
            _summary_row("Pre-GFFread", "mRNA_total", scan["total_mrna"], "mRNA records before filtering"),
            _summary_row("Pre-GFFread", "mRNA_kept", scan["kept_mrna"], "mRNA records supplied to GFFread"),
            _summary_row("Pre-GFFread", "mRNA_removed", int((scan["removed"]["type"] == "mRNA").sum()), "mRNA records removed before GFFread"),
            _summary_row("Post-GFFread", "gffread_protein_IDs", len(gffread_ids), "Unique transcript IDs in the GFFread FAA"),
            _summary_row("Post-GFFread", "mRNA_missing_after_gffread", len(post_missing), "Filtered-GFF mRNAs absent from the GFFread FAA"),
            _summary_row("Post-GFFread", "unexpected_gffread_IDs", len(unexpected_gffread), "GFFread IDs absent from the filtered GFF"),
            _summary_row("NCBI comparison", "gffread_IDs_absent_from_reference", len(gffread_not_in_reference_table), "GFFread IDs absent from the mitochondrial-free NCBI universe"),
            _summary_row("NCBI comparison", "reference_IDs_absent_from_gffread", len(reference_not_in_gffread), "NCBI transcript IDs absent from GFFread; does not alter the reference"),
        ]
    )
    reference_summary = pd.DataFrame(
        [
            _summary_row("Input", "original_gff", display_path(original_gff), "GFF used for exact identifier mapping"),
            _summary_row("Input", "reference_fasta", display_path(reference_fasta), "NCBI FAA defining the reference universe"),
            _summary_row("Mapping", "reference_proteins", len(reference), "Proteins in the NCBI FAA mapped through the GFF"),
            _summary_row("Mitochondrial removal", "rows_removed", len(mitochondrial_removed), f"RNA_ID starts with {mitochondrial_prefix}"),
            _summary_row("Prepared reference", "transcript_rows", len(transcript_reference), "One row per protein/transcript"),
            _summary_row("Prepared reference", "unique_RNA_ID", transcript_reference["RNA_ID"].nunique(), "Exact transcript identifiers"),
            _summary_row("Prepared reference", "locus_rows", len(locus_reference), "Exactly one row per locus_tag"),
        ]
    )
    gff_readme = pd.DataFrame(
        [
            {"section": "Purpose", "topic": "Universe", "description": "GFFread is audited only; it never defines or filters the NCBI reference universe."},
            {"section": "Inputs", "topic": "Original GFF", "description": display_path(original_gff)},
            {"section": "Inputs", "topic": "Filtered GFF", "description": display_path(filtered_gff)},
            {"section": "Inputs", "topic": "GFFread FAA", "description": display_path(gffread_fasta)},
            {"section": "Identifiers", "topic": "feature_id", "description": "Exact GFF ID attribute for any feature type."},
            {"section": "Identifiers", "topic": "RNA_ID", "description": "Exact ID attribute only for mRNA rows; blank for every other feature type."},
            {"section": "Identifiers", "topic": "gene_id", "description": "Gene feature ID or mRNA gene Parent; child rows are resolved through their transcript Parent when possible."},
            {"section": "Paths", "topic": "Documentation", "description": "All input paths in this report are relative to the project root."},
            {"section": "Interpretation", "topic": "Stages", "description": "Pre-GFFread rows were not supplied to GFFread; post-GFFread rows were supplied but emitted no protein."},
        ]
    )
    reference_readme = pd.DataFrame(
        [
            {"section": "Purpose", "topic": "Universe", "description": "The NCBI protein FAA is the reference universe; GFFread does not filter it."},
            {"section": "Inputs", "topic": "Original GFF", "description": display_path(original_gff)},
            {"section": "Inputs", "topic": "NCBI reference FAA", "description": display_path(reference_fasta)},
            {"section": "Identifiers", "topic": "Transcript table", "description": "One row per NCBI protein and exact RNA_ID, with gene_id and locus_tag mapped through the GFF."},
            {"section": "Identifiers", "topic": "Locus table", "description": "Exactly one row per locus_tag; transcript, protein, and product values are aggregated without loss."},
            {"section": "Products", "topic": "Source", "description": "Products come from NCBI FAA descriptions; organism suffixes in square brackets are removed."},
            {"section": "Paths", "topic": "Documentation", "description": "All input paths in this report are relative to the project root."},
        ]
    )

    gff_directory = output_directory / "audit" / "gff_gffread" / organism
    reference_directory = output_directory / "audit" / "reference" / organism
    prepared_directory = reference_directory / "prepared_reference"
    gff_directory.mkdir(parents=True, exist_ok=True)
    prepared_directory.mkdir(parents=True, exist_ok=True)
    outputs = {
        "gff_workbook": gff_directory / "gff_gffread_audit.xlsx",
        "gff_summary": gff_directory / "summary.tsv",
        "pre_gffread_removed": gff_directory / "pre_gffread_removed.tsv",
        "post_gffread_missing": gff_directory / "post_gffread_missing.tsv",
        "unexpected_gffread": gff_directory / "unexpected_gffread_ids.tsv",
        "gffread_not_in_reference": gff_directory / "gffread_not_in_reference.tsv",
        "reference_not_in_gffread": gff_directory / "reference_not_in_gffread.tsv",
        "reference_workbook": reference_directory / "ncbi_reference_audit.xlsx",
        "reference_summary": reference_directory / "summary.tsv",
        "mitochondrial_removed": reference_directory / "mitochondrial_removed.tsv",
        "transcript_reference": prepared_directory / "transcript_reference.tsv",
        "locus_reference": prepared_directory / "locus_reference.tsv",
    }
    sidecars = {
        "gff_summary": gff_summary,
        "pre_gffread_removed": scan["removed"],
        "post_gffread_missing": post_missing,
        "unexpected_gffread": unexpected_gffread,
        "gffread_not_in_reference": gffread_not_in_reference_table,
        "reference_not_in_gffread": reference_not_in_gffread,
        "reference_summary": reference_summary,
        "mitochondrial_removed": mitochondrial_removed,
        "transcript_reference": transcript_reference,
        "locus_reference": locus_reference,
    }
    for name, table in sidecars.items():
        table.to_csv(outputs[name], sep="\t", index=False)
    write_workbook(
        outputs["gff_workbook"],
        {
            "README": gff_readme,
            "Summary": gff_summary,
            "Pre-GFFread removed": scan["removed"],
            "Post-GFFread missing": post_missing,
            "Unexpected GFFread IDs": unexpected_gffread,
            "GFFread absent from reference": gffread_not_in_reference_table,
            "Reference absent from GFFread": reference_not_in_gffread,
        },
    )
    write_workbook(
        outputs["reference_workbook"],
        {
            "README": reference_readme,
            "Summary": reference_summary,
            "Mitochondrial removal": mitochondrial_removed,
            "Transcript reference": transcript_reference,
            "Locus reference": locus_reference,
        },
    )
    return {
        "gff_readme": gff_readme,
        "reference_readme": reference_readme,
        "gff_summary": gff_summary,
        "reference_summary": reference_summary,
        "pre_gffread_removed": scan["removed"],
        "post_gffread_missing": post_missing,
        "reference_not_in_gffread": reference_not_in_gffread,
        "mitochondrial_removed": mitochondrial_removed,
        "locus_map": scan["locus_map"],
        "transcript_reference": transcript_reference,
        "locus_reference": locus_reference,
        # Retain return aliases for callers while all new files use explicit contracts.
        "final_reference": transcript_reference,
        "locus_products": locus_reference,
        "outputs": outputs,
    }
