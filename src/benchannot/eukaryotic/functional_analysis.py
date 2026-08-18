"""Functional annotation comparison utilities for eukaryotic benchmarks.

The informative/non-informative rules and shift categories intentionally match
the prokaryotic BenchAnnot analysis. Description similarity and canonical-
protein assignment are separate analyses and are not inferred by these functions.
"""

from pathlib import Path
import re
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import pandas as pd


VAGUE_TERMS = (
    "hypothetical",
    "uncharacterized",
    "unknown",
    "putative",
    "predicted",
    "duf",
    "domain of unknown function",
)

FUNCTIONAL_SHIFT_ORDER = (
    "Both Informative",
    "Both Hypothetical",
    "Over-Annotation",
    "Under-Annotation",
)

FUNCTIONAL_SHIFT_COLORS = {
    "Both Informative": "#2a9d8f",
    "Both Hypothetical": "#264653",
    "Over-Annotation": "#e76f51",
    "Under-Annotation": "#f4a261",
}

_ISOFORM_SUFFIX = re.compile(
    r"\s*,\s*(?:isoform|transcript variant)\s+"
    r"([A-Za-z0-9][A-Za-z0-9._-]*)\s*$",
    flags=re.IGNORECASE,
)


def is_informative(product_text: object) -> bool:
    """Return whether a product description contains functional information.

    This deliberately uses the same definition as the prokaryotic analysis:
    missing text and descriptions containing a vague term are non-informative.
    The function does not assess whether an informative description is correct.
    """
    if pd.isna(product_text):
        return False

    text = str(product_text).lower().strip()
    if text in {"", "-"}:
        return False
    return not any(term in text for term in VAGUE_TERMS)


def classify_description(product_text: object) -> str:
    """Classify a description as ``Informative`` or ``Non-informative``."""
    return "Informative" if is_informative(product_text) else "Non-informative"


def separate_isoform_suffix(product_text: object) -> tuple[object, object]:
    """Separate a terminal RefSeq isoform label from a product description.

    The original value is preserved for missing inputs. Only the explicit
    terminal ``, isoform X`` and ``, transcript variant X`` conventions are
    removed; occurrences of those terms elsewhere are left unchanged.
    """
    if pd.isna(product_text):
        return product_text, pd.NA

    text = str(product_text).strip()
    match = _ISOFORM_SUFFIX.search(text)
    if match is None:
        return text, pd.NA
    return text[: match.start()].strip(), match.group(1)


def prepare_reference_descriptions(reference: pd.DataFrame) -> pd.DataFrame:
    """Return a reference table with clean product text and isoform metadata.

    Required columns are ``RNA_ID``, ``protein_id``, ``locus_tag``, and
    ``Reference``.
    ``Reference_original`` retains the exact input text, while ``Reference`` is
    stripped of a terminal isoform label for downstream comparison.
    """
    required = {"RNA_ID", "protein_id", "locus_tag", "Reference"}
    missing = required - set(reference.columns)
    if missing:
        raise ValueError(f"Reference is missing columns: {sorted(missing)}")
    if reference["RNA_ID"].isna().any() or reference["RNA_ID"].duplicated().any():
        raise ValueError("Reference RNA_ID values must be present and unique")
    if reference["protein_id"].isna().any() or reference["protein_id"].eq("").any():
        raise ValueError("Reference protein_id values must be present")
    if reference["locus_tag"].isna().any() or reference["locus_tag"].eq("").any():
        raise ValueError("Reference locus_tag values must be present")

    result = reference.copy()
    result["Reference_original"] = result["Reference"]
    separated = result["Reference"].apply(separate_isoform_suffix)
    result["Reference"] = separated.str[0]
    result["Reference_isoform"] = separated.str[1]
    return result


def build_transcript_annotation_table(
    reference: pd.DataFrame,
    tool_tables: Mapping[str, pd.DataFrame],
) -> pd.DataFrame:
    """Left-join unique tool annotations to the RefSeq transcript universe.

    Each tool table must contain ``RNA_ID`` and exactly one description
    column. The one-to-one merge contract prevents accidental row expansion.
    """
    result = prepare_reference_descriptions(reference)
    reference_ids = set(result["RNA_ID"])

    for tool_name, table in tool_tables.items():
        description_columns = [column for column in table if column != "RNA_ID"]
        if len(description_columns) != 1:
            raise ValueError(
                f"{tool_name} must have RNA_ID and one description column"
            )
        if table["RNA_ID"].isna().any() or table["RNA_ID"].duplicated().any():
            raise ValueError(f"{tool_name} RNA_ID values must be present and unique")
        extra_ids = sorted(set(table["RNA_ID"]) - reference_ids)
        if extra_ids:
            raise ValueError(
                f"{tool_name} RNA_ID values absent from reference: {extra_ids[:5]}"
            )
        result = result.merge(
            table,
            on="RNA_ID",
            how="left",
            validate="one_to_one",
        )
    return result


def select_first_reference_transcript(transcript_table: pd.DataFrame) -> pd.DataFrame:
    """Select the first reference-ordered transcript for every locus.

    This deterministic reduction is provided as a sensitivity analysis. The
    all-transcript analysis remains the primary isoform-aware result.
    """
    required = {"RNA_ID", "locus_tag"}
    missing = required - set(transcript_table.columns)
    if missing:
        raise ValueError(f"Transcript table is missing columns: {sorted(missing)}")
    return transcript_table.drop_duplicates(subset="locus_tag", keep="first").reset_index(
        drop=True
    )


def categorize_functional_shifts(
    annotation_table: pd.DataFrame,
    tools: Sequence[str],
    reference_column: str = "Reference",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Categorize tool descriptions relative to reference informativeness.

    ``Over-Annotation`` means the tool is informative while the reference is
    non-informative. ``Under-Annotation`` is the inverse. ``Both Informative``
    does not imply that the two descriptions encode the same function.

    Returns a summary table and a detailed table with one row per input record
    and tool.
    """
    required = {"RNA_ID", "protein_id", "locus_tag", reference_column, *tools}
    missing = required - set(annotation_table.columns)
    if missing:
        raise ValueError(f"Annotation table is missing columns: {sorted(missing)}")

    reference_info = annotation_table[reference_column].apply(is_informative)
    details = []
    for tool in tools:
        tool_info = annotation_table[tool].apply(is_informative)
        categories = pd.Series(index=annotation_table.index, dtype="object")
        categories[tool_info & reference_info] = "Both Informative"
        categories[~tool_info & ~reference_info] = "Both Hypothetical"
        categories[tool_info & ~reference_info] = "Over-Annotation"
        categories[~tool_info & reference_info] = "Under-Annotation"

        details.append(
            pd.DataFrame(
                {
                    "RNA_ID": annotation_table["RNA_ID"],
                    "protein_id": annotation_table["protein_id"],
                    "locus_tag": annotation_table["locus_tag"],
                    "Tool": tool,
                    "Reference_description": annotation_table[reference_column],
                    "Tool_description": annotation_table[tool],
                    "Reference_informative": reference_info,
                    "Tool_informative": tool_info,
                    "Category": categories,
                }
            )
        )

    detailed = pd.concat(details, ignore_index=True)
    summary = (
        detailed.groupby(["Tool", "Category"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=FUNCTIONAL_SHIFT_ORDER, fill_value=0)
    )
    return summary, detailed


def add_reference_baseline(
    summary: pd.DataFrame,
    annotation_table: pd.DataFrame,
    reference_column: str = "Reference",
) -> pd.DataFrame:
    """Add the reference informative/non-informative base rates to a summary."""
    reference_info = annotation_table[reference_column].apply(is_informative)
    reference_row = pd.DataFrame(
        {
            "Both Informative": [int(reference_info.sum())],
            "Both Hypothetical": [int((~reference_info).sum())],
            "Over-Annotation": [0],
            "Under-Annotation": [0],
        },
        index=["NCBI RefSeq reference"],
    )
    return pd.concat([reference_row, summary]).reindex(columns=FUNCTIONAL_SHIFT_ORDER)


def plot_functional_shifts(
    summary_with_reference: pd.DataFrame,
    title: str,
    output_path: Path,
    count_label: str,
) -> None:
    """Save the prokaryotic-style diverging functional-shift chart as PNG."""
    totals = summary_with_reference.sum(axis=1)
    if totals.nunique() != 1:
        raise ValueError("Every source must be evaluated against the same universe")
    total = int(totals.iloc[0])
    plot_table = summary_with_reference.copy()
    plot_table["Over-Annotation"] *= -1
    plot_table["Under-Annotation"] *= -1

    figure, axis = plt.subplots(figsize=(10, 6))
    axis.barh(
        plot_table.index,
        plot_table["Over-Annotation"],
        color=FUNCTIONAL_SHIFT_COLORS["Over-Annotation"],
        label="Over-Annotation",
    )
    axis.barh(
        plot_table.index,
        plot_table["Under-Annotation"],
        left=plot_table["Over-Annotation"],
        color=FUNCTIONAL_SHIFT_COLORS["Under-Annotation"],
        label="Under-Annotation",
    )
    axis.barh(
        plot_table.index,
        plot_table["Both Informative"],
        color=FUNCTIONAL_SHIFT_COLORS["Both Informative"],
        label="Agreement (Informative)",
    )
    axis.barh(
        plot_table.index,
        plot_table["Both Hypothetical"],
        left=plot_table["Both Informative"],
        color=FUNCTIONAL_SHIFT_COLORS["Both Hypothetical"],
        label="Agreement (Non-informative)",
    )
    axis.axvline(0, color="black", linewidth=1, linestyle="--", ymax=0.96, ymin=0.04)
    for spine in axis.spines.values():
        spine.set_visible(False)
    for container in axis.containers:
        labels = [f"{abs(value) / total * 100:.1f}%" if value else "" for value in container.datavalues]
        axis.bar_label(
            container,
            labels=labels,
            label_type="center",
            fontsize=7,
            fontweight="bold",
        )

    axis.set_xlabel(f"{count_label} Count")
    axis.set_ylabel("Annotation Tool")
    axis.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda value, _: f"{int(abs(value))}")
    )
    axis.set_title(title, fontsize=14)
    handles, labels = axis.get_legend_handles_labels()
    legend_order = [1, 0, 2, 3]
    axis.legend(
        [handles[index] for index in legend_order],
        [labels[index] for index in legend_order],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=4,
        frameon=False,
    )
    figure.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()
    plt.close(figure)
