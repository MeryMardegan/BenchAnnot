"""Module containing parsing logic for prokaryotic gene prediction tools"""

from typing import List, Optional, Tuple

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet

from benchannot.io.gff import parse_gff3


def clean_prokka_gff(input_path: str, output_path: str) -> str:
    """Remove the ##FASTA section from a Prokka GFF file.

    Parameters
    ----------
    input_path : str
        Path to the raw Prokka GFF3 file.
    output_path : str
        Path where the cleaned GFF3 file will be saved.

    Returns
    -------
    str
        The path to the cleaned GFF3 file.
    """
    with open(input_path, "r") as file, open(output_path, "w") as outfile:
        for line in file:
            if line.startswith("##FASTA"):
                break
            outfile.write(line)
    return output_path


def load_gff(file_path: str, source_name: str) -> pd.DataFrame:
    """Load a GFF3 file and add a Source column.

    Parameters
    ----------
    file_path : str
        Path to the GFF3 file.
    source_name : str
        Name of the annotation tool.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the GFF3 data with a Source column.
    """
    gff_df = parse_gff3(file_path, parse_attributes=True)
    gff_df["Source"] = source_name
    return gff_df


def standardize_gff(gff_df: pd.DataFrame, is_eggnog: bool = False) -> pd.DataFrame:
    """Standardize GFF3 files from different tools to have consistent column names and formats.

    Parameters
    ----------
    gff_df : pd.DataFrame
        DataFrame containing the GFF3 data with Source column already set.
    is_eggnog : bool, optional
        Whether the GFF is from eggNOG-mapper (which uses custom column names), by default False.

    Returns
    -------
    pd.DataFrame
        DataFrame containing only CDS features with standardized columns and a unified ORF ID.
    """
    if is_eggnog:
        gff_df = gff_df.rename(
            columns={"em_desc": "product", "em_Preferred_name": "gene"}
        )

    standard_df = gff_df[gff_df["Type"] == "CDS"].copy()
    standard_df["orf_id"] = (
        standard_df["Start"].astype(str)
        + "_"
        + standard_df["End"].astype(str)
        + "_"
        + standard_df["Strand"].astype(str)
    )
    return standard_df


def _calculate_f1_metrics(tp: int, fp: int, fn: int) -> Tuple[float, float, float]:
    """Helper function to calculate precision, recall, and F1-score.

    Parameters
    ----------
    tp : int
        Number of True Positives.
    fp : int
        Number of False Positives.
    fn : int
        Number of False Negatives.

    Returns
    -------
    tuple of float
        A tuple containing (precision, recall, f1_score).
    """
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1_score = (
        2 * (precision * recall) / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )

    return precision, recall, f1_score


def calculate_structural_metrics(
    all_gffs_df: pd.DataFrame, tools: List[str], ref_source: str = "Reference"
) -> pd.DataFrame:
    """Calculate True Positives, False Positives, False Negatives, Precision, Recall, and F1-Score.

    Parameters
    ----------
    all_gffs_df : pd.DataFrame
        Combined DataFrame containing standardized GFF data from all tools and the reference.
    tools : list of str
        List of annotation tool names to evaluate against the reference.
    ref_source : str, optional
        The source name of the ground truth reference annotations, by default 'Reference'.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the calculated structural metrics sorted by F1-Score in descending order.
    """
    reference_df = all_gffs_df[all_gffs_df["Source"] == ref_source]
    reference_orfs = set(reference_df["orf_id"])

    metrics_data = []

    for tool in tools:
        tool_df = all_gffs_df[all_gffs_df["Source"] == tool]
        tool_orfs = set(tool_df["orf_id"])

        tp = len(tool_orfs.intersection(reference_orfs))
        fp = len(tool_orfs - reference_orfs)
        fn = len(reference_orfs - tool_orfs)

        precision, recall, f1_score = _calculate_f1_metrics(tp, fp, fn)

        metrics_data.append(
            {
                "Annotation Tool": tool,
                "Exact Matches (TP)": tp,
                "Spurious/Shifted (FP)": fp,
                "Missed Genes (FN)": fn,
                "Precision": round(precision, 4),
                "Recall": round(recall, 4),
                "F1-Score": round(f1_score, 4),
            }
        )

    df_metrics = pd.DataFrame(metrics_data)
    return df_metrics.sort_values(by="F1-Score", ascending=False).reset_index(drop=True)


def plot_concordance_upset(
    all_gffs_df: pd.DataFrame,
    tools_order: List[str],
    title: str,
    save_path: Optional[str] = None,
) -> None:
    """Generate and optionally save an UpSet plot for gene prediction concordance.

    Parameters
    ----------
    all_gffs_df : pd.DataFrame
        Combined DataFrame containing standardized GFF data from all tools.
    tools_order : list of str
        List of tool names defining the order of categories in the UpSet plot.
    title : str
        Title of the UpSet plot.
    save_path : str, optional
        File path to save the generated figure. If None, the plot is not saved, by default None.
    """
    presence_matrix = pd.crosstab(all_gffs_df["orf_id"], all_gffs_df["Source"]) > 0
    presence_matrix = presence_matrix[tools_order].set_index(tools_order)

    fig = plt.figure(figsize=(10, 6))
    upset = UpSet(
        presence_matrix,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        facecolor="black",
    )

    upset.plot(fig=fig)
    plt.suptitle(title, fontsize=14, y=1.05)

    if save_path:
        fig.savefig(save_path, bbox_inches="tight", dpi=600)

    plt.show()


def export_concordant_cds(
    all_gffs_df: pd.DataFrame,
    tools: List[str],
    ref_source: str,
    output_path: str,
    concordant_only: bool = True,
) -> pd.DataFrame:
    """Export CDS with product annotations for all tools and reference.

    Parameters
    ----------
    all_gffs_df : pd.DataFrame
        Combined DataFrame containing standardized GFF data from all tools.
    tools : list of str
        List of annotation tool names to extract product data for.
    ref_source : str
        The source name of the ground truth reference annotations.
    output_path : str
        File path where the resulting CSV will be saved.
    concordant_only : bool, optional
        If True, only export CDS present in all tools and the reference.
        If False, export all reference CDS with left-joined tool annotations.
        By default True.

    Returns
    -------
    pd.DataFrame
        The merged DataFrame containing the CDS data with product annotations.
    """
    reference_df = all_gffs_df[all_gffs_df["Source"] == ref_source]

    if concordant_only:
        ref_orfs = set(reference_df["orf_id"])
        concordant_orfs = ref_orfs
        for tool in tools:
            tool_orfs = set(all_gffs_df[all_gffs_df["Source"] == tool]["orf_id"])
            concordant_orfs = concordant_orfs.intersection(tool_orfs)

        tp_df = reference_df[reference_df["orf_id"].isin(concordant_orfs)][
            ["orf_id", "Start", "End", "Strand"]
        ].copy()
    else:
        tp_df = reference_df[["orf_id", "Start", "End", "Strand"]].copy()

    for tool in tools:
        tool_sub = all_gffs_df[all_gffs_df["Source"] == tool][["orf_id", "product", "gene"]]
        tool_sub = tool_sub.rename(columns={"product": f"product_{tool}", "gene": f"gene_{tool}"})
        tp_df = pd.merge(tp_df, tool_sub, on="orf_id", how="left")

    ref_sub = reference_df[["orf_id", "product", "gene"]].rename(
        columns={"product": f"product_{ref_source}", "gene": f"gene_{ref_source}"}
    )
    tp_df = pd.merge(tp_df, ref_sub, on="orf_id", how="left")

    tp_df.to_csv(output_path, index=False)
    return tp_df