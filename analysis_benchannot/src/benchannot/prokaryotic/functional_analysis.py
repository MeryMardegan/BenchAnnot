"""Module containing functional annotation text processing and categorization."""

from typing import List
import pandas as pd


def is_informative(product_text: str) -> bool:
    """Determine if a gene product annotation is informative.

    Parameters
    ----------
    product_text : str
        The annotation text for a gene product.

    Returns
    -------
    bool
        True if the annotation is informative, False if vague or missing.
    """
    if pd.isna(product_text):
        return False

    text = str(product_text).lower().strip()
    vague_terms = [
        "hypothetical",
        "uncharacterized",
        "unknown",
        "putative",
        "predicted",
        "duf",
        "domain of unknown function",
    ]

    if any(term in text for term in vague_terms) or text in {"", "-"}:
        return False

    return True


def categorize_functional_shifts(
    df_func: pd.DataFrame, tools: List[str], ref_source: str = "Reference"
) -> pd.DataFrame:
    """Compare tool annotations against a reference to categorize functional shifts.

    Parameters
    ----------
    df_func : pd.DataFrame
        DataFrame containing product annotations for the reference and tools.
    tools : list of str
        List of tool names to evaluate.
    ref_source : str, optional
        Column suffix for the ground truth reference, by default 'Reference'.

    Returns
    -------
    pd.DataFrame
        A summary table cross-tabulating tools by their functional shift category.
    """
    all_sources = tools + [ref_source]
    for source in all_sources:
        col_name = f"product_{source}"
        info_col = f"is_info_{source}"
        df_func[info_col] = df_func[col_name].apply(is_informative)

    comparison_results = []

    for tool in tools:
        tool_info = df_func[f"is_info_{tool}"]
        ref_info = df_func[f"is_info_{ref_source}"]

        both_info = tool_info & ref_info
        both_hypo = ~tool_info & ~ref_info
        over_annot = tool_info & ~ref_info
        under_annot = ~tool_info & ref_info

        categories = pd.Series("Unknown", index=df_func.index)
        categories[both_info] = "Both Informative"
        categories[both_hypo] = "Both Hypothetical"
        categories[over_annot] = "Over-Annotation"
        categories[under_annot] = "Under-Annotation"

        for category in categories:
            comparison_results.append({"Tool": tool, "Category": category})

    df_plot = pd.DataFrame(comparison_results)
    summary_table = df_plot.groupby(["Tool", "Category"]).size().unstack(fill_value=0)

    col_order = [
        "Both Informative",
        "Both Hypothetical",
        "Over-Annotation",
        "Under-Annotation",
    ]
    return summary_table[col_order]
