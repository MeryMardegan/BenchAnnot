"""Utilities for loading and parsing GFF3 files into pandas DataFrames."""

from urllib.parse import unquote
import pandas as pd


def _parse_attributes(attributes: str) -> dict:
    """Parse a GFF3 attributes string (key=value;key2=value2)."""
    parsed = {}
    if not attributes or attributes == ".":
        return parsed

    for item in attributes.split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue

        key, value = item.split("=", 1)
        parsed[key] = unquote(value)

    return parsed


def parse_gff3(file_path: str, parse_attributes: bool = True) -> pd.DataFrame:
    """Load a GFF3 file and optionally expand attributes to individual columns."""
    column_names = [
        "Seqid",
        "Source",
        "Type",
        "Start",
        "End",
        "Score",
        "Strand",
        "Phase",
        "Attributes",
    ]

    gff_df = pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        header=None,
        names=column_names,
        dtype=str,
        na_filter=False,
    )

    if gff_df.empty:
        return gff_df

    gff_df["Start"] = pd.to_numeric(gff_df["Start"], errors="coerce")
    gff_df["End"] = pd.to_numeric(gff_df["End"], errors="coerce")

    if not parse_attributes:
        return gff_df

    attrs_df = gff_df["Attributes"].map(_parse_attributes).apply(pd.Series)
    attrs_df = attrs_df.replace("", pd.NA)

    return pd.concat([gff_df.drop(columns=["Attributes"]), attrs_df], axis=1)