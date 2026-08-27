import csv
import warnings
from pathlib import Path

import pandas as pd


KOFAM_COLUMNS = [
    "marker",
    "RNA_ID",
    "KO",
    "threshold",
    "score",
    "evalue",
    "Kofam",
]
INTERPRO_DATABASES = ("PANTHER", "NCBIfam", "SFLD", "Hamap", "PIRSF")
TOOL_ORDER = ("Kofam", "Pannzer", "EggNOG", "InterProScan")
SOURCE_ORDER = ("Reference", *TOOL_ORDER)


def _remove_mitochondrial_ids(
    table: pd.DataFrame, mitochondrial_prefix: str | None
) -> pd.DataFrame:
    if not mitochondrial_prefix:
        return table
    return table[~table["RNA_ID"].str.startswith(mitochondrial_prefix, na=False)]


def _validate_unique_rna_ids(table: pd.DataFrame, source: str) -> None:
    if "RNA_ID" not in table:
        raise ValueError(f"{source} is missing the RNA_ID column")
    if table["RNA_ID"].isna().any() or table["RNA_ID"].eq("").any():
        raise ValueError(f"{source} contains empty RNA_ID values")
    duplicated = table.loc[table["RNA_ID"].duplicated(), "RNA_ID"]
    if not duplicated.empty:
        raise ValueError(
            f"{source} contains duplicate RNA_ID values: {duplicated.head().tolist()}"
        )


def _set_processing_counts(
    result: pd.DataFrame, raw_rows: int, identified_rna_ids: set[str]
) -> pd.DataFrame:
    result.attrs["raw_rows"] = raw_rows
    result.attrs["identified_rna_ids"] = frozenset(identified_rna_ids)
    return result


def read_kofam(path: Path) -> pd.DataFrame:
    """Read KofamScan detail output and validate its seven-column schema."""
    path = Path(path)
    table = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        comment="#",
        header=None,
        engine="python",
        quotechar='"',
    )
    if table.shape[1] != len(KOFAM_COLUMNS):
        raise ValueError(
            f"{path} has {table.shape[1]} columns; expected {len(KOFAM_COLUMNS)}"
        )
    table.columns = KOFAM_COLUMNS
    if table["RNA_ID"].isna().any() or table["RNA_ID"].eq("").any():
        raise ValueError(f"{path} contains empty RNA_ID values")
    return _set_processing_counts(table, len(table), set(table["RNA_ID"]))


def read_fasta_ids(path: Path) -> set[str]:
    """Read unique first-token FASTA identifiers without normalizing them."""
    identifiers = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].strip().split(maxsplit=1)[0]
                if not identifier:
                    raise ValueError(f"{path} contains an empty FASTA identifier")
                identifiers.append(identifier)
    if not identifiers:
        raise ValueError(f"{path} contains no FASTA records")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"{path} contains duplicate FASTA identifiers")
    return set(identifiers)


def prepare_reference(path: Path) -> pd.DataFrame:
    """Read the prepared transcript reference without changing exact identifiers."""
    path = Path(path)
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"protein_id", "RNA_ID", "locus_tag", "product"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    result = table[["protein_id", "RNA_ID", "locus_tag", "product"]].rename(
        columns={"product": "Reference"}
    )
    _validate_unique_rna_ids(result, str(path))
    if result["protein_id"].eq("").any():
        raise ValueError(f"{path} contains empty protein_id values")
    if result["locus_tag"].eq("").any():
        raise ValueError(f"{path} contains empty locus_tag values")
    return result


def prepare_kofam(
    path: Path,
    significant_only: bool,
    mitochondrial_prefix: str | None = None,
) -> pd.DataFrame:
    """Return one Kofam row per RNA_ID, optionally retaining only marked hits."""
    table = read_kofam(path)
    raw_rows = len(table)
    identified_rna_ids = {
        value
        for value in table["RNA_ID"]
        if not mitochondrial_prefix or not value.startswith(mitochondrial_prefix)
    }
    if significant_only:
        table = table[table["marker"].eq("*")]
    result = (
        table.drop_duplicates(subset="RNA_ID", keep="first")
        .loc[:, ["RNA_ID", "Kofam"]]
        .copy()
    )
    result = _remove_mitochondrial_ids(result, mitochondrial_prefix).reset_index(
        drop=True
    )
    _validate_unique_rna_ids(result, str(path))
    return _set_processing_counts(result, raw_rows, identified_rna_ids)


def prepare_pannzer(
    path: Path, mitochondrial_prefix: str | None = None
) -> pd.DataFrame:
    """Select PANNZER description records (DE), one row per RNA_ID."""
    path = Path(path)
    table = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        comment="#",
        header=None,
        engine="python",
        quotechar='"',
    )
    if table.shape[1] != 6:
        raise ValueError(f"{path} has {table.shape[1]} columns; expected 6")
    if not table.empty and table.iloc[0, 0] == "qpid":
        table = table.iloc[1:].reset_index(drop=True)
    raw_rows = len(table)
    identified_rna_ids = {
        value
        for value in table[0]
        if not mitochondrial_prefix or not value.startswith(mitochondrial_prefix)
    }
    result = (
        table[table[1].eq("DE")]
        .drop_duplicates(subset=0, keep="first")[[0, 5]]
        .rename(columns={0: "RNA_ID", 5: "Pannzer"})
    )
    result = _remove_mitochondrial_ids(result, mitochondrial_prefix).reset_index(
        drop=True
    )
    _validate_unique_rna_ids(result, str(path))
    return _set_processing_counts(result, raw_rows, identified_rna_ids)


def prepare_eggnog(
    path: Path, mitochondrial_prefix: str | None = None
) -> pd.DataFrame:
    """Select the eggNOG description, one row per RNA_ID."""
    path = Path(path)
    table = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        comment="#",
        header=None,
        engine="python",
        quotechar='"',
    )
    if table.shape[1] < 8:
        raise ValueError(f"{path} has {table.shape[1]} columns; expected at least 8")
    raw_rows = len(table)
    identified_rna_ids = {
        value
        for value in table[0]
        if not mitochondrial_prefix or not value.startswith(mitochondrial_prefix)
    }
    result = (
        table.drop_duplicates(subset=0, keep="first")[[0, 7]]
        .rename(columns={0: "RNA_ID", 7: "EggNOG"})
        .copy()
    )
    result = _remove_mitochondrial_ids(result, mitochondrial_prefix).reset_index(
        drop=True
    )
    _validate_unique_rna_ids(result, str(path))
    return _set_processing_counts(result, raw_rows, identified_rna_ids)


def prepare_interproscan(
    path: Path, mitochondrial_prefix: str | None = None
) -> pd.DataFrame:
    """Select the lowest-e-value family assignment per RNA_ID."""
    path = Path(path)
    table = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        comment="#",
        header=None,
        engine="python",
        quoting=csv.QUOTE_NONE,
    )
    if table.shape[1] < 13:
        raise ValueError(f"{path} has {table.shape[1]} columns; expected at least 13")
    raw_rows = len(table)
    identified_rna_ids = {
        value
        for value in table[0]
        if not mitochondrial_prefix or not value.startswith(mitochondrial_prefix)
    }
    filtered = table[table[3].isin(INTERPRO_DATABASES)].copy()
    filtered["_evalue"] = pd.to_numeric(filtered[8], errors="coerce")
    filtered = filtered[filtered["_evalue"].notna()]
    result = (
        filtered.sort_values("_evalue", kind="stable")
        .drop_duplicates(subset=0, keep="first")[[0, 12]]
        .rename(columns={0: "RNA_ID", 12: "InterProScan"})
    )
    result = _remove_mitochondrial_ids(result, mitochondrial_prefix).reset_index(
        drop=True
    )
    _validate_unique_rna_ids(result, str(path))
    return _set_processing_counts(result, raw_rows, identified_rna_ids)


def merge_tool_with_reference(
    reference: pd.DataFrame, tool: pd.DataFrame
) -> pd.DataFrame:
    """Left-join one unique tool result per exact reference RNA_ID."""
    _validate_unique_rna_ids(reference, "reference")
    _validate_unique_rna_ids(tool, "tool table")
    missing_ids = sorted(set(tool["RNA_ID"]) - set(reference["RNA_ID"]))
    if missing_ids:
        raise ValueError(f"Tool RNA_ID values absent from reference: {missing_ids[:5]}")
    result = reference.merge(
        tool.assign(_hit=True),
        on="RNA_ID",
        how="left",
        validate="one_to_one",
    )
    result["_hit"] = result["_hit"].eq(True)
    return result


def build_functional_table(
    reference: pd.DataFrame, tool_tables: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    """Build the RNA-level downstream table using the supplied tool variants."""
    result = reference.copy()
    for source in TOOL_ORDER:
        if source not in tool_tables:
            raise ValueError(f"Missing downstream tool table: {source}")
        tool = tool_tables[source]
        description_columns = [column for column in tool if column != "RNA_ID"]
        if description_columns != [source]:
            raise ValueError(f"{source} must have columns RNA_ID and {source}")
        merged = merge_tool_with_reference(reference, tool)
        result[source] = merged[source]
        result[f"{source}_returned"] = merged["_hit"]
    return result


def build_rna_presence(
    reference: pd.DataFrame, tool_tables: dict[str, pd.DataFrame]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build exact RNA_ID/source rows and a fixed-order Boolean membership table."""
    if tuple(tool_tables) != TOOL_ORDER:
        raise ValueError(f"Tool source order must be {TOOL_ORDER}")
    rows = [reference[["RNA_ID"]].assign(Source="Reference")]
    for source in TOOL_ORDER:
        merged = merge_tool_with_reference(reference, tool_tables[source])
        rows.append(merged.loc[merged["_hit"], ["RNA_ID"]].assign(Source=source))

    presence = pd.concat(rows, ignore_index=True).drop_duplicates(
        subset=["RNA_ID", "Source"]
    )
    presence["Source"] = pd.Categorical(
        presence["Source"], categories=SOURCE_ORDER, ordered=True
    )
    presence = presence.sort_values(["RNA_ID", "Source"], kind="stable").reset_index(
        drop=True
    )
    presence["Source"] = presence["Source"].astype("string")
    membership = (
        pd.crosstab(presence["RNA_ID"], presence["Source"])
        .reindex(index=reference["RNA_ID"], columns=SOURCE_ORDER, fill_value=0)
        .astype(bool)
    )
    membership.index.name = "RNA_ID"
    return presence, membership


def build_cleaning_summary(
    organism: str,
    display_name: str,
    reference: pd.DataFrame,
    tool_tables: dict[str, pd.DataFrame],
    submitted_rna_ids: set[str] | None = None,
) -> pd.DataFrame:
    """Summarize GFFread input, tool-identified, and retained RNA_ID counts."""
    if tuple(tool_tables) != TOOL_ORDER:
        raise ValueError(f"Tool source order must be {TOOL_ORDER}")
    reference_ids = set(reference["RNA_ID"])
    submitted_rna_ids = reference_ids if submitted_rna_ids is None else submitted_rna_ids
    if not submitted_rna_ids:
        raise ValueError("submitted_rna_ids must not be empty")
    rows = []
    for source in TOOL_ORDER:
        table = tool_tables[source]
        _validate_unique_rna_ids(table, source)
        identified_ids = set(table.attrs.get("identified_rna_ids", table["RNA_ID"]))
        identified_submitted = identified_ids & submitted_rna_ids
        matched_mask = table["RNA_ID"].isin(reference_ids)
        matched = int(matched_mask.sum())
        matched_ids = set(table.loc[matched_mask, "RNA_ID"])
        represented_loci = int(
            reference.loc[reference["RNA_ID"].isin(matched_ids), "locus_tag"].nunique()
        )
        rows.append(
            {
                "organism": organism,
                "display_name": display_name,
                "source": source,
                "gffread_input_rna_ids": len(submitted_rna_ids),
                "raw_rows": int(table.attrs.get("raw_rows", len(table))),
                "tool_identified_rna_ids": len(identified_submitted),
                "identified_outside_submission": len(
                    identified_ids - submitted_rna_ids
                ),
                "cleaned_rows": len(table),
                "reference_rna_ids": len(reference),
                "retained_rna_ids": matched,
                "represented_loci": represented_loci,
                "percentage_reference": 100 * matched / len(reference),
            }
        )
    return pd.DataFrame(rows)


def save_table(table: pd.DataFrame, path: Path, *, index: bool = False) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=index)


def plot_cleaning_summary(summary: pd.DataFrame, output_path: Path) -> None:
    """Save a four-tool by two-species stacked RNA_ID processing figure."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    organisms = summary[["organism", "display_name"]].drop_duplicates()
    if len(organisms) != 2:
        raise ValueError("Cleaning summary figure requires exactly two organisms")
    expected = pd.MultiIndex.from_product(
        [organisms["organism"], TOOL_ORDER], names=["organism", "source"]
    )
    indexed = summary.set_index(["organism", "source"])
    if not indexed.index.is_unique or not expected.equals(indexed.index):
        indexed = indexed.reindex(expected)
    if indexed.isna().any().any():
        raise ValueError("Cleaning summary is missing an organism/tool combination")

    processing_notes = {
        "Kofam": "Significant (*) hit; first hit per RNA_ID",
        "Pannzer": "DE records; first description per RNA_ID",
        "EggNOG": "First annotation per RNA_ID",
        "InterProScan": "Selected family DBs; minimum numeric e-value",
    }
    colors = ["#264653", "#457b9d", "#2a9d8f"]
    fig, axes = plt.subplots(4, 2, figsize=(11, 14), layout="constrained")
    organism_rows = organisms.reset_index(drop=True)
    for row_index, source in enumerate(TOOL_ORDER):
        for column_index, organism_row in organism_rows.iterrows():
            organism = organism_row["organism"]
            row = indexed.loc[(organism, source)]
            ax = axes[row_index, column_index]
            input_count = int(row["gffread_input_rna_ids"])
            identified_count = int(row["tool_identified_rna_ids"])
            retained_count = int(row["retained_rna_ids"])
            not_identified = input_count - identified_count
            identified_not_retained = identified_count - retained_count
            if min(not_identified, identified_not_retained, retained_count) < 0:
                raise ValueError(
                    f"Invalid processing counts for {organism}/{source}: "
                    f"input={input_count}, identified={identified_count}, "
                    f"retained={retained_count}"
                )
            segments = [not_identified, identified_not_retained, retained_count]
            labels = [
                "Not\nidentified",
                "Identified but\nnot retained",
                "Retained\nRNA_ID",
            ]
            bottom = 0
            for value, color in zip(segments, colors, strict=True):
                bar = ax.bar([0], [value], bottom=bottom, color=color, width=0.15)
                if value:
                    ax.bar_label(
                        bar,
                        labels=[f"{value:,}"],
                        label_type="center",
                        fontsize=8,
                    )
                bottom += value
            ax.set_xticks([0], ["GFFread input\nRNA_ID"], fontsize=8)
            ax.set_ylim(0, input_count * 1.25)
            ax.text(
                0.01,
                0.97,
                f"Selection: {processing_notes[source]}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=7,
            )
            ax.spines[["top", "right"]].set_visible(False)
            if row_index == 0:
                italic_name = organism_row["display_name"].replace(" ", r"\ ")
                ax.set_title(rf"$\it{{{italic_name}}}$", fontsize=13, pad=12)
            if column_index == 0:
                ax.set_ylabel(f"{source}\nExact count", fontsize=10)
            else:
                ax.set_ylabel("Exact count", fontsize=9)
    legend = [
        Patch(color=colors[0], label="Not identified by tool"),
        Patch(color=colors[1], label="Identified but not retained"),
        Patch(color=colors[2], label="Retained RNA_ID"),
    ]
    fig.legend(
        handles=legend,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Tool-specific RNA_ID processing summary", fontsize=15, y=0.995)
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.91))
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_rna_upset(membership: pd.DataFrame, title: str, output_path: Path) -> None:
    """Save a black, fixed-source-order UpSet plot of exact RNA_ID elements."""
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_indicators

    if not membership.index.is_unique or membership.index.name != "RNA_ID":
        raise ValueError("UpSet membership must have unique exact RNA_ID elements")
    if tuple(membership.columns) != SOURCE_ORDER:
        raise ValueError(f"UpSet source order must be {SOURCE_ORDER}")
    upset_data = from_indicators(list(SOURCE_ORDER), membership)
    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts="%d",
        sort_by="cardinality",
        sort_categories_by="cardinality",
        element_size=34,
        facecolor="black",
    )
    figure = plt.figure(figsize=(11, 7), facecolor="white")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")
        upset.plot(fig=figure)
    figure.suptitle(title, fontsize=14, y=1.01)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def prepare_kofam_output(
    input_path: Path,
    output_path: Path,
    mitochondrial_prefix: str | None = None,
) -> pd.DataFrame:
    """Select the first significant Kofam hit per RNA_ID and write a TSV."""
    result = prepare_kofam(
        input_path,
        significant_only=True,
        mitochondrial_prefix=mitochondrial_prefix,
    )
    save_table(result, output_path)
    return result
