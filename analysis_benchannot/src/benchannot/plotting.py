"""Shared scientific plotting style and multi-organism figure renderers.

The identifiers and visual grammar intentionally match the prokaryotic
identifier-matching analysis: seaborn ``ticks``/``paper`` typography, thin
frames, white bar edges, italic species titles, a single top legend, and
300-dpi white-background PNGs. Diverging functional-shift bars use matplotlib
primitives because seaborn has no native diverging-stack mark, but they share
the same rcParams, palette, and label conventions.
"""

from __future__ import annotations

from math import ceil
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


TEAL = "#2a9d8f"
NAVY = "#264653"
ORANGE = "#e76f51"
SAND = "#f4a261"
BLUE = "#457b9d"
RUST = "#c8553d"

CATEGORICAL_COLORS = (TEAL, ORANGE, NAVY, SAND, BLUE, RUST)

FUNCTIONAL_SHIFT_ORDER = (
    "Both Informative",
    "Both Hypothetical",
    "Over-Annotation",
    "Under-Annotation",
)

FUNCTIONAL_SHIFT_LABELS = {
    "Both Informative": "Agreement (Informative)",
    "Both Hypothetical": "Agreement (Non-informative)",
    "Over-Annotation": "Over-Annotation",
    "Under-Annotation": "Under-Annotation",
}

FUNCTIONAL_SHIFT_COLORS = {
    "Both Informative": TEAL,
    "Both Hypothetical": NAVY,
    "Over-Annotation": ORANGE,
    "Under-Annotation": SAND,
}

STAGE_COLORS = (NAVY, BLUE, TEAL, RUST)

TOOL_DISPLAY = {
    "Reference": "NCBI RefSeq reference",
    "NCBI RefSeq reference": "NCBI RefSeq reference",
    "Prokka": "Prokka",
    "Bakta": "Bakta",
    "PGAP": "PGAP",
    "eggNOG": "eggNOG",
    "EggNOG": "eggNOG",
    "Kofam": "Kofam",
    "Pannzer": "Pannzer",
    "InterProScan": "InterProScan",
}

PLOT_NAMES = {
    "bacillus_subtilis": r"$\it{B.\ subtilis}$",
    "escherichia_coli": r"$\it{E.\ coli}$",
    "saccharomyces_cerevisiae": r"$\it{S.\ cerevisiae}$",
    "drosophila_melanogaster": r"$\it{D.\ melanogaster}$",
}


def apply_style() -> None:
    """Apply the shared BenchAnnot plotting style."""
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.6,
            "axes.edgecolor": "#555555",
            "axes.labelcolor": "#333333",
            "text.color": "#333333",
            "xtick.color": "#333333",
            "ytick.color": "#333333",
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "axes.grid": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )


def display_name(name: str) -> str:
    """Return the canonical display label for a tool or reference source."""
    return TOOL_DISPLAY.get(str(name), str(name))


def italic_name(name: str) -> str:
    """Render a species name in mathtext italics, preserving already-styled text."""
    text = str(name)
    if text.startswith("$"):
        return text
    return r"$\it{" + text.replace(" ", r"\ ") + "}$"


def panel_title(
    axis: plt.Axes,
    name: str,
    *,
    panel_index: int | None = None,
    fontsize: float = 8.5,
    pad: float = 8,
) -> None:
    """Set an italic panel title with optional panel letter (A, B, etc.)."""
    title_text = italic_name(PLOT_NAMES.get(name, name))
    if panel_index is not None and panel_index < 26:
        letter = chr(ord("A") + panel_index)
        title_text = f"{letter}) {title_text}"
    axis.set_title(title_text, fontsize=fontsize, pad=pad)


def denominator_note(axis: plt.Axes, text: str) -> None:
    """Place a small denominator note in the upper-left of a panel."""
    axis.text(
        0.0,
        1.01,
        text,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=7,
        color="#666666",
    )


def style_axis(axis: plt.Axes, *, hide_left: bool = False) -> None:
    """Apply the shared axis frame style."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    if hide_left:
        axis.spines["left"].set_visible(False)
    for spine in axis.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("#555555")
    axis.tick_params(length=2.5, width=0.6, colors="#333333")
    axis.set_axisbelow(True)
    axis.grid(False)


def save_figure(
    figure: plt.Figure,
    output_path: Path | str,
    *,
    show: bool = True,
) -> Path:
    """Save a figure at publication resolution and optionally display it."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved {output_path}")
    if show:
        plt.show()
    plt.close(figure)
    return output_path


def _panel_grid(
    panel_count: int,
    ncols: int | None,
    *,
    width: float,
    height: float,
    sharey: bool = True,
    wspace: float | None = None,
    hspace: float | None = None,
):
    ncols = ncols or panel_count
    nrows = ceil(panel_count / ncols)
    gridspec_kw = {}
    if wspace is not None:
        gridspec_kw["wspace"] = wspace
    if hspace is not None:
        gridspec_kw["hspace"] = hspace
    figure, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(width * ncols, height * nrows),
        sharey=sharey,
        squeeze=False,
        gridspec_kw=gridspec_kw if gridspec_kw else None,
    )
    return figure, axes.flatten()


def _top_legend(figure: plt.Figure, axis: plt.Axes, *, ncol: int = 4) -> None:
    handles, labels = axis.get_legend_handles_labels()
    legend = axis.get_legend()
    if legend is not None:
        legend.remove()
    if not handles:
        return
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=min(ncol, len(labels)),
        frameon=False,
        fontsize=8,
    )


def plot_rate_panels(
    panels: Mapping[str, pd.DataFrame],
    *,
    ylabel: str,
    title: str | None = None,
    source_order: Sequence[str] | None = None,
    metric_order: Sequence[str] | None = None,
    palette: Mapping[str, str] | Sequence[str] | None = None,
    ylim: tuple[float, float] | None = None,
    value_format: str = "{:.1f}",
    panel_notes: Mapping[str, str] | None = None,
    ncols: int | None = None,
    panel_width: float = 3.4,
    panel_height: float = 3.2,
    bar_width: float = 0.65,
    sharey: bool = True,
    wspace: float | None = None,
    label_fontsize: float = 6,
    label_rotation: float = 0,
    closed_frame: bool = False,
    legend_labels: Sequence[str] | None = None,
    legend_bbox: tuple[float, float] = (0.5, 1.04),
    save_path: Path | str | None = None,
    show: bool = True,
) -> Path | None:
    """Render side-by-side grouped-rate panels from tidy panel tables.

    Each panel table must contain ``Source``, ``Metric``, and ``Value`` columns.
    Panel keys become italic species titles.
    """
    if not panels:
        raise ValueError("At least one panel is required")
    figure, axes = _panel_grid(
        len(panels),
        ncols,
        width=panel_width,
        height=panel_height,
        sharey=sharey,
        wspace=wspace,
    )
    for index, (panel_key, table) in enumerate(panels.items()):
        axis = axes[index]
        sns.barplot(
            data=table,
            x="Source",
            y="Value",
            hue="Metric",
            order=list(source_order) if source_order is not None else None,
            hue_order=list(metric_order) if metric_order is not None else None,
            palette=palette,
            width=bar_width,
            saturation=1.0,
            gap=0,
            edgecolor="white",
            linewidth=0.4,
            ax=axis,
        )
        for container in axis.containers:
            labels = [
                value_format.format(value) if value else ""
                for value in container.datavalues
            ]
            axis.bar_label(
                container,
                labels=labels,
                padding=2,
                fontsize=label_fontsize,
                rotation=label_rotation,
                color="#333333",
            )
        if ylim is not None:
            axis.set_ylim(*ylim)
        elif label_rotation != 0:
            current_ymax = table["Value"].max()
            if pd.notna(current_ymax) and current_ymax > 0:
                axis.set_ylim(0, current_ymax * 1.25)
        axis.set_xlabel("")
        axis.set_ylabel(ylabel if index % (ncols or len(panels)) == 0 else "")
        panel_title(axis, panel_key, panel_index=index)
        if panel_notes and panel_key in panel_notes:
            denominator_note(axis, panel_notes[panel_key])
        if closed_frame:
            for spine in axis.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.6)
                spine.set_color("#555555")
            axis.tick_params(length=2, width=0.5, colors="#333333")
            axis.set_axisbelow(True)
            axis.grid(False)
        else:
            style_axis(axis)
        if axis.get_legend() is not None:
            axis.get_legend().remove()
    handles, labels = axes[0].get_legend_handles_labels()
    display_labels = list(legend_labels) if legend_labels is not None else labels
    if handles:
        figure.legend(
            handles,
            display_labels,
            loc="upper center",
            bbox_to_anchor=legend_bbox,
            ncol=min(len(display_labels), 4),
            frameon=False,
            fontsize=8,
        )
    if title:
        figure.suptitle(title, fontsize=11, y=1.08)
    figure.subplots_adjust(top=0.88, bottom=0.15)
    if save_path is None:
        if show:
            plt.show()
        plt.close(figure)
        return None
    return save_figure(figure, save_path, show=show)


def plot_stacked_panels(
    panels: Mapping[str, pd.DataFrame],
    *,
    categories: Sequence[str],
    colors: Sequence[str],
    labels: Sequence[str] | None = None,
    ylabel: str,
    title: str | None = None,
    ncols: int = 1,
    bar_width: float = 0.62,
    save_path: Path | str | None = None,
    show: bool = True,
) -> Path | None:
    """Render stacked processing-category panels (organisms as rows by default).

    Each panel table is indexed by source and contains one column per stack
    category, in stacking order.
    """
    if not panels:
        raise ValueError("At least one panel is required")
    labels = labels or categories
    figure, axes = _panel_grid(len(panels), ncols, width=5.0, height=4.2)
    for index, (panel_key, table) in enumerate(panels.items()):
        axis = axes[index]
        missing = set(categories) - set(table.columns)
        if missing:
            raise ValueError(f"{panel_key} is missing categories: {sorted(missing)}")
        positions = np.arange(len(table))
        totals = table[list(categories)].sum(axis=1).to_numpy()
        bottom = np.zeros(len(table))
        for category, color, label in zip(categories, colors, labels, strict=True):
            values = table[category].to_numpy()
            bars = axis.bar(
                positions,
                values,
                bottom=bottom,
                color=color,
                label=label,
                width=bar_width,
                edgecolor="white",
                linewidth=0.4,
            )
            for bar, value, total in zip(bars, values, totals, strict=True):
                if value <= 0 or total <= 0:
                    continue
                fraction = value / total
                if fraction < 0.04:
                    continue
                axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    f"{int(value):,}",
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="white" if fraction >= 0.15 else "#222222",
                )
            bottom += values
        axis.set_xticks(
            positions,
            [display_name(source) for source in table.index],
            fontsize=8,
        )
        axis.set_ylim(0, max(totals) * 1.15 if len(totals) else 1)
        axis.set_ylabel(ylabel if index % ncols == 0 else "")
        panel_title(axis, panel_key, panel_index=index)
        style_axis(axis)
        if index % ncols != 0:
            axis.set_ylabel("")
    _top_legend(figure, axes[0], ncol=len(labels))
    if title:
        figure.suptitle(title, fontsize=11, y=1.06)
    figure.tight_layout()
    if save_path is None:
        if show:
            plt.show()
        plt.close(figure)
        return None
    return save_figure(figure, save_path, show=show)


def plot_functional_shift_panels(
    panels: Mapping[str, pd.DataFrame],
    *,
    count_label: str,
    save_path: Path | str | None = None,
    show: bool = True,
) -> Path | None:
    """Render side-by-side diverging functional-shift panels.

    Each panel table is indexed by source (reference first) and contains the
    four functional-shift categories. Every source must share the same total.
    """
    if not panels:
        raise ValueError("At least one panel is required")
    figure, axes = _panel_grid(len(panels), None, width=5.4, height=3.8)
    for index, (panel_key, summary) in enumerate(panels.items()):
        axis = axes[index]
        missing = set(FUNCTIONAL_SHIFT_ORDER) - set(summary.columns)
        if missing:
            raise ValueError(f"{panel_key} is missing categories: {sorted(missing)}")
        totals = summary[list(FUNCTIONAL_SHIFT_ORDER)].sum(axis=1)
        if totals.nunique() != 1:
            raise ValueError("Every source must be evaluated against the same universe")
        total = int(totals.iloc[0])
        plot_table = summary[list(FUNCTIONAL_SHIFT_ORDER)].copy()
        plot_table["Over-Annotation"] *= -1
        plot_table["Under-Annotation"] *= -1

        over = axis.barh(
            plot_table.index,
            plot_table["Over-Annotation"],
            color=FUNCTIONAL_SHIFT_COLORS["Over-Annotation"],
            label=FUNCTIONAL_SHIFT_LABELS["Over-Annotation"],
        )
        under = axis.barh(
            plot_table.index,
            plot_table["Under-Annotation"],
            left=plot_table["Over-Annotation"],
            color=FUNCTIONAL_SHIFT_COLORS["Under-Annotation"],
            label=FUNCTIONAL_SHIFT_LABELS["Under-Annotation"],
        )
        informative = axis.barh(
            plot_table.index,
            plot_table["Both Informative"],
            color=FUNCTIONAL_SHIFT_COLORS["Both Informative"],
            label=FUNCTIONAL_SHIFT_LABELS["Both Informative"],
        )
        hypothetical = axis.barh(
            plot_table.index,
            plot_table["Both Hypothetical"],
            left=plot_table["Both Informative"],
            color=FUNCTIONAL_SHIFT_COLORS["Both Hypothetical"],
            label=FUNCTIONAL_SHIFT_LABELS["Both Hypothetical"],
        )
        axis.axvline(0, color="#555555", linewidth=0.8, linestyle="--")
        # Annotate disagreement side (Over-Annotation and Under-Annotation)
        for row_idx, source in enumerate(plot_table.index):
            over_val = abs(summary.loc[source, "Over-Annotation"])
            under_val = abs(summary.loc[source, "Under-Annotation"])
            over_pct = over_val / total * 100
            under_pct = under_val / total * 100
            y = row_idx

            # If both disagreement segments are narrow (< 6%), combine them at the stack tip
            if over_pct < 6 and under_pct < 6:
                tip = -(over_val + under_val)
                parts = []
                if over_pct > 0:
                    parts.append(f"{over_pct:.1f}%")
                if under_pct > 0:
                    parts.append(f"{under_pct:.1f}%")
                if parts:
                    axis.annotate(
                        " | ".join(parts),
                        xy=(tip, y),
                        xytext=(-4, 0),
                        textcoords="offset points",
                        ha="right",
                        va="center",
                        fontsize=7,
                        fontweight="bold",
                        clip_on=False,
                    )
            else:
                # Label Over-Annotation
                is_yeast_interpro = (
                    source == "InterProScan"
                    and panel_key in {"saccharomyces_cerevisiae", "Saccharomyces cerevisiae"}
                )
                if over_pct >= 6 or is_yeast_interpro:
                    axis.text(
                        -over_val / 2,
                        y,
                        f"{over_pct:.1f}%",
                        ha="center",
                        va="center",
                        fontsize=7,
                        fontweight="bold",
                        color="#333333",
                    )
                elif over_pct > 0:
                    axis.annotate(
                        f"{over_pct:.1f}%",
                        xy=(0, y),
                        xytext=(-4, 0),
                        textcoords="offset points",
                        ha="right",
                        va="center",
                        fontsize=7,
                        fontweight="bold",
                        clip_on=False,
                    )

                # Label Under-Annotation
                if under_pct >= 6:
                    axis.text(
                        -(over_val + under_val / 2),
                        y,
                        f"{under_pct:.1f}%",
                        ha="center",
                        va="center",
                        fontsize=7,
                        fontweight="bold",
                        color="#333333",
                    )
                elif under_pct > 0:
                    axis.annotate(
                        f"{under_pct:.1f}%",
                        xy=(-(over_val + under_val), y),
                        xytext=(-4, 0),
                        textcoords="offset points",
                        ha="right",
                        va="center",
                        fontsize=7,
                        fontweight="bold",
                        clip_on=False,
                    )

        # Annotate agreement side (Both Informative and Both Hypothetical)
        agreement_containers = (
            ("Both Informative", informative),
            ("Both Hypothetical", hypothetical),
        )
        for category, container in agreement_containers:
            for bar, value in zip(container.patches, container.datavalues, strict=True):
                if not value:
                    continue
                percentage = abs(value) / total * 100
                if percentage >= 6:
                    axis.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_y() + bar.get_height() / 2,
                        f"{percentage:.1f}%",
                        ha="center",
                        va="center",
                        fontsize=7,
                        fontweight="bold",
                        color=(
                            "white"
                            if category == "Both Hypothetical"
                            else "#333333"
                        ),
                    )
                    continue
                axis.annotate(
                    f"{percentage:.1f}%",
                    xy=(bar.get_x() + bar.get_width(), bar.get_y() + bar.get_height() / 2),
                    xytext=(4, 0),
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    fontsize=7,
                    fontweight="bold",
                    clip_on=False,
                )
        axis.set_xlabel(
            f"{count_label} count\n"
            "\u2190 Disagreement        Agreement \u2192"
        )
        axis.set_ylabel("Annotation tool" if index == 0 else "")
        axis.xaxis.set_major_formatter(
            plt.FuncFormatter(lambda value, _: f"{int(abs(value)):,}")
        )
        axis.set_yticks(range(len(plot_table.index)))
        axis.set_yticklabels(
            [
                "Reference" if str(source) == "Reference" else display_name(source)
                for source in plot_table.index
            ]
        )
        panel_title(axis, panel_key, panel_index=index)
        denominator_note(axis, f"n = {total:,}")
        style_axis(axis, hide_left=True)
        axis.tick_params(axis="y", length=0, labelleft=index == 0, pad=12)
    handles, labels = axes[0].get_legend_handles_labels()
    legend_order = [1, 0, 2, 3]
    figure.legend(
        [handles[index] for index in legend_order],
        [labels[index] for index in legend_order],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=4,
        frameon=False,
        fontsize=8,
    )
    figure.tight_layout()
    if save_path is None:
        if show:
            plt.show()
        plt.close(figure)
        return None
    return save_figure(figure, save_path, show=show)



UPSET_LABELS = {"Reference": "Reference"}


def plot_upset_panels(
    memberships: Mapping[str, pd.DataFrame],
    *,
    source_order: Sequence[str],
    ncols: int = 2,
    element_size: float = 24,
    count_label: str = "CDS",
    panel_dpi: int = 300,
    save_path: Path | str | None = None,
    show: bool = True,
) -> Path | None:
    """Merge standard black UpSet plots side by side, one italic species label each.

    Each membership table is boolean (one row per element, one column per
    source). Every panel is a regular UpSet figure (same look as the original
    plots, with a smaller ``element_size``), rendered then tiled compactly.
    ``count_label`` names the counted element on the axes (e.g. "CDS", "RNA_ID").
    No figure title; each species label is centered over its panel.
    """
    from io import BytesIO
    import warnings

    from upsetplot import UpSet, from_indicators

    if not memberships:
        raise ValueError("At least one panel is required")
    labels = [UPSET_LABELS.get(source, display_name(source)) for source in source_order]
    images = []
    for idx, (key, membership) in enumerate(memberships.items()):
        table = membership[list(source_order)].astype(bool)
        table.columns = labels
        upset = UpSet(
            from_indicators(labels, table),
            subset_size="count",
            show_counts="%d",
            sort_by="cardinality",
            facecolor="black",
            element_size=element_size,
        )
        panel = plt.figure(facecolor="white")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")
            upset_axes = upset.plot(fig=panel)
        upset_axes["intersections"].set_ylabel(f"Intersection size ({count_label})")
        if upset_axes.get("totals") is not None:
            upset_axes["totals"].set_xlabel(f"Total {count_label}")
        letter = f"{chr(ord('A') + idx)}) " if len(memberships) > 1 and idx < 26 else ""
        panel.suptitle(f"{letter}{italic_name(PLOT_NAMES.get(key, key))}", fontsize=9, y=1.0)
        buffer = BytesIO()
        panel.savefig(buffer, format="png", dpi=panel_dpi, bbox_inches="tight", facecolor="white")
        plt.close(panel)
        buffer.seek(0)
        images.append(plt.imread(buffer, format="png"))

    ncols = max(1, min(ncols, len(images)))
    nrows = ceil(len(images) / ncols)
    width = sum(img.shape[1] for img in images[:ncols]) / panel_dpi
    height = max(img.shape[0] for img in images) / panel_dpi * nrows
    figure, axes = plt.subplots(
        nrows, ncols, figsize=(width, height), squeeze=False, facecolor="white",
        gridspec_kw={"wspace": 0.02, "hspace": 0.02},
    )
    for axis in axes.flat:
        axis.axis("off")
    for axis, image in zip(axes.flat, images):
        axis.imshow(image, interpolation="none")
    figure.subplots_adjust(left=0, right=1, top=1, bottom=0)
    if save_path is None:
        if show:
            plt.show()
        plt.close(figure)
        return None
    return save_figure(figure, save_path, show=show)
