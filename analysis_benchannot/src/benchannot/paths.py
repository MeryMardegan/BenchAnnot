from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Mapping


Dataset = Literal["origin", "reproduced"]

ANALYSIS_ROOT = Path(__file__).resolve().parents[2]
REPOSITORY_ROOT = ANALYSIS_ROOT.parent


@dataclass(frozen=True)
class ProjectPaths:
    dataset: Dataset
    repository_root: Path
    analysis_root: Path
    data_root: Path
    dataset_root: Path
    genome_eukaryote: Path
    genome_prokaryote: Path
    eukaryote_tools: Path
    prokaryote_tools: Path
    analysis_output: Path


def project_paths(dataset: Dataset) -> ProjectPaths:
    if dataset not in {"origin", "reproduced"}:
        raise ValueError(
            f"Invalid dataset '{dataset}'. "
            "Expected 'origin' or 'reproduced'."
        )

    data_root = REPOSITORY_ROOT / "data"
    dataset_root = data_root / dataset

    return ProjectPaths(
        dataset=dataset,
        repository_root=REPOSITORY_ROOT,
        analysis_root=ANALYSIS_ROOT,
        data_root=data_root,
        dataset_root=dataset_root,
        genome_eukaryote=data_root / "genome_eukaryote",
        genome_prokaryote=data_root / "genome_prokaryote",
        eukaryote_tools=dataset_root / "eukaryote_output_tools",
        prokaryote_tools=dataset_root / "prokaryote_output_tools",
        analysis_output=ANALYSIS_ROOT / "2_run" / "output" / dataset,
    )


def repository_relative(path: Path | str) -> str:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = REPOSITORY_ROOT / candidate
    resolved = candidate.resolve(strict=False)

    try:
        return resolved.relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError as exc:
        raise ValueError(
            f"Path is outside the repository: {resolved}"
        ) from exc


def validate_required_inputs(
    inputs: Mapping[str, Path | str],
    *,
    context: str,
) -> None:
    missing: list[tuple[str, Path]] = []
    invalid: list[tuple[str, Path]] = []

    for label, path in inputs.items():
        candidate = Path(path)
        if not candidate.is_absolute():
            candidate = REPOSITORY_ROOT / candidate

        if not candidate.exists():
            missing.append((label, candidate))
        elif not candidate.is_file():
            invalid.append((label, candidate))

    if not missing and not invalid:
        return

    lines = [f"Invalid required inputs for {context}:"]

    if missing:
        lines.append("Missing:")
        lines.extend(
            f"- {label}: {repository_relative(path)}"
            for label, path in missing
        )

    if invalid:
        lines.append("Not regular files:")
        lines.extend(
            f"- {label}: {repository_relative(path)}"
            for label, path in invalid
        )

    raise FileNotFoundError("\n".join(lines))
