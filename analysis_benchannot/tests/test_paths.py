import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from benchannot.paths import (
    ANALYSIS_ROOT,
    REPOSITORY_ROOT,
    project_paths,
    repository_relative,
    validate_required_inputs,
)


class ProjectPathsTests(unittest.TestCase):
    def test_resolves_origin_layout(self):
        paths = project_paths("origin")

        self.assertEqual(paths.dataset, "origin")
        self.assertEqual(paths.repository_root, REPOSITORY_ROOT)
        self.assertEqual(paths.analysis_root, ANALYSIS_ROOT)
        self.assertEqual(paths.data_root, REPOSITORY_ROOT / "data")
        self.assertEqual(paths.dataset_root, REPOSITORY_ROOT / "data" / "origin")
        self.assertEqual(
            paths.genome_eukaryote,
            REPOSITORY_ROOT / "data" / "genome_eukaryote",
        )
        self.assertEqual(
            paths.genome_prokaryote,
            REPOSITORY_ROOT / "data" / "genome_prokaryote",
        )
        self.assertEqual(
            paths.eukaryote_tools,
            REPOSITORY_ROOT / "data" / "origin" / "eukaryote_output_tools",
        )
        self.assertEqual(
            paths.prokaryote_tools,
            REPOSITORY_ROOT / "data" / "origin" / "prokaryote_output_tools",
        )
        self.assertEqual(
            paths.analysis_output,
            ANALYSIS_ROOT / "2_run" / "output" / "origin",
        )

    def test_resolves_reproduced_layout(self):
        paths = project_paths("reproduced")

        self.assertEqual(
            paths.dataset_root,
            REPOSITORY_ROOT / "data" / "reproduced",
        )
        self.assertEqual(
            paths.eukaryote_tools,
            REPOSITORY_ROOT / "data" / "reproduced" / "eukaryote_output_tools",
        )
        self.assertEqual(
            paths.prokaryote_tools,
            REPOSITORY_ROOT / "data" / "reproduced" / "prokaryote_output_tools",
        )
        self.assertEqual(
            paths.analysis_output,
            ANALYSIS_ROOT / "2_run" / "output" / "reproduced",
        )

    def test_rejects_invalid_dataset(self):
        with self.assertRaisesRegex(
            ValueError,
            "Invalid dataset 'other'. Expected 'origin' or 'reproduced'.",
        ):
            project_paths("other")  # type: ignore[arg-type]

    def test_resolution_is_independent_of_current_working_directory(self):
        expected = project_paths("origin")
        previous_cwd = Path.cwd()

        with TemporaryDirectory() as directory:
            try:
                os.chdir(directory)
                actual = project_paths("origin")
            finally:
                os.chdir(previous_cwd)

        self.assertEqual(actual, expected)

    def test_does_not_create_directories(self):
        with patch.object(Path, "mkdir") as mkdir:
            project_paths("reproduced")

        mkdir.assert_not_called()

    def test_datasets_have_distinct_output_directories(self):
        origin = project_paths("origin")
        reproduced = project_paths("reproduced")

        self.assertNotEqual(origin.analysis_output, reproduced.analysis_output)


class RepositoryRelativeTests(unittest.TestCase):
    def test_formats_existing_repository_path(self):
        self.assertEqual(
            repository_relative(REPOSITORY_ROOT / "main.nf"),
            "main.nf",
        )

    def test_formats_nonexistent_repository_path(self):
        path = REPOSITORY_ROOT / "data" / "origin" / "not-created.tsv"

        self.assertEqual(
            repository_relative(path),
            "data/origin/not-created.tsv",
        )

    def test_resolves_relative_path_from_repository_root(self):
        previous_cwd = Path.cwd()

        with TemporaryDirectory() as directory:
            try:
                os.chdir(directory)
                relative = repository_relative("data/origin/file.tsv")
            finally:
                os.chdir(previous_cwd)

        self.assertEqual(relative, "data/origin/file.tsv")

    def test_rejects_path_outside_repository(self):
        with TemporaryDirectory() as directory:
            with self.assertRaisesRegex(
                ValueError,
                "Path is outside the repository:",
            ):
                repository_relative(directory)


class ValidateRequiredInputsTests(unittest.TestCase):
    def test_accepts_existing_regular_files(self):
        validate_required_inputs(
            {"Workflow entrypoint": REPOSITORY_ROOT / "main.nf"},
            context="test analysis",
        )

    def test_reports_missing_input(self):
        missing = REPOSITORY_ROOT / "data" / "__missing_test_input__.tsv"

        with self.assertRaises(FileNotFoundError) as raised:
            validate_required_inputs(
                {"Tool output": missing},
                context="test analysis",
            )

        message = str(raised.exception)
        self.assertIn("Invalid required inputs for test analysis:", message)
        self.assertIn("Missing:", message)
        self.assertIn(
            f"- Tool output: {repository_relative(missing)}",
            message,
        )

    def test_reports_non_file_input(self):
        invalid = REPOSITORY_ROOT / "data"

        with self.assertRaises(FileNotFoundError) as raised:
            validate_required_inputs(
                {"Tool output": invalid},
                context="test analysis",
            )

        message = str(raised.exception)
        self.assertIn("Not regular files:", message)
        self.assertIn(
            f"- Tool output: {repository_relative(invalid)}",
            message,
        )

    def test_aggregates_missing_and_non_file_inputs(self):
        missing = REPOSITORY_ROOT / "data" / "__missing_test_input__.tsv"
        invalid = REPOSITORY_ROOT / "data"

        with self.assertRaises(FileNotFoundError) as raised:
            validate_required_inputs(
                {
                    "Missing output": missing,
                    "Directory output": invalid,
                },
                context="combined preflight",
            )

        message = str(raised.exception)
        self.assertIn("Missing:", message)
        self.assertIn("Not regular files:", message)
        self.assertIn("- Missing output:", message)
        self.assertIn("- Directory output:", message)

    def test_resolves_relative_inputs_from_repository_root(self):
        previous_cwd = Path.cwd()
        with TemporaryDirectory() as directory:
            try:
                os.chdir(directory)
                validate_required_inputs(
                    {"Workflow entrypoint": "main.nf"},
                    context="relative input",
                )
            finally:
                os.chdir(previous_cwd)

    def test_does_not_create_missing_inputs(self):
        missing = (
            REPOSITORY_ROOT
            / "data"
            / "__missing_test_directory__"
            / "input.tsv"
        )

        with self.assertRaises(FileNotFoundError):
            validate_required_inputs(
                {"Tool output": missing},
                context="side effect check",
            )

        self.assertFalse(missing.parent.exists())
