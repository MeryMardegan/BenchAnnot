from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from benchannot.eukaryotic.prepare_output_tool import (
    SOURCE_ORDER,
    TOOL_ORDER,
    build_cleaning_summary,
    build_functional_table,
    build_rna_presence,
    prepare_interproscan,
    prepare_kofam,
    prepare_kofam_output,
    prepare_pannzer,
    prepare_reference,
    read_fasta_ids,
    read_kofam,
)


class PrepareKofamOutputTests(unittest.TestCase):
    def test_reads_unique_exact_fasta_ids(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "submitted.faa"
            path.write_text(">rna-1 description\nMA\n>rna-2\nMK\n", encoding="utf-8")

            self.assertEqual(read_fasta_ids(path), {"rna-1", "rna-2"})

    def test_pannzer_header_is_not_an_rna_id(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "pannzer.txt"
            path.write_text(
                "qpid\ttype\tscore\tPPV\tid\tdesc\n"
                "rna-1\tDE\t1\t1\t1\tdescription\n",
                encoding="utf-8",
            )

            result = prepare_pannzer(path)

            self.assertEqual(result["RNA_ID"].tolist(), ["rna-1"])
            self.assertEqual(result.attrs["identified_rna_ids"], {"rna-1"})

    def test_selects_first_significant_hit(self):
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "kofam.txt"
            input_path.write_text(
                "# header\n"
                "\trna-1\tK00001\t10\t9\t0.1\t\"non-significant\"\n"
                "*\trna-1\tK00002\t10\t20\t1e-6\t\"first significant\"\n"
                "*\trna-1\tK00003\t10\t19\t2e-6\t\"second significant\"\n"
                "*\trna-Q1\tK00004\t10\t20\t1e-6\t\"mitochondrial\"\n",
                encoding="utf-8",
            )
            output_path = Path(directory) / "prepared.tsv"

            result = prepare_kofam_output(input_path, output_path, "rna-Q")

            expected = pd.DataFrame(
                {"RNA_ID": ["rna-1"], "Kofam": ["first significant"]}
            )
            pd.testing.assert_frame_equal(result, expected)
            pd.testing.assert_frame_equal(
                pd.read_csv(output_path, sep="\t"), expected
            )

    def test_rejects_wrong_column_count(self):
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "invalid.txt"
            input_path.write_text("*\trna-1\tK00001\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "expected 7"):
                read_kofam(input_path)

    def test_reference_retains_protein_provenance_and_uses_rna_key(self):
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "reference.tsv"
            input_path.write_text(
                "protein_id\tRNA_ID\tlocus_tag\tproduct\n"
                "XP_1\trna-1\tgene-1\tkinase\n",
                encoding="utf-8",
            )

            result = prepare_reference(input_path)

            self.assertEqual(
                result.columns.tolist(),
                ["protein_id", "RNA_ID", "locus_tag", "Reference"],
            )
            self.assertEqual(result.loc[0, "protein_id"], "XP_1")

    def test_functional_kofam_excludes_non_significant_candidates(self):
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "kofam.txt"
            input_path.write_text(
                "\trna-1\tK00001\t10\t30\t1e-20\tnon-significant\n"
                "*\trna-1\tK00002\t10\t20\t1e-10\tsignificant\n",
                encoding="utf-8",
            )

            result = prepare_kofam(input_path, significant_only=True)

            self.assertEqual(result.loc[0, "Kofam"], "significant")

    def test_interproscan_uses_stable_numeric_minimum_evalue(self):
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "interpro.tsv"

            def row(rna_id, database, evalue, description):
                fields = [
                    rna_id,
                    "md5",
                    "100",
                    database,
                    "signature",
                    "signature description",
                    "1",
                    "100",
                    evalue,
                    "T",
                    "date",
                    "IPR000001",
                    description,
                ]
                return "\t".join(fields)

            input_path.write_text(
                "\n".join(
                    [
                        row("rna-1", "PANTHER", "2.7E-167", "numeric minimum"),
                        row("rna-1", "PIRSF", "1.5E-15", "lexical minimum"),
                        row("rna-2", "PANTHER", "1e-20", "first tied row"),
                        row("rna-2", "PIRSF", "1e-20", "second tied row"),
                        row("rna-3", "Pfam", "1e-100", "excluded database"),
                        row("rna-3", "Hamap", "1e-5", "included database"),
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            result = prepare_interproscan(input_path).set_index("RNA_ID")

            self.assertEqual(result.loc["rna-1", "InterProScan"], "numeric minimum")
            self.assertEqual(result.loc["rna-2", "InterProScan"], "first tied row")
            self.assertEqual(result.loc["rna-3", "InterProScan"], "included database")

    def test_presence_uses_exact_rna_ids_in_fixed_source_order(self):
        reference = pd.DataFrame(
            {
                "protein_id": ["p1", "p2", "p3"],
                "RNA_ID": ["rna-1", "rna-2", "rna-3"],
                "locus_tag": ["gene-1", "gene-1", "gene-2"],
                "Reference": ["product 1", "product 1", "product 2"],
            }
        )
        tools = {
            "Kofam": pd.DataFrame({"RNA_ID": ["rna-2"], "Kofam": [""]}),
            "Pannzer": pd.DataFrame({"RNA_ID": ["rna-1"], "Pannzer": ["-"]}),
            "EggNOG": pd.DataFrame({"RNA_ID": ["rna-3"], "EggNOG": [None]}),
            "InterProScan": pd.DataFrame(
                {"RNA_ID": ["rna-2"], "InterProScan": ["family"]}
            ),
        }

        presence, membership = build_rna_presence(reference, tools)

        self.assertEqual(tuple(membership.columns), SOURCE_ORDER)
        self.assertEqual(membership.index.tolist(), reference["RNA_ID"].tolist())
        self.assertTrue(membership.loc["rna-2", "Kofam"])
        self.assertFalse(membership.loc["rna-1", "Kofam"])
        self.assertEqual(int(membership["Kofam"].sum()), 1)
        self.assertIn(["rna-3", "EggNOG"], presence.astype(str).values.tolist())

    def test_first_reference_transcript_presence_can_be_built_from_selected_ids(self):
        reference = pd.DataFrame(
            {
                "protein_id": ["p1", "p2", "p3"],
                "RNA_ID": ["rna-1", "rna-2", "rna-3"],
                "locus_tag": ["gene-1", "gene-1", "gene-2"],
                "Reference": ["first", "second", "third"],
            }
        )
        tools = {
            source: pd.DataFrame(
                {"RNA_ID": ["rna-1", "rna-2", "rna-3"], source: ["hit"] * 3}
            )
            for source in TOOL_ORDER
        }
        first_reference = reference.drop_duplicates("locus_tag", keep="first")
        first_ids = set(first_reference["RNA_ID"])
        first_tools = {
            source: table[table["RNA_ID"].isin(first_ids)].copy()
            for source, table in tools.items()
        }

        _, membership = build_rna_presence(first_reference, first_tools)

        self.assertEqual(membership.index.tolist(), ["rna-1", "rna-3"])
        self.assertEqual(int(membership["Reference"].sum()), 2)
        self.assertEqual(int(membership["Kofam"].sum()), 2)

    def test_summary_counts_rows_even_when_descriptions_are_empty(self):
        reference = pd.DataFrame(
            {
                "protein_id": ["p1", "p2"],
                "RNA_ID": ["rna-1", "rna-2"],
                "locus_tag": ["gene-1", "gene-1"],
                "Reference": ["product", "product"],
            }
        )
        tools = {
            source: pd.DataFrame({"RNA_ID": ["rna-2"], source: [value]})
            for source, value in zip(
                ("Kofam", "Pannzer", "EggNOG", "InterProScan"),
                ("", "-", None, "family"),
                strict=True,
            )
        }
        tools["Kofam"].attrs["raw_rows"] = 3

        summary = build_cleaning_summary("test", "Test species", reference, tools)
        functional = build_functional_table(reference, tools)

        self.assertEqual(summary["gffread_input_rna_ids"].tolist(), [2, 2, 2, 2])
        self.assertEqual(summary["tool_identified_rna_ids"].tolist(), [1, 1, 1, 1])
        self.assertEqual(summary["identified_outside_submission"].tolist(), [0, 0, 0, 0])
        self.assertEqual(summary["retained_rna_ids"].tolist(), [1, 1, 1, 1])
        self.assertEqual(summary.loc[0, "raw_rows"], 3)
        self.assertEqual(functional.loc[1, "RNA_ID"], "rna-2")
        for source in tools:
            self.assertIn(source, functional)
            self.assertIn(f"{source}_returned", functional)
            self.assertFalse(functional.loc[0, f"{source}_returned"])
            self.assertTrue(functional.loc[1, f"{source}_returned"])
