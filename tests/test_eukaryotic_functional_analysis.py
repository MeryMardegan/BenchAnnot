import unittest

import pandas as pd

from benchannot.eukaryotic.functional_analysis import (
    add_reference_baseline,
    build_transcript_annotation_table,
    categorize_functional_shifts,
    is_informative,
    select_first_reference_transcript,
    separate_isoform_suffix,
)
from benchannot.prokaryotic.functional_analysis import is_informative as prok_is_informative


class EukaryoticFunctionalAnalysisTests(unittest.TestCase):
    def test_informative_rule_matches_prokaryotic_rule(self):
        descriptions = [
            None,
            "",
            "-",
            "hypothetical protein",
            "uncharacterized protein",
            "putative kinase",
            "predicted kinase",
            "domain of unknown function",
            "DNA helicase",
        ]
        self.assertEqual(
            [is_informative(value) for value in descriptions],
            [prok_is_informative(value) for value in descriptions],
        )
        self.assertFalse(is_informative("-"))
        self.assertFalse(is_informative("predicted kinase"))

    def test_separates_only_terminal_isoform_suffix(self):
        self.assertEqual(
            separate_isoform_suffix("protein kinase, isoform X2"),
            ("protein kinase", "X2"),
        )
        self.assertEqual(
            separate_isoform_suffix("protein kinase, transcript variant B"),
            ("protein kinase", "B"),
        )
        self.assertEqual(
            separate_isoform_suffix("acetyl-CoA synthetase isoform"),
            ("acetyl-CoA synthetase isoform", pd.NA),
        )

    def test_builds_transcript_table_and_selects_first_reference_transcript(self):
        reference = pd.DataFrame(
            {
                "RNA_ID": ["rna-1", "rna-2", "rna-3"],
                "protein_id": ["NP_1.1", "NP_2.1", "NP_3.1"],
                "locus_tag": ["gene-1", "gene-1", "gene-2"],
                "Reference": [
                    "protein alpha, isoform B",
                    "protein alpha, isoform A",
                    "hypothetical protein",
                ],
            }
        )
        kofam = pd.DataFrame(
            {"RNA_ID": ["rna-1", "rna-3"], "Kofam": ["alpha", "kinase"]}
        )

        all_transcripts = build_transcript_annotation_table(reference, {"Kofam": kofam})
        representative = select_first_reference_transcript(all_transcripts)

        self.assertEqual(all_transcripts["Reference"].tolist()[0], "protein alpha")
        self.assertEqual(all_transcripts["Reference_isoform"].tolist()[0], "B")
        self.assertEqual(representative["RNA_ID"].tolist(), ["rna-1", "rna-3"])

    def test_categorizes_functional_shifts(self):
        table = pd.DataFrame(
            {
                "RNA_ID": ["rna-1", "rna-2", "rna-3", "rna-4"],
                "protein_id": ["NP_1.1", "NP_2.1", "NP_3.1", "NP_4.1"],
                "locus_tag": ["g1", "g2", "g3", "g4"],
                "Reference": [
                    "enzyme",
                    "hypothetical protein",
                    "enzyme",
                    "uncharacterized protein",
                ],
                "Kofam": [
                    "enzyme family",
                    "kinase",
                    None,
                    "hypothetical protein",
                ],
            }
        )

        summary, details = categorize_functional_shifts(table, ["Kofam"])
        with_reference = add_reference_baseline(summary, table)

        self.assertEqual(summary.loc["Kofam", "Both Informative"], 1)
        self.assertEqual(summary.loc["Kofam", "Both Hypothetical"], 1)
        self.assertEqual(summary.loc["Kofam", "Over-Annotation"], 1)
        self.assertEqual(summary.loc["Kofam", "Under-Annotation"], 1)
        self.assertEqual(len(details), 4)
        self.assertTrue((with_reference.sum(axis=1) == 4).all())


if __name__ == "__main__":
    unittest.main()
