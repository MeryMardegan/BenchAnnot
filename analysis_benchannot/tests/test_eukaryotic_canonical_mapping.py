import json
from hashlib import sha256
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from benchannot.eukaryotic.canonical_mapping import (
    build_canonical_functional_table,
    build_coverage_table,
    load_uniprot_cache,
    resolve_canonical_loci,
)


def snapshot(records):
    response = {"results": records}
    digest = sha256(
        json.dumps(response, ensure_ascii=True, separators=(",", ":"), sort_keys=True).encode(
            "ascii"
        )
    ).hexdigest()
    return {
        "metadata": {
            "proteome_id": "UP000000803",
            "request_url": "https://example.test/request",
            "response_url": "https://example.test/response",
            "retrieved_at": "2026-08-04T00:00:00Z",
            "response_sha256": digest,
            "record_count": len(records),
        },
        "response": response,
    }


def record(accession, refseq, sequence, *, reviewed=True, displayed=None):
    value = {
        "primaryAccession": accession,
        "entryType": (
            "UniProtKB reviewed (Swiss-Prot)"
            if reviewed
            else "UniProtKB unreviewed (TrEMBL)"
        ),
        "sequence": {"value": sequence},
        "uniProtKBCrossReferences": [{"database": "RefSeq", "id": refseq}],
    }
    if displayed is not None:
        value["isoformSequenceStatus"] = displayed
    return value


class CanonicalMappingTests(unittest.TestCase):
    def test_requires_reviewed_exact_version_and_sequence(self):
        reference = pd.DataFrame(
            {
                "RNA_ID": ["rna-1", "rna-2", "rna-3"],
                "protein_id": ["NP_1.1", "NP_2.1", "NP_3.1"],
                "locus_tag": ["g1", "g2", "g3"],
            }
        )
        payload = snapshot(
            [
                record("P1", "NP_1.1", "AAAA"),
                record("P2", "NP_2.1", "XXXX"),
                record("P3", "NP_3.1", "CCCC", reviewed=False),
                record("P4", "NP_3", "CCCC"),
            ]
        )
        mapping, selection = resolve_canonical_loci(
            reference,
            {"NP_1.1": "AAAA", "NP_2.1": "BBBB", "NP_3.1": "CCCC"},
            payload,
        )

        statuses = mapping.set_index("locus_tag")["resolution_status"].to_dict()
        self.assertEqual(statuses["g1"], "resolved_canonical")
        self.assertEqual(statuses["g2"], "sequence_conflict")
        self.assertEqual(statuses["g3"], "no_reviewed_mapping")
        self.assertEqual(selection["RNA_ID"].tolist(), ["rna-1"])

    def test_displayed_isoform_and_equivalent_transcript_tie(self):
        reference = pd.DataFrame(
            {
                "RNA_ID": ["rna-first", "rna-second"],
                "protein_id": ["NP_1.1", "NP_2.1"],
                "locus_tag": ["g1", "g1"],
            }
        )
        payload = snapshot(
            [
                record("P1-1", "NP_1.1", "AAAA", displayed="Displayed"),
                record("P1-2", "NP_2.1", "AAAA", displayed="Displayed"),
            ]
        )
        mapping, selection = resolve_canonical_loci(
            reference,
            {"NP_1.1": "AAAA", "NP_2.1": "AAAA"},
            payload,
            submitted_rna_ids={"rna-second"},
        )

        self.assertTrue((mapping["resolution_status"] == "equivalent_transcript_tie").all())
        self.assertEqual(selection.loc[0, "RNA_ID"], "rna-first")
        self.assertEqual(selection.loc[0, "submission_status"], "canonical_not_submitted")
        self.assertEqual(selection.loc[0, "equivalent_RNA_IDs"], "rna-first;rna-second")

    def test_multiple_exact_sequences_are_ambiguous(self):
        reference = pd.DataFrame(
            {
                "RNA_ID": ["rna-1", "rna-2"],
                "protein_id": ["NP_1.1", "NP_2.1"],
                "locus_tag": ["g1", "g1"],
            }
        )
        payload = snapshot(
            [record("P1", "NP_1.1", "AAAA"), record("P2", "NP_2.1", "BBBB")]
        )
        mapping, selection = resolve_canonical_loci(
            reference, {"NP_1.1": "AAAA", "NP_2.1": "BBBB"}, payload
        )
        self.assertTrue(
            (mapping["resolution_status"] == "ambiguous_reviewed_isoforms").all()
        )
        self.assertTrue(selection.empty)

    def test_canonical_denominator_keeps_missing_tool_annotations(self):
        functional = pd.DataFrame(
            {
                "RNA_ID": ["rna-1", "rna-2"],
                "protein_id": ["NP_1.1", "NP_2.1"],
                "locus_tag": ["g1", "g2"],
                "Reference": ["enzyme", "kinase"],
                "Kofam": ["-", "kinase"],
                "Kofam_returned": [True, True],
            }
        )
        selection = pd.DataFrame(
            {
                "locus_tag": ["g1", "g2"],
                "RNA_ID": ["rna-1", "rna-2"],
                "protein_id": ["NP_1.1", "NP_2.1"],
                "submission_status": ["canonical_not_submitted", "submitted"],
            }
        )
        canonical = build_canonical_functional_table(functional, selection)
        coverage = build_coverage_table(functional, selection, ["Kofam"])

        self.assertEqual(len(canonical), 2)
        self.assertEqual(canonical.loc[0, "Kofam"], "-")
        self.assertEqual(coverage.loc[1, "Recovered canonical RNA_ID"], 2)

    def test_cache_hash_is_verified(self):
        payload = snapshot([])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cache.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            loaded = load_uniprot_cache(path, "UP000000803")
            self.assertEqual(loaded["metadata"]["record_count"], 0)
            payload["metadata"]["response_sha256"] = "0" * 64
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                load_uniprot_cache(path)


if __name__ == "__main__":
    unittest.main()
