from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from benchannot.eukaryotic.audit_prepare_reference import (
    build_locus_reference,
    build_reference_audit,
    read_reference_fasta,
    scan_original_gff,
)


class ReferenceAuditTests(unittest.TestCase):
    def test_locus_reference_collapses_isoforms_without_discarding_values(self):
        transcripts = pd.DataFrame(
            [
                {"protein_id": "NP_1", "RNA_ID": "rna-1", "gene_id": "gene-1", "locus_tag": "L1", "product": "isoform A"},
                {"protein_id": "NP_2", "RNA_ID": "rna-2", "gene_id": "gene-1", "locus_tag": "L1", "product": "isoform B"},
            ]
        )

        loci = build_locus_reference(transcripts)

        self.assertEqual(len(loci), 1)
        self.assertEqual(loci.loc[0, "transcript_count"], 2)
        self.assertEqual(loci.loc[0, "RNA_IDs"], "rna-1 | rna-2")
        self.assertEqual(loci.loc[0, "products"], "isoform A | isoform B")

    def test_reads_reference_descriptions_and_rejects_empty_sequences(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "reference.faa"
            path.write_text(
                ">NP_1 protein alpha [Example species]\nMAAA\n"
                ">NP_2 protein beta, isoform B [Example species]\nMBBB\n",
                encoding="ascii",
            )

            result = read_reference_fasta(path)

            self.assertEqual(result["protein_id"].tolist(), ["NP_1", "NP_2"])
            self.assertEqual(
                result["product"].tolist(),
                ["protein alpha", "protein beta, isoform B"],
            )

    def test_reference_is_independent_of_gffread_inventory(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            gff = root / "reference.gff"
            filtered = root / "filtered.gff"
            reference_fasta = root / "reference.faa"
            gffread_fasta = root / "gffread.faa"
            output = root / "output"
            rows = [
                "chr1\tRefSeq\tmRNA\t1\t100\t.\t+\t.\tID=rna-1;Parent=gene-1;locus_tag=gene-1;product=mRNA product 1",
                "chr1\tRefSeq\tCDS\t1\t100\t.\t+\t0\tID=cds-NP_1;Parent=rna-1;protein_id=NP_1;locus_tag=gene-1;product=CDS product 1",
                "chr1\tRefSeq\tmRNA\t101\t200\t.\t+\t.\tID=rna-2;Parent=gene-2;locus_tag=gene-2;product=mRNA product 2",
                "chr1\tRefSeq\tCDS\t101\t200\t.\t+\t0\tID=cds-NP_2;Parent=rna-2;protein_id=NP_2;locus_tag=gene-2;product=CDS product 2",
                "mito\tRefSeq\tmRNA\t1\t30\t.\t+\t.\tID=rna-Q1;Parent=gene-Q;locus_tag=gene-Q;product=mitochondrial",
                "mito\tRefSeq\tCDS\t1\t30\t.\t+\t0\tID=cds-YP_1;Parent=rna-Q1;protein_id=YP_1;locus_tag=gene-Q;product=mitochondrial",
            ]
            text = "##gff-version 3\n" + "\n".join(rows) + "\n"
            gff.write_text(text, encoding="utf-8")
            filtered.write_text(text, encoding="utf-8")
            reference_fasta.write_text(
                ">NP_1 reference protein 1 [Species]\nMAAA\n"
                ">NP_2 reference protein 2 [Species]\nMBBB\n"
                ">YP_1 mitochondrial protein [Species]\nMCCC\n",
                encoding="ascii",
            )
            gffread_fasta.write_text(">rna-1 submitted\nMAAA\n", encoding="ascii")

            result = build_reference_audit(
                "example",
                gff,
                reference_fasta,
                filtered,
                gffread_fasta,
                "rna-Q",
                output,
            )

            self.assertEqual(
                result["final_reference"]["RNA_ID"].tolist(), ["rna-1", "rna-2"]
            )
            self.assertEqual(
                result["final_reference"]["product"].tolist(),
                ["reference protein 1", "reference protein 2"],
            )
            self.assertEqual(
                result["reference_not_in_gffread"]["RNA_ID"].tolist(), ["rna-2"]
            )
            self.assertEqual(
                result["mitochondrial_removed"]["RNA_ID"].tolist(), ["rna-Q1"]
            )
            self.assertEqual(
                result["transcript_reference"]["gene_id"].tolist(),
                ["gene-1", "gene-2"],
            )
            self.assertEqual(len(result["locus_reference"]), 2)
            self.assertTrue(
                (output / "audit/gff_gffread/example/gff_gffread_audit.xlsx").is_file()
            )
            self.assertTrue(
                (output / "audit/reference/example/ncbi_reference_audit.xlsx").is_file()
            )
            prepared = output / "audit/reference/example/prepared_reference"
            self.assertTrue((prepared / "transcript_reference.tsv").is_file())
            self.assertTrue((prepared / "locus_reference.tsv").is_file())
            empty_sidecar = output / "audit/gff_gffread/example/unexpected_gffread_ids.tsv"
            self.assertEqual(empty_sidecar.read_text(encoding="utf-8"), "RNA_ID\n")

    def test_feature_ids_do_not_populate_rna_id_for_non_mrna_rows(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "reference.gff"
            path.write_text(
                "chr1\tRefSeq\tmRNA\t1\t100\t.\t+\t.\t"
                "ID=rna-1;Parent=gene-1;locus_tag=L1\n"
                "chr1\tRefSeq\tCDS\t1\t100\t.\t?\t0\t"
                "ID=cds-1;Parent=rna-1;protein_id=NP_1;locus_tag=L1\n",
                encoding="utf-8",
            )

            removed = scan_original_gff(path)["removed"]

            self.assertEqual(removed.loc[0, "feature_id"], "cds-1")
            self.assertEqual(removed.loc[0, "RNA_ID"], "")
            self.assertEqual(removed.loc[0, "gene_id"], "gene-1")


if __name__ == "__main__":
    unittest.main()
