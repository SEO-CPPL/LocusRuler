import unittest

from locus_ruler.fragmentation import (
    ASSEMBLY_SPLIT_CANDIDATE,
    classify_fragmentation,
)
from locus_ruler.ruler import _consolidate_piece


class FragmentationTests(unittest.TestCase):
    def test_non_wrapping_piece_preserves_contig_length_for_edge_checks(self):
        piece = _consolidate_piece([{
            "sseqid": "contig_a",
            "sstart": "338",
            "send": "1372",
            "qstart": "9303",
            "qend": "10335",
            "length": "1033",
            "pident": "81.7",
            "slen": "33149",
        }])
        self.assertEqual(piece["_contig_len"], 33149)
        self.assertNotIn("_wraps_origin", piece)

    def test_complementary_edge_pieces_are_assembly_split_candidates(self):
        pieces = [
            {
                "sseqid": "contig_a",
                "sstart": "9140",
                "send": "161",
                "qstart": "1",
                "qend": "8978",
                "_contig_len": 19527,
            },
            {
                "sseqid": "contig_b",
                "sstart": "338",
                "send": "1372",
                "qstart": "9303",
                "qend": "10335",
                "_contig_len": 33149,
            },
        ]
        result = classify_fragmentation(
            pieces,
            assembly_level="Scaffold",
            coverage=0.969,
            settings={
                "q_comp_threshold": 0.70,
                "q_overlap_threshold": 0.30,
                "edge_margin_bp": 500,
            },
        )
        self.assertEqual(result["fragmentation_type"], ASSEMBLY_SPLIT_CANDIDATE)
        self.assertEqual(result["contig_edge_support"], "Y")


if __name__ == "__main__":
    unittest.main()
