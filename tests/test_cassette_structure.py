import unittest

from locus_ruler.cassette_structure import (
    CASSETTE_SPAN_MULTIPLE,
    _accepted_pieces,
    _colocated_pieces,
    _display_genes,
    _display_settings,
    _ordered_genes,
    _signatures,
    _style_for_gene,
    _structure_id,
)


class CassetteStructureTests(unittest.TestCase):
    def setUp(self):
        self.piece_rows = [
            {
                "genome_acc": "genome",
                "locus_id": "locus",
                "piece_idx": "0",
                "contig": "contig_a",
                "strand": "+",
                "q_start": "1",
                "q_end": "500",
            },
            {
                "genome_acc": "genome",
                "locus_id": "locus",
                "piece_idx": "1",
                "contig": "contig_b",
                "strand": "-",
                "q_start": "501",
                "q_end": "800",
            },
        ]
        self.diagnostics = [
            {
                "piece_idx": "0",
                "zone": "flank_left",
                "gene_contig": "contig_a",
                "gene_start": "1",
                "gene_end": "90",
                "locus_tag": "left_outer",
                "family": "left_flank",
                "state": "INTACT",
            },
            {
                "piece_idx": "0",
                "zone": "cluster",
                "gene_contig": "contig_a",
                "gene_start": "100",
                "gene_end": "200",
                "locus_tag": "gene_a",
                "family": "family_a",
                "state": "INTACT",
            },
            {
                "piece_idx": "1",
                "zone": "cluster",
                "gene_contig": "contig_b",
                "gene_start": "900",
                "gene_end": "990",
                "locus_tag": "gene_b",
                "family": "family_b",
                "state": "INTACT",
            },
            {
                "piece_idx": "1",
                "zone": "cluster",
                "gene_contig": "contig_b",
                "gene_start": "700",
                "gene_end": "800",
                "locus_tag": "gene_c",
                "family": "",
                "product": "Novel protein",
                "state": "",
            },
            {
                "piece_idx": "1",
                "zone": "flank_right",
                "gene_contig": "contig_b",
                "gene_start": "600",
                "gene_end": "690",
                "locus_tag": "right_outer",
                "family": "right_flank",
                "state": "INTACT",
            },
            {
                "piece_idx": "1",
                "zone": "flank_right",
                "gene_contig": "contig_b",
                "gene_start": "500",
                "gene_end": "590",
                "locus_tag": "untyped_local_context",
                "family": "",
                "label": "[unassigned]",
                "state": "",
            },
            {
                "piece_idx": "not_accepted",
                "zone": "cluster",
                "gene_contig": "orphan",
                "gene_start": "1",
                "gene_end": "100",
                "locus_tag": "cargo_like_hit",
                "family": "known_cargo",
                "state": "INTACT",
            },
        ]

    def test_only_main_locus_ruler_pieces_define_membership(self):
        pieces = _accepted_pieces(self.piece_rows, "locus")["genome"]
        genes = _ordered_genes(self.diagnostics, pieces)
        self.assertEqual(
            [gene["locus_tag"] for gene in genes],
            [
                "left_outer",
                "gene_a",
                "gene_b",
                "gene_c",
                "right_outer",
                "untyped_local_context",
            ],
        )
        self.assertNotIn("cargo_like_hit", {gene["locus_tag"] for gene in genes})

    def test_sibling_piece_is_included_in_query_order(self):
        pieces = _accepted_pieces(self.piece_rows, "locus")["genome"]
        genes = _ordered_genes(self.diagnostics, pieces)
        structure, annotation, state, assembly = _signatures(genes)
        self.assertEqual(structure, "family_a>family_b>unassigned")
        self.assertEqual(
            assembly,
            "family_a>||>family_b>unassigned",
        )
        self.assertIn("unassigned[novel_protein]", annotation)
        self.assertIn("family_b:INTACT", state)

    def test_display_keeps_typed_flanks_but_excludes_untyped_piece_context(self):
        pieces = _accepted_pieces(self.piece_rows, "locus")["genome"]
        genes = _ordered_genes(self.diagnostics, pieces)
        displayed = _display_genes(genes)
        self.assertEqual(
            [gene["locus_tag"] for gene in displayed],
            ["left_outer", "gene_a", "gene_b", "gene_c", "right_outer"],
        )
        self.assertEqual(
            [gene["locus_tag"] for gene in displayed if gene["is_cassette_member"] == "Y"],
            ["gene_a", "gene_b", "gene_c"],
        )

    def test_single_piece_display_keeps_untyped_bracketing_flank_context(self):
        pieces = _accepted_pieces(self.piece_rows[:1], "locus")["genome"]
        diagnostics = [
            self.diagnostics[0],
            self.diagnostics[1],
            {
                "piece_idx": "0",
                "zone": "flank_right",
                "gene_contig": "contig_a",
                "gene_start": "210",
                "gene_end": "300",
                "locus_tag": "untyped_right_flank",
                "family": "",
                "label": "[unassigned]",
                "state": "",
            },
        ]
        displayed = _display_genes(_ordered_genes(diagnostics, pieces))
        self.assertEqual(
            [gene["locus_tag"] for gene in displayed],
            ["left_outer", "gene_a", "untyped_right_flank"],
        )

    def test_product_style_is_display_only_fallback_for_untyped_flank(self):
        display = _display_settings({
            "cassette_display": {
                "product_styles": [
                    {"pattern": "(?i)^example flank$", "style": "right_flank"}
                ],
                "family_styles": {
                    "right_flank": {"color": "#123456", "label": ""}
                },
            }
        })
        gene = {
            "zone": "flank_right",
            "family": "",
            "label": "[unassigned]",
            "product": "Example flank",
        }
        style = _style_for_gene(gene, display)
        self.assertEqual(style["key"], "right_flank")
        self.assertEqual(style["color"], "#123456")

    def test_structure_id_is_independent_of_assembly_boundary(self):
        signature = "family_a>family_b>unassigned"
        self.assertEqual(_structure_id(signature), _structure_id(signature))
        self.assertTrue(_structure_id(signature).startswith("CS_"))

    def test_adjacent_non_intact_fragments_collapse_but_intact_copies_remain(self):
        pieces = _accepted_pieces(self.piece_rows, "locus")["genome"]
        fragments = [
            {
                "piece_idx": "0",
                "zone": "cluster",
                "gene_contig": "contig_a",
                "gene_start": str(start),
                "gene_end": str(start + 90),
                "locus_tag": f"fragment_{start}",
                "family": "one_family",
                "state": state,
            }
            for start, state in ((100, "PSEUDOGENE"), (210, "PSEUDOGENE"))
        ]
        fragment_signature = _signatures(_ordered_genes(fragments, pieces))[0]
        self.assertEqual(fragment_signature, "one_family")
        for row in fragments:
            row["state"] = "INTACT"
        copy_signature = _signatures(_ordered_genes(fragments, pieces))[0]
        self.assertEqual(copy_signature, "one_family>one_family")

    def test_two_intact_cds_under_one_alignment_are_one_split_gene(self):
        """Both CDS complete, so both read INTACT, yet the anchor aligned once
        across the pair: the annotation split a gene rather than duplicating it.
        A real duplication aligns twice and disagrees on coverage or identity."""
        pieces = _accepted_pieces(self.piece_rows, "locus")["genome"]

        def halves(cov_a, pid_a, cov_b, pid_b):
            return [
                {
                    "piece_idx": "0",
                    "zone": "cluster",
                    "gene_contig": "contig_a",
                    "gene_start": str(start),
                    "gene_end": str(start + 90),
                    "locus_tag": f"half_{start}",
                    "family": "one_family",
                    "state": "INTACT",
                    "tblastn_cov": cov,
                    "tblastn_pid": pid,
                }
                for start, cov, pid in ((100, cov_a, pid_a), (210, cov_b, pid_b))
            ]

        shared = _signatures(_ordered_genes(halves("0.966", "77.6", "0.966", "77.6"), pieces))[0]
        self.assertEqual(shared, "one_family")

        distinct = _signatures(_ordered_genes(halves("0.98", "81.2", "0.72", "64.5"), pieces))[0]
        self.assertEqual(distinct, "one_family>one_family")


# The span is now resolved per locus, so the tests pin one explicitly.
SPAN = 60_000


def _piece(contig, lo, hi, covered=1, edge="N"):
    return {
        "contig": contig,
        "s_lo": str(lo),
        "s_hi": str(hi),
        "n_covered_genes": str(covered),
        "aligned_bp": str(hi - lo),
        "contig_edge_proximal": edge,
    }


class ColocatedPieceTests(unittest.TestCase):
    """Cassette membership is subject-space: only pieces physically near the
    body count, so a lenient rescue elsewhere in the genome cannot inflate the
    cassette signature."""

    def setUp(self):
        self.body = _piece("contig_a", 100_000, 140_000, covered=8)

    def test_single_piece_passes_through(self):
        kept, dropped = _colocated_pieces({"0": self.body}, SPAN)
        self.assertEqual(sorted(kept), ["0"])
        self.assertEqual(dropped, [])

    def test_nearby_same_contig_piece_joins_and_distant_one_is_dropped(self):
        near = _piece("contig_a", 150_000, 155_000)
        far = _piece("contig_a", 2_000_000, 2_005_000)
        kept, dropped = _colocated_pieces(
            {"0": self.body, "1": near, "2": far}, SPAN
        )
        self.assertEqual(sorted(kept), ["0", "1"])
        self.assertEqual(dropped, ["2"])

    def test_gap_boundary_is_inclusive(self):
        edge_of_span = _piece(
            "contig_a", 140_000 + SPAN, 145_000 + SPAN
        )
        just_beyond = _piece(
            "contig_a",
            140_001 + SPAN,
            145_001 + SPAN,
        )
        kept, _ = _colocated_pieces({"0": self.body, "1": edge_of_span}, SPAN)
        self.assertEqual(sorted(kept), ["0", "1"])
        kept, dropped = _colocated_pieces({"0": self.body, "1": just_beyond}, SPAN)
        self.assertEqual(sorted(kept), ["0"])
        self.assertEqual(dropped, ["1"])

    def test_other_contig_joins_only_when_both_sides_touch_a_contig_edge(self):
        other = _piece("contig_b", 1, 5_000)
        kept, dropped = _colocated_pieces({"0": self.body, "1": other}, SPAN)
        self.assertEqual(sorted(kept), ["0"])
        self.assertEqual(dropped, ["1"])

        split_body = _piece("contig_a", 100_000, 140_000, covered=8, edge="Y")
        split_other = _piece("contig_b", 1, 5_000, edge="Y")
        kept, dropped = _colocated_pieces({"0": split_body, "1": split_other}, SPAN)
        self.assertEqual(sorted(kept), ["0", "1"])
        self.assertEqual(dropped, [])

    def test_body_is_the_piece_covering_most_reference_genes(self):
        small_first = _piece("contig_a", 1_000, 2_000, covered=1)
        big_elsewhere = _piece("contig_b", 500_000, 560_000, covered=9)
        kept, dropped = _colocated_pieces(
            {"0": small_first, "1": big_elsewhere}, SPAN
        )
        self.assertEqual(sorted(kept), ["1"])
        self.assertEqual(dropped, ["0"])


if __name__ == "__main__":
    unittest.main()
