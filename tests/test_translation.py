"""The CDS fallback must translate the way NCBI table 11 does."""

import unittest
from itertools import product

from locus_ruler.db_utils import CODON_TABLE, START_CODONS, revcomp, translate_dna

BASES = "TCAG"
ALL_CODONS = ["".join(c) for c in product(BASES, repeat=3)]

# NCBI translation table 11, written out independently of CODON_TABLE so the
# two have to agree rather than sharing a source.
STANDARD = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)


class CodonTableTests(unittest.TestCase):
    def test_every_codon_is_present(self):
        self.assertEqual(len(CODON_TABLE), 64)
        for codon in ALL_CODONS:
            self.assertIn(codon, CODON_TABLE)

    def test_mapping_matches_ncbi_table_11(self):
        for codon, expected in zip(ALL_CODONS, STANDARD):
            with self.subTest(codon=codon):
                self.assertEqual(CODON_TABLE[codon], expected)

    def test_the_three_stops_are_the_expected_ones(self):
        stops = {c for c, aa in CODON_TABLE.items() if aa == "*"}
        self.assertEqual(stops, {"TAA", "TAG", "TGA"})

    def test_initiators_are_the_table_11_set(self):
        self.assertEqual(
            START_CODONS,
            frozenset({"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}),
        )


class TranslationTests(unittest.TestCase):
    def test_an_atg_start_is_methionine(self):
        self.assertEqual(translate_dna("ATGAAATTTTAA"), "MKF")

    def test_a_gtg_start_is_also_methionine(self):
        """NCBI writes M here; a plain codon lookup would write V."""
        self.assertEqual(translate_dna("GTGAAATTTTAA"), "MKF")

    def test_a_ttg_start_is_also_methionine(self):
        self.assertEqual(translate_dna("TTGAAATTTTAA"), "MKF")

    def test_an_internal_gtg_stays_valine(self):
        """Only the first codon is an initiator."""
        self.assertEqual(translate_dna("ATGGTGTTTTAA"), "MVF")

    def test_a_non_initiator_start_is_left_alone(self):
        self.assertEqual(translate_dna("AAATTTTAA"), "KF")

    def test_the_terminal_stop_is_dropped(self):
        self.assertEqual(translate_dna("ATGAAATAA"), "MK")

    def test_an_internal_stop_becomes_x(self):
        self.assertEqual(translate_dna("ATGTAAAAATAA"), "MXK")

    def test_an_unknown_codon_becomes_x(self):
        self.assertEqual(translate_dna("ATGNNNAAATAA"), "MXK")

    def test_a_trailing_partial_codon_is_ignored(self):
        self.assertEqual(translate_dna("ATGAAAT"), "MK")

    def test_cds_false_leaves_the_first_codon_alone(self):
        self.assertEqual(translate_dna("GTGAAATTTTAA", cds=False), "VKF")


class RevcompTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(revcomp("ATGC"), "GCAT")

    def test_case_is_preserved(self):
        self.assertEqual(revcomp("atgc"), "gcat")

    def test_round_trip(self):
        seq = "ATGAAATTTGGGCCC"
        self.assertEqual(revcomp(revcomp(seq)), seq)


if __name__ == "__main__":
    unittest.main()
