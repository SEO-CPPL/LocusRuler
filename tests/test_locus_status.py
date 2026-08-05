import unittest

from locus_ruler.gene_state import classify_gene_state
from locus_ruler.locus_status import (
    apply_locus_coverage_floor,
    classify_locus_status,
    normalize_status_role,
    summarize_status_detail,
)


ANCHORS = [
    {"locus_tag": "core1", "role": "inner", "status_role": "CORE"},
    {"locus_tag": "core2", "role": "inner", "status_role": "CORE"},
    {"locus_tag": "assoc", "role": "aux", "status_role": "ASSOCIATED"},
    {"locus_tag": "shared", "role": "aux", "status_role": "SHARED"},
    {"locus_tag": "flank", "role": "flank_L", "status_role": "CONTEXT"},
]


class LocusStatusTests(unittest.TestCase):
    def classify(self, **states):
        return classify_locus_status(states, ANCHORS)[0]

    def test_complete(self):
        self.assertEqual(
            self.classify(
                core1="INTACT", core2="INTACT",
                assoc="INTACT", shared="INTACT",
            ),
            "COMPLETE",
        )

    def test_associated_genes_do_not_change_final_status(self):
        self.assertEqual(
            self.classify(
                core1="INTACT", core2="INTACT",
                assoc="ABSENT", shared="INTACT",
            ),
            "COMPLETE",
        )

    def test_absent_uses_only_core_genes(self):
        self.assertEqual(
            self.classify(
                core1="ABSENT", core2="PSEUDOGENE",
                assoc="INTACT", shared="INTACT",
            ),
            "ABSENT",
        )

    def test_shared_only_does_not_make_locus_present(self):
        self.assertEqual(
            self.classify(
                core1="ABSENT", core2="ABSENT",
                assoc="ABSENT", shared="INTACT",
            ),
            "ABSENT",
        )

    def test_shared_absence_makes_intact_core_conditional(self):
        self.assertEqual(
            self.classify(
                core1="INTACT", core2="INTACT",
                assoc="INTACT", shared="ABSENT",
            ),
            "CONDITIONAL",
        )

    def test_shared_decay_makes_intact_core_conditional(self):
        self.assertEqual(
            self.classify(
                core1="INTACT", core2="INTACT",
                assoc="INTACT", shared="PSEUDOGENE",
            ),
            "CONDITIONAL",
        )

    def test_core_decay_takes_priority_over_shared(self):
        self.assertEqual(
            self.classify(
                core1="INTACT", core2="ABSENT",
                assoc="INTACT", shared="INTACT",
            ),
            "DECAYED",
        )

    def test_low_coverage_decayed_is_demoted_to_absent(self):
        self.assertEqual(
            apply_locus_coverage_floor("DECAYED", 0.29, 0.30),
            "ABSENT",
        )
        self.assertEqual(
            apply_locus_coverage_floor("DECAYED", 0.30, 0.30),
            "DECAYED",
        )
        self.assertEqual(
            apply_locus_coverage_floor("COMPLETE", 0.10, 0.30),
            "COMPLETE",
        )

    def test_dna_only_is_not_gene_presence(self):
        self.assertEqual(
            classify_gene_state(0.0, 0.0, True, False),
            "ABSENT",
        )

    def test_invalid_status_role_fails(self):
        with self.assertRaises(ValueError):
            normalize_status_role("BGC", "inner")

    def test_core_role_is_required(self):
        with self.assertRaises(ValueError):
            classify_locus_status(
                {"assoc": "INTACT"},
                [{"locus_tag": "assoc", "role": "aux", "status_role": "ASSOCIATED"}],
            )

    def test_output_counts_exclude_context(self):
        _, detail = classify_locus_status(
            {
                "core1": "INTACT",
                "core2": "PSEUDOGENE",
                "assoc": "ABSENT",
                "shared": "INTACT",
                "flank": "PSEUDOGENE",
            },
            ANCHORS,
        )
        self.assertEqual(
            summarize_status_detail(detail),
            {"genes_found": 3, "genes_expected": 4, "pseudo_count": 1},
        )


if __name__ == "__main__":
    unittest.main()
