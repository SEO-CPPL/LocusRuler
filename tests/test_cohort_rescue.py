"""Coverage measured against the cohort-agreed query window."""

import unittest

from locus_ruler.cohort_rescue import (
    build_spans,
    cohort_coverage,
    consensus_span,
    describe,
)


def hit(qstart, qend, pident=55.0, bitscore=200.0):
    return {"qstart": str(qstart), "qend": str(qend),
            "pident": str(pident), "bitscore": str(bitscore)}


class ConsensusSpanTests(unittest.TestCase):
    def test_finds_the_window_the_cohort_agrees_on(self):
        hits = [hit(190, 320) for _ in range(20)]
        self.assertEqual(consensus_span(hits), (190, 320, 131))

    def test_small_cohort_yields_nothing(self):
        self.assertIsNone(consensus_span([hit(190, 320) for _ in range(9)]))

    def test_minority_outliers_do_not_widen_the_window(self):
        hits = [hit(200, 300) for _ in range(18)] + [hit(1, 320) for _ in range(2)]
        span = consensus_span(hits)
        self.assertEqual(span[:2], (200, 300))

    def test_full_length_cohort_gives_a_full_length_window(self):
        """Where everyone aligns end to end the rule must be a no-op."""
        hits = [hit(1, 322) for _ in range(30)]
        self.assertEqual(consensus_span(hits), (1, 322, 322))


class CohortCoverageTests(unittest.TestCase):
    def test_partial_hit_scores_high_against_a_domain_window(self):
        span = (190, 320, 131)
        self.assertAlmostEqual(cohort_coverage(hit(214, 320), span), 107 / 131, places=3)

    def test_disjoint_hit_scores_zero(self):
        self.assertEqual(cohort_coverage(hit(1, 100), (190, 320, 131)), 0.0)

    def test_n_terminal_hit_against_a_full_length_window_stays_low(self):
        """The short-hit case: a short hit to a long consensus."""
        self.assertAlmostEqual(cohort_coverage(hit(1, 100), (1, 298, 298)),
                               100 / 298, places=3)


class BuildSpansTests(unittest.TestCase):
    def _run(self, n_confident, qrange, tag="ANCHOR", qlen=323):
        blast = {
            f"GCF_{i:09d}.1": {tag: [hit(*qrange)]}
            for i in range(n_confident)
        }
        return build_spans(blast, {tag: qlen}, min_coverage=0.40, min_identity=40.0)

    def test_builds_from_passing_hits_only(self):
        spans = self._run(20, (1, 200))          # 200/323 = 0.62, passes
        self.assertIn("ANCHOR", spans)

    def test_anchor_with_no_passing_hits_gets_no_span(self):
        spans = self._run(20, (1, 100))          # 100/323 = 0.31, fails
        self.assertEqual(spans, {})

    def test_unknown_qlen_is_skipped(self):
        blast = {"GCF_000000001.1": {"OTHER": [hit(1, 200)]}}
        self.assertEqual(build_spans(blast, {}, 0.40, 40.0), {})

    def test_one_vote_per_genome(self):
        """Many weak HSPs in one genome must not outvote the cohort."""
        tag = "ANCHOR"
        blast = {
            f"GCF_{i:09d}.1": {tag: [hit(190, 320)]} for i in range(15)
        }
        blast["GCF_999999999.1"] = {tag: [hit(1, 322, bitscore=1.0) for _ in range(50)]}
        spans = build_spans(blast, {tag: 323}, 0.40, 40.0)
        self.assertEqual(spans[tag][:2], (190, 320))


class EndToEndShapeTests(unittest.TestCase):
    """Rescue a domain-only ortholog without rescuing a mere fragment."""

    def test_rescues_the_domain_only_ortholog(self):
        tag = "ANCHOR"
        blast = {f"G{i}": {tag: [hit(190, 320)]} for i in range(20)}
        spans = build_spans(blast, {tag: 323}, 0.40, 40.0)
        failing = hit(214, 320)                       # 107/323 = 0.33, rejected
        self.assertLess((320 - 214 + 1) / 323, 0.40)
        self.assertGreaterEqual(cohort_coverage(failing, spans[tag]), 0.40)

    def test_does_not_rescue_a_fragment_of_a_full_length_family(self):
        tag = "PROBE"
        blast = {f"G{i}": {tag: [hit(1, 298)]} for i in range(20)}
        spans = build_spans(blast, {tag: 320}, 0.40, 40.0)
        failing = hit(1, 100)                         # 100/320 = 0.31, rejected
        self.assertLess(cohort_coverage(failing, spans[tag]), 0.40)


class WhyPositionIsStillRequiredTests(unittest.TestCase):
    """A shared domain cannot separate a locus gene from a paralog elsewhere,
    so cohort coverage must never be the only evidence for a rescue."""

    def test_a_paralog_also_scores_well_on_the_domain(self):
        tag = "ANCHOR"
        blast = {f"G{i}": {tag: [hit(190, 320)]} for i in range(20)}
        span = build_spans(blast, {tag: 323}, 0.40, 40.0)[tag]

        locus_gene = hit(214, 320)        # raw 0.33
        paralog = hit(250, 320)           # raw 0.22, a shorter domain-only match

        self.assertGreaterEqual(cohort_coverage(locus_gene, span), 0.40)
        self.assertGreaterEqual(
            cohort_coverage(paralog, span), 0.40,
            "the domain alone does not discriminate -- position must",
        )


class DescribeTests(unittest.TestCase):
    def test_reports_only_anchors_whose_consensus_is_short(self):
        lines = describe({"A": (190, 320, 131), "B": (1, 322, 322)},
                         {"A": 323, "B": 323})
        self.assertEqual(len(lines), 1)
        self.assertIn("A:", lines[0])
        self.assertIn("41%", lines[0])


if __name__ == "__main__":
    unittest.main()
