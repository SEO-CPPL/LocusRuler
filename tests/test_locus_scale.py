"""Windows derive from the reference locus unless pinned explicitly."""

import unittest

from locus_ruler.locus_scale import (
    DEFAULT_REFERENCE_BP,
    DERIVE,
    describe,
    reference_length,
    resolve,
)


class ReferenceLengthTests(unittest.TestCase):
    def test_reads_the_recorded_cluster_length(self):
        self.assertEqual(reference_length({"_auto": {"cluster_ref_bp": 10335}}), 10335)

    def test_falls_back_when_the_field_is_absent(self):
        for cfg in ({}, {"_auto": {}}, {"_auto": {"cluster_ref_bp": 0}},
                    {"_auto": {"cluster_ref_bp": None}}):
            with self.subTest(cfg=cfg):
                self.assertEqual(reference_length(cfg), DEFAULT_REFERENCE_BP)

    def test_garbage_does_not_raise(self):
        self.assertEqual(reference_length({"_auto": {"cluster_ref_bp": "n/a"}}),
                         DEFAULT_REFERENCE_BP)


class ResolveTests(unittest.TestCase):
    def test_negative_derives_from_the_reference(self):
        self.assertEqual(resolve(DERIVE, 10_000, 0.50), 5_000)
        self.assertEqual(resolve(DERIVE, 10_000, 0.20), 2_000)

    def test_explicit_values_are_honored(self):
        self.assertEqual(resolve(7_500, 10_000, 0.50), 7_500)
        self.assertEqual(resolve(0, 10_000, 0.50), 0)

    def test_missing_configuration_derives(self):
        self.assertEqual(resolve(None, 10_000, 0.50), 5_000)

    def test_minimum_is_respected(self):
        self.assertEqual(resolve(DERIVE, 100, 0.50, minimum=500), 500)

    def test_a_derived_window_is_a_constant_fraction(self):
        """A fixed bp figure is not, which is why the fraction exists."""
        small, large = 3_992, 35_779
        self.assertAlmostEqual(resolve(DERIVE, small, 0.50) / small, 0.50, places=2)
        self.assertAlmostEqual(resolve(DERIVE, large, 0.50) / large, 0.50, places=2)
        self.assertGreater(5000 / small, 1.0)
        self.assertLess(5000 / large, 0.2)

    def test_a_distant_hit_still_fits_the_derived_window(self):
        self.assertGreaterEqual(resolve(DERIVE, 10_335, 0.50), 4_813)


class DescribeTests(unittest.TestCase):
    def test_derived_windows_report_their_share(self):
        line = describe("gap", DERIVE, 5_168, 10_335)
        self.assertIn("5,168 bp", line)
        self.assertIn("50%", line)
        self.assertIn("10,335 bp reference locus", line)

    def test_pinned_windows_say_so(self):
        self.assertIn("fixed by settings", describe("gap", 5_000, 5_000, 10_335))


if __name__ == "__main__":
    unittest.main()
