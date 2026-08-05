"""Where a Pfam search is cut, and what may be filtered afterwards."""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from locus_ruler.domain_recovery import (
    profiles_carry_ga,
    run_hmmsearch_if_configured,
    used_gathering_thresholds,
)

# Two records, cut down to the lines the check reads.
WITH_GA = """\
HMMER3/f [3.3.2]
NAME  ABC_tran
ACC   PF00005.29
GA    22.60 22.60;
//
HMMER3/f [3.3.2]
NAME  ABC_membrane
ACC   PF00664.25
GA    22.80 22.80;
//
"""

WITHOUT_GA = """\
HMMER3/f [3.3.2]
NAME  custom_profile
//
"""


class GatheringThresholdTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)

    def _hmm(self, text):
        path = self.dir / "profiles.hmm"
        path.write_text(text, encoding="utf-8")
        return path

    def test_pfam_profiles_carry_a_threshold(self):
        self.assertTrue(profiles_carry_ga(self._hmm(WITH_GA)))

    def test_a_hand_built_profile_does_not(self):
        self.assertFalse(profiles_carry_ga(self._hmm(WITHOUT_GA)))

    def test_one_profile_short_disqualifies_the_file(self):
        """hmmsearch refuses --cut_ga for the whole run, not just that model."""
        self.assertFalse(profiles_carry_ga(self._hmm(WITH_GA + WITHOUT_GA)))

    def test_an_empty_file_is_not_treated_as_qualifying(self):
        self.assertFalse(profiles_carry_ga(self._hmm("")))

    def test_a_missing_file_is_reported_not_raised(self):
        self.assertFalse(profiles_carry_ga(self.dir / "absent.hmm"))


class SearchCommandTests(unittest.TestCase):
    """The flag has to reach hmmsearch, and the marker has to record it."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)
        self.faa = self.dir / "combined_proteins.faa"
        self.faa.write_text(">a\nMKF\n", encoding="utf-8")

    def _run(self, hmm_text, dom_cfg):
        hmm = self.dir / "profiles.hmm"
        hmm.write_text(hmm_text, encoding="utf-8")
        settings = {"paths": {"root": str(self.dir)}, "tools": {}}
        cfg = {"hmm_files": [str(hmm)], **dom_cfg}
        with mock.patch("subprocess.run") as ran:
            out = run_hmmsearch_if_configured(
                settings, str(self.faa), self.dir / "work", "demo",
                domain_cfg=cfg)
        return out, ran.call_args[0][0]

    def test_pfam_profiles_are_cut_at_the_gathering_threshold(self):
        out, cmd = self._run(WITH_GA, {})
        self.assertIn("--cut_ga", cmd)
        self.assertTrue(used_gathering_thresholds(out))

    def test_the_evalue_flags_are_dropped_when_it_is_used(self):
        """Two thresholds on one search would make the looser one meaningless."""
        _, cmd = self._run(WITH_GA, {})
        self.assertNotIn("-E", cmd)
        self.assertNotIn("--domE", cmd)

    def test_a_profile_without_one_falls_back_to_evalues(self):
        out, cmd = self._run(WITHOUT_GA, {})
        self.assertNotIn("--cut_ga", cmd)
        self.assertIn("-E", cmd)
        self.assertFalse(used_gathering_thresholds(out))

    def test_it_can_be_turned_off(self):
        out, cmd = self._run(WITH_GA, {"cut_ga": False})
        self.assertNotIn("--cut_ga", cmd)
        self.assertFalse(used_gathering_thresholds(out))

    def test_an_external_domtblout_is_not_assumed_to_be_cut(self):
        """Nothing records how someone else's search was thresholded."""
        stray = self.dir / "elsewhere.domtblout"
        stray.write_text("", encoding="utf-8")
        self.assertFalse(used_gathering_thresholds(stray))

    def test_no_domtblout_at_all_reads_as_not_cut(self):
        self.assertFalse(used_gathering_thresholds(None))


if __name__ == "__main__":
    unittest.main()
