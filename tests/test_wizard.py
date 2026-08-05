"""Wizard prompts, navigation, and handoff to the pipeline."""

import argparse
import csv
import os
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from locus_ruler import discover, wizard
from locus_ruler.build_config import check_anchors_csv
from locus_ruler.locus_status import VALID_STATUS_ROLES
from locus_ruler.setup import ASSEMBLY_LEVELS, query_ncbi
from locus_ruler.wizard import (
    BACK,
    _LEVELS,
    _ask,
    _parse_levels,
    _pick_contig_cap,
    _pick_levels,
    _prompt,
    _runner,
)


class RunnerTests(unittest.TestCase):
    def test_never_relies_on_a_bare_name(self):
        first = _runner()[0]
        self.assertNotEqual(
            first, "locus-ruler",
            "a bare command name is only found when the venv is on PATH",
        )

    def test_resolves_to_something_that_exists(self):
        runner = _runner()
        self.assertTrue(Path(runner[0]).exists(), f"{runner[0]} is not a file")

    def test_module_fallback_is_importable(self):
        runner = _runner()
        if runner[0] == sys.executable:
            self.assertEqual(runner[1:], ["-m", "locus_ruler.run"])


class ColorTests(unittest.TestCase):
    """Escape codes belong on a terminal and nowhere else."""

    def test_no_escapes_when_not_a_tty(self):
        # pytest captures stdout, so isatty() is False here by construction
        self.assertNotIn("\033", discover.paint("text", "bold", "cyan"))

    def test_unknown_style_is_ignored(self):
        self.assertEqual(discover.paint("text", "chartreuse"), "text")

    def test_no_styles_returns_input(self):
        self.assertEqual(discover.paint("text"), "text")

    def test_honours_no_color(self):
        saved = discover._USE_COLOR
        try:
            discover._USE_COLOR = True
            self.assertIn("\033", discover.paint("x", "bold"))
            os.environ["NO_COLOR"] = "1"
            # The flag is computed at import, so this checks the path is switchable.
            discover._USE_COLOR = False
            self.assertNotIn("\033", discover.paint("x", "bold"))
        finally:
            discover._USE_COLOR = saved
            os.environ.pop("NO_COLOR", None)


class GroupingTests(unittest.TestCase):
    def _hit(self, contig, start, end):
        return {"locus_tag": f"t{start}", "contig": contig, "start": start,
                "end": end, "strand": "+", "product": "p", "gene_name": ""}

    def test_nearby_hits_form_one_candidate(self):
        hits = [self._hit("c1", 1000, 2000), self._hit("c1", 5000, 6000)]
        self.assertEqual(len(discover.group(hits)), 1)

    def test_distant_hits_split(self):
        hits = [self._hit("c1", 1000, 2000), self._hit("c1", 500_000, 501_000)]
        self.assertEqual(len(discover.group(hits)), 2)

    def test_different_contigs_never_merge(self):
        hits = [self._hit("c1", 1000, 2000), self._hit("c2", 1100, 2100)]
        self.assertEqual(len(discover.group(hits)), 2)


class SearchTests(unittest.TestCase):
    """One search box, three ways of naming the same gene."""

    GENES = [
        ("REF_RS00100", "c1", 1000, 2000, "+", "ABC transporter substrate-binding protein", "abcA"),
        ("REF_RS00105", "c1", 2100, 3000, "+", "short-chain dehydrogenase", "abcD"),
        ("REF_RS00110", "c1", 3100, 4000, "+", "outer membrane porin", ""),
        ("OTHER_RS999", "c1", 900_000, 901_000, "+", "hypothetical protein", ""),
    ]

    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.db = Path(self._dir.name) / "t.db"
        con = sqlite3.connect(self.db)
        con.execute("CREATE TABLE proteins (genome_acc TEXT, locus_tag TEXT, "
                    "contig TEXT, start INT, end INT, strand TEXT, "
                    "product TEXT, gene_name TEXT)")
        con.executemany("INSERT INTO proteins VALUES ('GCF_1',?,?,?,?,?,?,?)",
                        self.GENES)
        con.commit()
        con.close()
        self.addCleanup(self._dir.cleanup)

    def _find(self, term):
        return [h["locus_tag"] for h in discover.search(self.db, "GCF_1", term)]

    def test_finds_by_product_word(self):
        self.assertEqual(self._find("transporter"), ["REF_RS00100"])

    def test_finds_by_gene_name(self):
        self.assertEqual(self._find("abcD"), ["REF_RS00105"])

    def test_finds_by_exact_locus_tag(self):
        """Someone who already knows the gene should not have to describe it."""
        self.assertEqual(self._find("REF_RS00110"), ["REF_RS00110"])

    def test_locus_tag_prefix_finds_the_run(self):
        self.assertEqual(self._find("REF_RS001"),
                         ["REF_RS00100", "REF_RS00105", "REF_RS00110"])

    def test_locus_tag_search_is_case_insensitive(self):
        self.assertEqual(self._find("ref_rs00110"), ["REF_RS00110"])

    def test_a_tag_from_another_genome_is_not_a_match(self):
        self.assertEqual(discover.search(self.db, "GCF_2", "REF_RS00110"), [])


class CandidateListingTests(unittest.TestCase):
    """You have to be able to see the gene you searched for."""

    def _hit(self, tag, start, product="hypothetical protein"):
        return {"locus_tag": tag, "contig": "c1", "start": start,
                "end": start + 900, "strand": "+", "product": product,
                "gene_name": ""}

    def test_the_matched_locus_tag_is_shown(self):
        groups = discover.group([self._hit("REF_RS14250", 1000)])
        self.assertIn("REF_RS14250", "\n".join(discover.format_candidates(groups)))

    def test_every_tag_is_shown_up_to_the_limit(self):
        hits = [self._hit(f"T{i}", 1000 + i * 1000)
                for i in range(discover.SHOWN_HITS)]
        text = "\n".join(discover.format_candidates(discover.group(hits)))
        for hit in hits:
            self.assertIn(hit["locus_tag"], text)
        self.assertNotIn("more", text)

    def test_a_long_list_is_truncated_and_says_so(self):
        hits = [self._hit(f"T{i}", 1000 + i * 1000) for i in range(10)]
        text = "\n".join(discover.format_candidates(discover.group(hits)))
        self.assertIn(f"and {10 - discover.SHOWN_HITS} more", text)

    def test_tags_are_paired_with_their_own_product(self):
        hits = [self._hit("T1", 1000, "GeneA family synthase"),
                self._hit("T2", 3000, "GeneC family transporter")]
        text = "\n".join(discover.format_candidates(discover.group(hits)))
        for line in text.splitlines():
            if "T1" in line:
                self.assertIn("GeneA", line)
            if "T2" in line:
                self.assertIn("GeneC", line)


class AskTests(unittest.TestCase):
    """Enter should take the obvious answer, but only where one was offered."""

    def test_enter_takes_the_default(self):
        with mock.patch("builtins.input", return_value=""):
            self.assertEqual(_ask("? ", range(1, 3), default=1), 1)

    def test_enter_without_a_default_still_aborts(self):
        with mock.patch("builtins.input", return_value=""):
            with self.assertRaises(SystemExit):
                _ask("? ", range(1, 3))

    def test_an_explicit_number_beats_the_default(self):
        with mock.patch("builtins.input", return_value="2"):
            self.assertEqual(_ask("? ", range(1, 3), default=1), 2)

    def test_b_asks_to_go_back(self):
        for word in ("b", "B", "back"):
            with mock.patch("builtins.input", return_value=word):
                self.assertIs(_ask("? ", range(1, 3)), BACK)

    def test_back_beats_the_default(self):
        """Otherwise b would be swallowed wherever Enter means something."""
        with mock.patch("builtins.input", return_value="b"):
            self.assertIs(_ask("? ", range(1, 3), default=1), BACK)

    def test_text_prompts_take_back_only_when_asked(self):
        with mock.patch("builtins.input", return_value="b"):
            self.assertIs(_prompt("? ", allow_back=True), BACK)
            self.assertEqual(_prompt("? "), "b")


class WalkTests(unittest.TestCase):
    """Going back has to undo the answer, not just the position."""

    def _walk_with(self, scripted):
        """Drive _walk with handlers that pop their next scripted result."""
        calls = []

        def handler_for(name):
            def handler(state, ctx):
                calls.append(name)
                return scripted[name].pop(0)
            return handler

        steps = [(name, handler_for(name)) for name in scripted]
        with mock.patch.object(wizard, "_STEPS", steps):
            code = wizard._walk({"given": {}})
        return code, calls

    def test_a_straight_run_visits_each_step_once(self):
        code, calls = self._walk_with(
            {"one": ["a"], "two": ["b"], "action": [0]})
        self.assertEqual(calls, ["one", "two", "action"])
        self.assertEqual(code, 0)

    def test_back_returns_to_the_previous_step(self):
        code, calls = self._walk_with(
            {"one": ["a", "a2"], "two": [BACK, "b"], "action": [0]})
        self.assertEqual(calls, ["one", "two", "one", "two", "action"])
        self.assertEqual(code, 0)

    def test_back_from_the_first_step_re_asks_it(self):
        code, calls = self._walk_with(
            {"one": [BACK, "a"], "action": [0]})
        self.assertEqual(calls, ["one", "one", "action"])

    def test_the_action_step_returns_its_exit_code(self):
        code, _ = self._walk_with({"one": ["a"], "action": [3]})
        self.assertEqual(code, 3)

    def test_downstream_answers_are_discarded(self):
        """A new reference genome must not keep a locus chosen in the old one."""
        seen = []

        def one(state, ctx):
            return "target-b" if seen else "target-a"

        def two(state, ctx):
            seen.append(dict(state))
            return BACK if len(seen) == 1 else "locus"

        def action(state, ctx):
            seen.append(dict(state))
            return 0

        with mock.patch.object(
                wizard, "_STEPS",
                [("one", one), ("two", two), ("action", action)]):
            wizard._walk({"given": {}})
        # the second visit to `two` must not have seen a stale `two`
        self.assertNotIn("two", seen[1])


class AnchorCheckTests(unittest.TestCase):
    """A hand-edited anchor table must not fail silently."""

    HEADER = ("role,locus_tag,family,status_role,exception,lenient,"
              "pfam,pfam_split,product,contig,start,end,strand")
    ROW = "cluster_L,T1,famA,CORE,FALSE,FALSE,,FALSE,GeneA,c1,1,900,+"

    def _write(self, text):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "x_anchors.csv"
        path.write_text(text, encoding="utf-8")
        return path

    def test_a_clean_table_has_no_problems(self):
        self.assertEqual(check_anchors_csv(
            self._write(f"{self.HEADER}\n{self.ROW}\n")), [])

    def test_a_renamed_key_column_is_caught(self):
        """This is the failure that discards every curated value in silence."""
        header = self.HEADER.replace("locus_tag", "locusTag")
        problems = check_anchors_csv(self._write(f"{header}\n{self.ROW}\n"))
        self.assertTrue(any("locus_tag" in p for p in problems), problems)

    def test_an_invalid_status_role_is_caught(self):
        row = self.ROW.replace(",CORE,", ",KORE,")
        problems = check_anchors_csv(self._write(f"{self.HEADER}\n{row}\n"))
        self.assertTrue(any("KORE" in p for p in problems), problems)

    def test_a_yes_in_a_boolean_column_is_caught(self):
        """'Y' is not an error downstream, it just quietly means FALSE."""
        row = self.ROW.replace(",FALSE,FALSE,", ",Y,FALSE,")
        problems = check_anchors_csv(self._write(f"{self.HEADER}\n{row}\n"))
        self.assertTrue(any("FALSE" in p for p in problems), problems)

    def test_an_emptied_table_is_caught(self):
        self.assertTrue(check_anchors_csv(self._write(f"{self.HEADER}\n")))

    def test_a_blank_locus_tag_is_caught(self):
        row = self.ROW.replace("cluster_L,T1,", "cluster_L,,")
        problems = check_anchors_csv(self._write(f"{self.HEADER}\n{row}\n"))
        self.assertTrue(any("locus_tag" in p for p in problems), problems)

    def test_a_pfam_column_holding_a_product_name_is_caught(self):
        row = self.ROW.replace(",FALSE,,FALSE,", ",FALSE,ABC transporter,FALSE,")
        problems = check_anchors_csv(self._write(f"{self.HEADER}\n{row}\n"))
        self.assertTrue(any("PF accession" in p for p in problems), problems)

    def test_split_without_a_pfam_is_caught(self):
        """TRUE here reads as a setting, but nothing acts on it."""
        row = self.ROW.replace(",FALSE,,FALSE,", ",FALSE,,TRUE,")
        problems = check_anchors_csv(self._write(f"{self.HEADER}\n{row}\n"))
        self.assertTrue(any("pfam_split" in p for p in problems), problems)

    def test_a_pfam_with_split_is_accepted(self):
        row = self.ROW.replace(",FALSE,,FALSE,", ",FALSE,PF00005,TRUE,")
        self.assertEqual(
            check_anchors_csv(self._write(f"{self.HEADER}\n{row}\n")), [])

    def test_a_missing_file_is_reported_not_raised(self):
        problems = check_anchors_csv(Path("/nonexistent/x_anchors.csv"))
        self.assertTrue(problems)


class RoleCurationTests(unittest.TestCase):
    """Setting a role is a question, so it should not require a file format."""

    HEADER = ("role,locus_tag,family,status_role,exception,lenient,"
              "pfam,pfam_split,product,contig,start,end,strand")
    ROWS = [
        "flank_L,F1,,CONTEXT,FALSE,FALSE,,FALSE,ATPase,c1,1,900,+",
        "cluster_L,T1,famA,CORE,FALSE,FALSE,,FALSE,GeneA,c1,1000,1900,+",
        "inner,T2,famA,CORE,FALSE,FALSE,,FALSE,GeneB,c1,2000,2900,+",
        "cluster_R,T3,famB,CORE,FALSE,FALSE,,FALSE,GeneC,c1,3000,3900,+",
        "flank_R,F2,,CONTEXT,FALSE,FALSE,,FALSE,MFS,c1,4000,4900,+",
    ]

    def _write(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "x_anchors.csv"
        path.write_text("\n".join([self.HEADER, *self.ROWS]) + "\n",
                        encoding="utf-8")
        return path

    def _rows(self, path):
        with path.open(newline="", encoding="utf-8") as handle:
            return {r["locus_tag"]: r for r in csv.DictReader(handle)}

    def _run(self, path, answers):
        with mock.patch("builtins.input", side_effect=answers):
            wizard._curate_roles(path)
        return self._rows(path)

    def test_only_cluster_genes_are_asked_about(self):
        """Three prompts for three cluster genes; the flanks are not offered."""
        path = self._write()
        with mock.patch("builtins.input", side_effect=["", "", ""]) as asked:
            wizard._curate_roles(path)
        self.assertEqual(asked.call_count, 3)

    def test_enter_keeps_the_generated_value(self):
        rows = self._run(self._write(), ["", "", ""])
        self.assertEqual([rows[t]["status_role"] for t in ("T1", "T2", "T3")],
                         ["CORE", "CORE", "CORE"])

    def test_a_number_sets_the_role(self):
        rows = self._run(self._write(), ["", "", "3"])
        self.assertEqual(rows["T3"]["status_role"], "ASSOCIATED")

    def test_a_name_prefix_works_too(self):
        rows = self._run(self._write(), ["sh", "", ""])
        self.assertEqual(rows["T1"]["status_role"], "SHARED")

    def test_l_toggles_lenient_without_advancing(self):
        rows = self._run(self._write(), ["l", "", "", ""])
        self.assertEqual(rows["T1"]["lenient"], "TRUE")
        self.assertEqual(rows["T1"]["status_role"], "CORE")

    def test_x_toggles_exception(self):
        rows = self._run(self._write(), ["x", "", "", ""])
        self.assertEqual(rows["T1"]["exception"], "TRUE")

    def test_a_toggle_is_reversible(self):
        rows = self._run(self._write(), ["l", "l", "", "", ""])
        self.assertEqual(rows["T1"]["lenient"], "FALSE")

    def test_p_goes_back_a_gene(self):
        # p, not b: b means "back to the menu" everywhere else in the wizard.
        rows = self._run(self._write(), ["1", "p", "2", "", ""])
        self.assertEqual(rows["T1"]["status_role"], "SHARED")

    def test_flank_rows_are_left_alone(self):
        rows = self._run(self._write(), ["4", "4", "4"])  # 4 = IGNORE
        self.assertEqual(rows["F1"]["status_role"], "CONTEXT")
        self.assertEqual(rows["F2"]["status_role"], "CONTEXT")

    def test_the_saved_table_still_passes_the_check(self):
        path = self._write()
        self._run(path, ["2", "x", "", "3"])
        self.assertEqual(check_anchors_csv(path), [])

    def test_0_clears_both_flags_at_once(self):
        """A dedicated reset, since toggling is easy to lose track of."""
        rows = self._run(self._write(), ["l", "x", "0", "", "", ""])
        self.assertEqual(rows["T1"]["lenient"], "FALSE")
        self.assertEqual(rows["T1"]["exception"], "FALSE")

    def test_0_does_not_touch_the_role(self):
        rows = self._run(self._write(), ["2", "l", "0", "", "", ""])
        self.assertEqual(rows["T1"]["status_role"], "SHARED")

    def test_b_exits_without_visiting_the_rest(self):
        path = self._write()
        with mock.patch("builtins.input", side_effect=["b"]) as asked:
            wizard._curate_roles(path)
        self.assertEqual(asked.call_count, 1)

    def test_b_saves_edits_made_before_it(self):
        """Leaving early must not lose what was already set."""
        rows = self._run(self._write(), ["2", "b"])
        self.assertEqual(rows["T1"]["status_role"], "SHARED")

    def test_b_leaves_untouched_genes_at_their_default(self):
        rows = self._run(self._write(), ["2", "b"])
        self.assertEqual(rows["T2"]["status_role"], "CORE")

    def test_every_offered_role_is_a_valid_one(self):
        for name, _, _ in wizard._ROLE_HELP:
            self.assertIn(name, VALID_STATUS_ROLES)

    def test_every_scored_role_is_offered(self):
        # CONTEXT is absent on purpose: it scores like IGNORE and is the flank default.
        offered = {name for name, _, _ in wizard._ROLE_HELP}
        self.assertEqual(offered, VALID_STATUS_ROLES - {"CONTEXT"})


class FamilyCurationTests(unittest.TestCase):
    """Family curation walks the same rows as role curation, but sets a
    different, free-text column -- its own fixture, not a shared one."""

    HEADER = RoleCurationTests.HEADER
    ROWS = RoleCurationTests.ROWS

    def _write(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "x_anchors.csv"
        path.write_text("\n".join([self.HEADER, *self.ROWS]) + "\n",
                        encoding="utf-8")
        return path

    def _rows(self, path):
        with path.open(newline="", encoding="utf-8") as handle:
            return {r["locus_tag"]: r for r in csv.DictReader(handle)}

    def _run(self, path, answers):
        with mock.patch("builtins.input", side_effect=answers):
            wizard._curate_families(path)
        return self._rows(path)

    def test_enter_keeps_the_generated_value(self):
        rows = self._run(self._write(), ["", "", ""])
        self.assertEqual(rows["T1"]["family"], "famA")

    def test_typed_text_overwrites_it(self):
        rows = self._run(self._write(), ["renamed", "", ""])
        self.assertEqual(rows["T1"]["family"], "renamed")

    def test_p_goes_back_a_gene(self):
        rows = self._run(self._write(), ["one", "p", "two", "", ""])
        self.assertEqual(rows["T1"]["family"], "two")

    def test_b_exits_without_visiting_the_rest(self):
        path = self._write()
        with mock.patch("builtins.input", side_effect=["b"]) as asked:
            wizard._curate_families(path)
        self.assertEqual(asked.call_count, 1)

    def test_b_saves_edits_made_before_it(self):
        rows = self._run(self._write(), ["renamed", "b"])
        self.assertEqual(rows["T1"]["family"], "renamed")
        self.assertEqual(rows["T2"]["family"], "famA")


class PfamCurationTests(unittest.TestCase):
    """The Pfam columns are the ones a typo disables in silence, so the
    prompt has to reject anything that is not an accession."""

    HEADER = RoleCurationTests.HEADER
    ROWS = RoleCurationTests.ROWS

    def _write(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "x_anchors.csv"
        path.write_text("\n".join([self.HEADER, *self.ROWS]) + "\n",
                        encoding="utf-8")
        return path

    def _rows(self, path):
        with path.open(newline="", encoding="utf-8") as handle:
            return {r["locus_tag"]: r for r in csv.DictReader(handle)}

    def _run(self, path, answers):
        with mock.patch("builtins.input", side_effect=answers):
            wizard._curate_pfam(path)
        return self._rows(path)

    def test_only_cluster_genes_are_asked_about(self):
        path = self._write()
        with mock.patch("builtins.input", side_effect=["", "", ""]) as asked:
            wizard._curate_pfam(path)
        self.assertEqual(asked.call_count, 3)

    def test_an_accession_is_stored(self):
        rows = self._run(self._write(), ["PF00005", "", ""])
        self.assertEqual(rows["T1"]["pfam"], "PF00005")

    def test_several_accessions_are_normalized(self):
        """Any separator in, one canonical form out."""
        rows = self._run(self._write(), ["PF00664, PF00005", "", ""])
        self.assertEqual(rows["T1"]["pfam"], "PF00005+PF00664")

    def test_a_version_suffix_is_dropped(self):
        rows = self._run(self._write(), ["PF00005.29", "", ""])
        self.assertEqual(rows["T1"]["pfam"], "PF00005")

    def test_a_non_accession_is_refused_and_re_asked(self):
        """A product name would be stored, then silently match nothing."""
        rows = self._run(self._write(), ["ABC transporter", "PF00005", "", ""])
        self.assertEqual(rows["T1"]["pfam"], "PF00005")

    def test_s_toggles_split_without_advancing(self):
        rows = self._run(self._write(), ["s", "PF00005", "", ""])
        self.assertEqual(rows["T1"]["pfam_split"], "TRUE")
        self.assertEqual(rows["T1"]["pfam"], "PF00005")

    def test_a_toggle_is_reversible(self):
        rows = self._run(self._write(), ["s", "s", "", "", ""])
        self.assertEqual(rows["T1"]["pfam_split"], "FALSE")

    def test_0_clears_the_gene(self):
        rows = self._run(self._write(), ["PF00005", "p", "s", "0", "", ""])
        self.assertEqual(rows["T1"]["pfam"], "")
        self.assertEqual(rows["T1"]["pfam_split"], "FALSE")

    def test_enter_leaves_a_gene_alone(self):
        rows = self._run(self._write(), ["", "", ""])
        self.assertEqual(rows["T1"]["pfam"], "")

    def test_p_goes_back_a_gene(self):
        rows = self._run(self._write(), ["PF00005", "p", "PF00664", "", ""])
        self.assertEqual(rows["T1"]["pfam"], "PF00664")

    def test_b_saves_edits_made_before_it(self):
        rows = self._run(self._write(), ["PF00005", "b"])
        self.assertEqual(rows["T1"]["pfam"], "PF00005")

    def test_the_saved_table_still_passes_the_check(self):
        path = self._write()
        self._run(path, ["s", "PF00005", "", ""])
        self.assertEqual(check_anchors_csv(path), [])

    def test_split_left_without_a_pfam_is_caught_before_the_run(self):
        """The prompt allows the order s-then-accession, so it cannot refuse
        the toggle; the pre-run check is what catches the half-set gene."""
        path = self._write()
        self._run(path, ["s", "", "", ""])
        self.assertTrue(any("pfam_split" in p
                            for p in check_anchors_csv(path)))


class AuxSearchTests(unittest.TestCase):
    """A search that finds nothing, or finds too much, must not throw the
    user back to the menu with their answer lost."""

    def _hit(self, tag):
        return {"locus_tag": tag, "genome_acc": "GCF_012345678.1",
                "contig": "c1", "start": 1, "end": 900, "strand": "+",
                "product": "hypothetical protein",
                "species": "Genus species", "strain": "S1"}

    def test_no_match_asks_again_instead_of_giving_up(self):
        with mock.patch.object(wizard, "search_any_genome",
                               side_effect=[[], [self._hit("T9")]]):
            with mock.patch("builtins.input", side_effect=["zz", "gene", "1"]):
                picked = wizard._pick_aux_gene(Path("db"))
        self.assertEqual(picked["locus_tag"], "T9")

    def test_more_hits_than_fit_are_still_offered(self):
        """The old code refused outright, which made a common word a dead end."""
        many = [self._hit(f"T{i}") for i in range(wizard._AUX_SHOWN + 1)]
        with mock.patch.object(wizard, "search_any_genome", return_value=many):
            with mock.patch("builtins.input", side_effect=["gene", "1"]):
                picked = wizard._pick_aux_gene(Path("db"))
        self.assertEqual(picked["locus_tag"], "T0")

    def test_b_at_the_list_returns_to_the_search_prompt(self):
        with mock.patch.object(wizard, "search_any_genome",
                               return_value=[self._hit("T9")]):
            with mock.patch("builtins.input", side_effect=["gene", "b", ""]):
                self.assertIsNone(wizard._pick_aux_gene(Path("db")))

    def test_enter_at_the_search_prompt_cancels(self):
        with mock.patch("builtins.input", side_effect=[""]):
            self.assertIsNone(wizard._pick_aux_gene(Path("db")))


class EditorTests(unittest.TestCase):
    def test_editor_variable_wins(self):
        with mock.patch.dict(os.environ, {"EDITOR": "nano"}, clear=False):
            os.environ.pop("VISUAL", None)
            self.assertEqual(wizard._editor(), ["nano"])

    def test_visual_beats_editor(self):
        with mock.patch.dict(os.environ,
                             {"VISUAL": "vim", "EDITOR": "nano"}, clear=False):
            self.assertEqual(wizard._editor(), ["vim"])

    def test_arguments_in_the_variable_survive(self):
        """$EDITOR is a command line, not a program name."""
        with mock.patch.dict(os.environ, {"VISUAL": "code -w"}, clear=False):
            self.assertEqual(wizard._editor(), ["code", "-w"])

    def test_it_is_only_ever_a_hint(self):
        """The wizard names an editor; it never launches one."""
        self.assertFalse(hasattr(wizard, "_edit_file"))


class ConfirmTests(unittest.TestCase):
    def test_enter_can_mean_yes(self):
        with mock.patch("builtins.input", return_value=""):
            self.assertTrue(wizard._confirm("? ", default=True))

    def test_enter_still_means_no_by_default(self):
        with mock.patch("builtins.input", return_value=""):
            self.assertFalse(wizard._confirm("? "))

    def test_an_explicit_no_beats_a_yes_default(self):
        with mock.patch("builtins.input", return_value="n"):
            self.assertFalse(wizard._confirm("? ", default=True))


class DatasetDirTests(unittest.TestCase):
    """input/ pairs with output_root's output/ -- two flat top-level
    directories from a fresh checkout, not a third `data/` unrelated to
    either name.
    """

    def test_lands_under_input(self):
        self.assertEqual(wizard._dataset_dir("genus"),
                         Path("input") / "genus")

    def test_matches_the_example_targets_own_layout(self):
        # settings.example.toml's [[targets]] block points at input/example/
        self.assertEqual(wizard._dataset_dir("example"),
                         Path("input/example"))

    def test_settings_example_toml_agrees_with_this_default(self):
        """The shipped example is the one artifact that can silently drift
        from the wizard's own default -- nothing re-derives one from the
        other, so a rename on one side and not the other would not error,
        just point a first run at the wrong directory.
        """
        example = Path(__file__).resolve().parents[1] / "settings.example.toml"
        text = example.read_text(encoding="utf-8")
        self.assertIn('output_root = "output"', text)
        self.assertIn('gff_dir = "input/example/gff"', text)
        self.assertNotIn("data/example", text)
        self.assertNotIn('"outputs/', text)


class PickTargetTests(unittest.TestCase):
    """settings.example.toml ships with an `example` target already
    declared, pointing at files that do not exist until someone downloads
    them -- so "a [[targets]] block exists" and "there is something to
    analyze" are different questions, and only the second one should route
    around the download offer.
    """

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.base = Path(self._tmp.name)

    def _ready_target(self, name):
        d = self.base / name
        (d / "gff").mkdir(parents=True)
        (d / "genomes").mkdir(parents=True)
        (d / f"{name}.db").write_text("x")
        (d / "combined_proteins.faa").write_text("x")
        return {"name": name, "db": str(d / f"{name}.db"),
               "gff_dir": str(d / "gff"), "fna_dir": str(d / "genomes"),
               "faa": str(d / "combined_proteins.faa")}

    def _unready_target(self, name):
        d = self.base / name
        return {"name": name, "db": str(d / f"{name}.db"),
               "gff_dir": str(d / "gff"), "fna_dir": str(d / "genomes"),
               "faa": str(d / "combined_proteins.faa")}

    def test_a_ready_target_is_offered(self):
        settings = {"targets": [self._ready_target("real")]}
        with mock.patch("builtins.input", return_value=""):
            got = wizard._pick_target(settings, None, Path("settings.toml"))
        self.assertEqual(got, "real")

    def test_declared_but_unbuilt_is_treated_as_no_dataset(self):
        """The exact shipped-example scenario: a target exists in settings,
        nothing exists on disk."""
        settings = {"targets": [self._unready_target("example")]}
        with mock.patch.object(wizard, "_build_dataset",
                               return_value="built") as build:
            got = wizard._pick_target(settings, None, Path("settings.toml"))
        build.assert_called_once()
        self.assertEqual(got, "built")

    def test_says_which_declared_target_is_not_built(self):
        settings = {"targets": [self._unready_target("example")]}
        with mock.patch.object(wizard, "_build_dataset", return_value=None), \
             mock.patch("builtins.print") as printed:
            with self.assertRaises(SystemExit):
                wizard._pick_target(settings, None, Path("settings.toml"))
        seen = " ".join(str(c.args[0]) for c in printed.call_args_list if c.args)
        self.assertIn("example", seen)

    def test_an_unready_target_never_shows_up_in_the_menu(self):
        settings = {"targets": [self._ready_target("real"),
                                self._unready_target("example")]}
        with mock.patch("builtins.input", return_value=""), \
             mock.patch("builtins.print") as printed:
            got = wizard._pick_target(settings, None, Path("settings.toml"))
        self.assertEqual(got, "real")
        seen = " ".join(str(c.args[0]) for c in printed.call_args_list if c.args)
        self.assertIn("real", seen)
        self.assertNotIn("example", seen)

    def test_a_ready_target_is_never_routed_through_build(self):
        settings = {"targets": [self._ready_target("real")]}
        with mock.patch("builtins.input", return_value=""), \
             mock.patch.object(wizard, "_build_dataset") as build:
            wizard._pick_target(settings, None, Path("settings.toml"))
        build.assert_not_called()

    def test_an_explicit_target_skips_the_readiness_check_entirely(self):
        settings = {"targets": [self._unready_target("example")]}
        got = wizard._pick_target(settings, "example", Path("settings.toml"))
        self.assertEqual(got, "example")


class PickSettingsTests(unittest.TestCase):
    """A fresh checkout has only settings.example.toml."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._cwd = os.getcwd()
        os.chdir(self._tmp.name)
        self.addCleanup(os.chdir, self._cwd)

    def test_an_existing_settings_toml_is_used_untouched(self):
        Path("settings.toml").write_text("mine", encoding="utf-8")
        with mock.patch("builtins.input", side_effect=AssertionError(
                "should not be asked -- settings.toml already exists")):
            got = wizard._pick_settings(None)
        self.assertEqual(got, Path("settings.toml"))
        self.assertEqual(got.read_text(encoding="utf-8"), "mine")

    def test_offers_to_copy_the_example_when_nothing_else_exists(self):
        Path("settings.example.toml").write_text("[paths]\nroot = \".\"\n",
                                                  encoding="utf-8")
        with mock.patch("builtins.input", return_value=""):  # Enter = yes
            got = wizard._pick_settings(None)
        self.assertEqual(got, Path("settings.toml"))
        self.assertTrue(Path("settings.toml").exists())
        self.assertIn("root", Path("settings.toml").read_text(encoding="utf-8"))

    def test_declining_the_copy_falls_through_to_the_command_hint(self):
        Path("settings.example.toml").write_text("x", encoding="utf-8")
        with mock.patch("builtins.input", return_value="n"):
            with self.assertRaises(SystemExit) as caught:
                wizard._pick_settings(None)
        self.assertIn("cp settings.example.toml", str(caught.exception))
        self.assertFalse(Path("settings.toml").exists())

    def test_no_example_either_exits_with_the_same_hint_as_before(self):
        with self.assertRaises(SystemExit) as caught:
            wizard._pick_settings(None)
        self.assertIn("cp settings.example.toml", str(caught.exception))

    def test_an_explicit_settings_path_skips_all_of_this(self):
        explicit = Path("elsewhere.toml")
        explicit.write_text("x", encoding="utf-8")
        got = wizard._pick_settings(str(explicit))
        self.assertEqual(got, explicit)


class OutdirTests(unittest.TestCase):
    """The config (and everything step 1 writes beside it) must default to
    where run.py itself will put the results, not the current directory --
    otherwise the anchor table and the tables/diagnostics/ it produces end
    up in two unrelated places, connected only by memory of which shell the
    wizard was started from.
    """

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.output_root = Path(self._tmp.name) / "outputs" / "locus_ruler"

    def _ctx(self, outdir=None, locus_id=None):
        return {
            "settings": {"paths": {"output_root": str(self.output_root)}},
            "args": argparse.Namespace(outdir=outdir, locus_id=locus_id),
        }

    def test_matches_where_run_py_puts_the_results(self):
        # run.py: output_root = Path(settings["paths"]["output_root"])
        got = wizard._default_outdir(self._ctx(), "example", "locus_x")
        self.assertEqual(got, self.output_root / "example" / "locus_x")

    def test_write_config_uses_it_when_outdir_is_not_given(self):
        state = {"target": "example", "accession": "GCF_1.1", "term": "aero",
                 "bounds": {"flank_L": "A", "flank_R": "B"}}
        out = wizard._write_config(state, self._ctx())
        # output_root/<target>/<locus_id>/, matching run.py's locus_out
        self.assertEqual(
            out,
            self.output_root / "example" / "aero_GCF_1_1" / "aero_GCF_1_1.json")
        self.assertTrue(out.exists())

    def test_an_explicit_outdir_still_wins(self):
        explicit = Path(self._tmp.name) / "elsewhere"
        state = {"target": "example", "accession": "GCF_1.1", "term": "aero",
                 "bounds": {"flank_L": "A", "flank_R": "B"}}
        out = wizard._write_config(state, self._ctx(outdir=str(explicit)))
        self.assertEqual(out, explicit / "aero_GCF_1_1.json")


class StepActionWiringTests(unittest.TestCase):
    """Regression: choosing [2] once raised TypeError at runtime --
    `_step_action` called `_curate_then_run` without the `db` it needs for
    aux search, and no test drove that call site, so a live database was
    the first thing to notice.
    """

    def _ctx(self):
        return {
            "settings": {"targets": [{"name": "example", "db": "x.db"}]},
            "settings_path": Path("settings.toml"),
            "args": argparse.Namespace(locus_id=None, outdir=None),
        }

    def test_curate_then_run_receives_a_db_path(self):
        with mock.patch.object(wizard, "_write_config",
                               return_value=Path("out.json")), \
             mock.patch.object(wizard, "_runner", return_value=["lr"]), \
             mock.patch.object(wizard, "_ask", return_value=2), \
             mock.patch.object(wizard, "_db_for", return_value=Path("x.db")), \
             mock.patch.object(wizard, "_curate_then_run") as curate:
            wizard._step_action({"target": "example"}, self._ctx())
        self.assertEqual(len(curate.call_args.args), 3,
                         "_curate_then_run is called with the wrong arity")
        self.assertEqual(curate.call_args.args[2], Path("x.db"))

    def test_stop_at_anchors_does_not_need_a_db(self):
        """Its arity is different on purpose; this pins that difference."""
        with mock.patch.object(wizard, "_write_config",
                               return_value=Path("out.json")), \
             mock.patch.object(wizard, "_runner", return_value=["lr"]), \
             mock.patch.object(wizard, "_ask", return_value=3), \
             mock.patch.object(wizard, "_stop_at_anchors") as stop:
            wizard._step_action({"target": "example"}, self._ctx())
        self.assertEqual(len(stop.call_args.args), 2)


class GivenTests(unittest.TestCase):
    def test_a_flag_is_used_once_then_forgotten(self):
        """Coming back to a question must ask it, not reinstate the flag."""
        ctx = {"given": {"target": "example"}}
        self.assertEqual(wizard._given(ctx, "target"), "example")
        self.assertIsNone(wizard._given(ctx, "target"))


class AssemblyLevelTests(unittest.TestCase):
    """The wizard's menu and the downloader's flag have to mean the same thing."""

    _NAMES = [name for name, _ in _LEVELS]

    def test_the_menu_matches_the_downloader(self):
        self.assertEqual(self._NAMES, list(ASSEMBLY_LEVELS))

    def test_enter_keeps_every_level(self):
        with mock.patch("builtins.input", return_value=""):
            self.assertEqual(_pick_levels(), list(ASSEMBLY_LEVELS))

    def test_picking_one_returns_one_level(self):
        with mock.patch("builtins.input", return_value="1"):
            self.assertEqual(_pick_levels(), ["complete"])

    def test_several_numbers_are_all_kept(self):
        with mock.patch("builtins.input", return_value="1 3"):
            self.assertEqual(_pick_levels(), ["complete", "scaffold"])

    def test_a_non_adjacent_combination_is_allowed(self):
        """Complete plus raw contigs, skipping the middle, is a real request."""
        self.assertEqual(_parse_levels("1 4", self._NAMES),
                         ["complete", "contig"])

    def test_commas_work_too(self):
        self.assertEqual(_parse_levels("1,2", self._NAMES),
                         ["complete", "chromosome"])

    def test_names_work_instead_of_numbers(self):
        self.assertEqual(_parse_levels("scaffold complete", self._NAMES),
                         ["complete", "scaffold"])

    def test_an_unambiguous_prefix_is_enough(self):
        self.assertEqual(_parse_levels("comp scaf", self._NAMES),
                         ["complete", "scaffold"])

    def test_the_result_is_always_in_menu_order(self):
        self.assertEqual(_parse_levels("4 2", self._NAMES),
                         ["chromosome", "contig"])

    def test_duplicates_collapse(self):
        self.assertEqual(_parse_levels("1 1 complete", self._NAMES),
                         ["complete"])

    def test_a_typo_rejects_the_whole_line(self):
        """Silently dropping a token would narrow the download unnoticed."""
        self.assertEqual(_parse_levels("1 scafold", self._NAMES), [])

    def test_an_out_of_range_number_is_rejected(self):
        self.assertEqual(_parse_levels("1 9", self._NAMES), [])

    def test_a_typo_is_re_asked(self):
        with mock.patch("builtins.input", side_effect=["nonsense", "2"]):
            self.assertEqual(_pick_levels(), ["chromosome"])

    def test_a_narrowed_query_asks_ncbi_for_less(self):
        """The level choice has to reach the URL, not just the printout."""
        seen = {}

        class _Resp:
            def __enter__(self_inner):
                return self_inner

            def __exit__(self_inner, *a):
                return False

            def read(self_inner):
                return b'{"reports": []}'

        def _capture(req, timeout=None):
            seen["url"] = req.full_url
            return _Resp()

        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch("urllib.request.urlopen", _capture):
                query_ncbi("Genus", 0, set(), Path(tmp), ["complete"])
        self.assertIn("assembly_level=complete_genome", seen["url"])
        self.assertNotIn("assembly_level=contig", seen["url"])
        self.assertNotIn("assembly_level=scaffold", seen["url"])

    def test_the_default_query_asks_for_all_four(self):
        seen = {}

        class _Resp:
            def __enter__(self_inner):
                return self_inner

            def __exit__(self_inner, *a):
                return False

            def read(self_inner):
                return b'{"reports": []}'

        def _capture(req, timeout=None):
            seen["url"] = req.full_url
            return _Resp()

        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch("urllib.request.urlopen", _capture):
                query_ncbi("Genus", 0, set(), Path(tmp))
        for ncbi_name in ASSEMBLY_LEVELS.values():
            self.assertIn(f"assembly_level={ncbi_name}", seen["url"])


class ContigCapTests(unittest.TestCase):
    def test_enter_takes_the_shown_default(self):
        with mock.patch("builtins.input", return_value=""):
            self.assertEqual(_pick_contig_cap(100), 100)

    def test_zero_is_accepted_as_no_limit(self):
        with mock.patch("builtins.input", return_value="0"):
            self.assertEqual(_pick_contig_cap(), 0)

    def test_a_number_is_taken_verbatim(self):
        with mock.patch("builtins.input", return_value="5"):
            self.assertEqual(_pick_contig_cap(), 5)

    def test_nonsense_is_re_asked_not_silently_defaulted(self):
        with mock.patch("builtins.input", side_effect=["lots", "-3", "12"]):
            self.assertEqual(_pick_contig_cap(), 12)


if __name__ == "__main__":
    unittest.main()
