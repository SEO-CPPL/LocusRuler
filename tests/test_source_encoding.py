"""Source files stay UTF-8, free of CP949 round-trip damage."""

import re
import unittest
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent / "locus_ruler"

# What UTF-8 punctuation turns into after a CP949 round-trip:
DAMAGE_SIGNATURES = {
    "box drawing eaten by CP949": re.compile("\\?\u0080"),
    "orphan U+0080 control char": re.compile("\u0080"),
    "multiplication sign eaten by CP949": re.compile("\ud69e"),
    "punctuation collapsed to '??'": re.compile(r"\?\?"),
    "stray Hangul from a split sequence": re.compile(
        "[\ubb5b\ubb65\ubb68\ubb6a\uace3\ubd3a]"
    ),
}


def _python_sources():
    return sorted(PKG.glob("*.py"))


class SourceEncodingTests(unittest.TestCase):
    def test_every_module_is_valid_utf8(self):
        for path in _python_sources():
            with self.subTest(module=path.name):
                try:
                    path.read_bytes().decode("utf-8")
                except UnicodeDecodeError as exc:
                    self.fail(f"{path.name} is not valid UTF-8: {exc}")

    def test_no_mojibake_signature_in_sources(self):
        for path in _python_sources():
            text = path.read_text(encoding="utf-8")
            for lineno, line in enumerate(text.splitlines(), 1):
                for label, pattern in DAMAGE_SIGNATURES.items():
                    with self.subTest(module=path.name, line=lineno, damage=label):
                        self.assertIsNone(
                            pattern.search(line),
                            f"{path.name}:{lineno} looks like {label}: "
                            f"{line.strip()[:90]!r}",
                        )

    def test_no_control_characters_outside_tab_and_newline(self):
        allowed = {"\t", "\n", "\r"}
        for path in _python_sources():
            text = path.read_text(encoding="utf-8")
            bad = {
                ch for ch in text
                if ch not in allowed and (ord(ch) < 0x20 or 0x7F <= ord(ch) <= 0x9F)
            }
            with self.subTest(module=path.name):
                self.assertEqual(
                    bad, set(),
                    f"{path.name} contains control characters: "
                    f"{sorted(hex(ord(c)) for c in bad)}",
                )


class GeneStateDocstringTests(unittest.TestCase):
    """The decision table documents thresholds the published results depend on,
    so it must state its operators explicitly, not a damaged placeholder."""

    def test_decision_table_states_its_operators(self):
        from locus_ruler.gene_state import classify_gene_state

        doc = classify_gene_state.__doc__ or ""
        self.assertIn("cov >= intact_cov AND pid >= intact_pid", doc)
        self.assertIn("min_cov <= cov < intact_cov", doc)
        for state in ("PSEUDOGENE", "INTACT", "DIVERGENT", "ABSENT"):
            self.assertIn(f"-> {state}", doc)


if __name__ == "__main__":
    unittest.main()
