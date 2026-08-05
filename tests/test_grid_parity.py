"""loci.xlsx and cassette_matrix.xlsx share their layout constants."""

import unittest

from locus_ruler.cassette_structure import (
    MATRIX_FIXED_HEADERS,
    MATRIX_FIXED_WIDTHS,
)
from locus_ruler.duplication import _load_ruler_genome_meta
from locus_ruler.ruler import SUMMARY_COLS
from locus_ruler.writers import (
    IDENTITY_HEADERS,
    IDENTITY_WIDTHS,
    LOCI_FIXED_HEADERS,
    LOCI_FIXED_WIDTHS,
    STATUS_ALIASES,
    freeze_after_status,
)

_GRIDS = {
    "loci.xlsx": (LOCI_FIXED_HEADERS, LOCI_FIXED_WIDTHS),
    "cassette_matrix.xlsx": (MATRIX_FIXED_HEADERS, MATRIX_FIXED_WIDTHS),
}


def _loci_headers() -> list[str]:
    return LOCI_FIXED_HEADERS


def _cassette_headers() -> list[str]:
    return MATRIX_FIXED_HEADERS


class IdentityColumnTests(unittest.TestCase):
    def test_both_grids_lead_with_the_same_identity_columns(self):
        n = len(IDENTITY_HEADERS)
        self.assertEqual(_loci_headers()[:n], IDENTITY_HEADERS)
        self.assertEqual(_cassette_headers()[:n], IDENTITY_HEADERS)

    def test_assembly_level_is_on_both_grids(self):
        """The column the user found on one sheet and not the other."""
        self.assertIn("assembly_level", _loci_headers())
        self.assertIn("assembly_level", _cassette_headers())

    def test_every_identity_column_has_a_width(self):
        for name in IDENTITY_HEADERS:
            self.assertIn(name, IDENTITY_WIDTHS)

    def test_every_fixed_column_has_a_width_on_both_grids(self):
        for grid, (headers, widths) in _GRIDS.items():
            for name in headers:
                self.assertIn(name, widths, f"{name} has no width on {grid}")

    def test_shared_columns_have_the_same_width(self):
        """Side by side, the identity block must line up."""
        for name in IDENTITY_HEADERS:
            self.assertEqual(LOCI_FIXED_WIDTHS[name], MATRIX_FIXED_WIDTHS[name])

    def test_the_status_column_is_the_same_width(self):
        self.assertEqual(LOCI_FIXED_WIDTHS["status"],
                         MATRIX_FIXED_WIDTHS["locus_status"])

    def test_status_follows_the_identity_block_on_both(self):
        for headers in (_loci_headers(), _cassette_headers()):
            self.assertIn(headers[len(IDENTITY_HEADERS)], STATUS_ALIASES)


class IdentitySourceTests(unittest.TestCase):
    """A header alone proves nothing: the cassette grid reads its identity
    columns out of genome_summary.csv, which reads them out of the DB."""

    def test_the_summary_csv_carries_every_identity_column(self):
        for name in IDENTITY_HEADERS:
            self.assertIn(name, SUMMARY_COLS,
                          f"cassette_matrix reads {name} from genome_summary.csv")

    def test_the_db_loader_fills_every_identity_column(self):
        """assembly_level had a header on both grids and a value on neither."""
        import sqlite3
        import tempfile
        from pathlib import Path

        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        db = Path(tmp.name) / "t.db"
        con = sqlite3.connect(db)
        con.execute("CREATE TABLE genomes (accession TEXT, species TEXT, "
                    "strain TEXT, assembly_level TEXT)")
        con.execute("INSERT INTO genomes VALUES "
                    "('GCF_012345678.1', 'Genus species', 'S1', 'Complete Genome')")
        con.commit()
        con.close()

        meta = _load_ruler_genome_meta(str(db), ["GCF_012345678.1"])
        row = meta["GCF_012345678.1"]
        for name in IDENTITY_HEADERS:
            if name == "genome_acc":
                continue
            self.assertTrue(row.get(name),
                            f"{name} came back empty from the DB loader")


class FreezeTests(unittest.TestCase):
    def test_freeze_lands_just_past_status(self):
        # genome_acc species strain assembly_level status | coverage ...
        self.assertEqual(freeze_after_status(_loci_headers()), "F2")

    def test_both_grids_freeze_at_the_same_column(self):
        self.assertEqual(freeze_after_status(_loci_headers()),
                         freeze_after_status(_cassette_headers()))

    def test_it_stops_short_of_the_fixed_block(self):
        """Freezing every fixed column was the old behavior; status is the cut."""
        from openpyxl.utils import column_index_from_string

        for grid, (headers, _) in _GRIDS.items():
            frozen = column_index_from_string(freeze_after_status(headers)[:-1])
            self.assertLess(frozen, len(headers) + 1,
                            f"{grid} still freezes the whole fixed block")

    def test_the_status_column_itself_stays_visible(self):
        from openpyxl.utils import column_index_from_string

        for grid, (headers, _) in _GRIDS.items():
            frozen = column_index_from_string(freeze_after_status(headers)[:-1])
            status = next(i for i, n in enumerate(headers, start=1)
                          if n in STATUS_ALIASES)
            self.assertEqual(frozen, status + 1, grid)

    def test_it_follows_a_reordered_header_list(self):
        moved = ["species", "status", "genome_acc"]
        self.assertEqual(freeze_after_status(moved), "C2")

    def test_it_accepts_either_status_spelling(self):
        self.assertEqual(freeze_after_status(["a", "status"]),
                         freeze_after_status(["a", "locus_status"]))


if __name__ == "__main__":
    unittest.main()
