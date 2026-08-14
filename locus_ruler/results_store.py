#!/usr/bin/env python3
"""Per-genome HSP-level fields kept on disk instead of in the run's memory.

Holding them all made peak memory track the genome count, which killed a
13,956-genome run partway through. Step 4 writes them here as one JSON line
per genome, and later stages read back only the line they are working on.
"""

import json
from collections.abc import Mapping
from pathlib import Path

HEAVY_FIELDS = (
    "_pieces",
    "_orphans",
    "_assembly_hsps",
    "_assembly_excluded_hsps",
    "_internal_breaks",
)


def drop_heavy(res: dict) -> None:
    for field in HEAVY_FIELDS:
        res.pop(field, None)


def split_heavy(res: dict) -> tuple[dict, dict]:
    """The same result as (light, heavy), without mutating the original."""
    heavy = {k: res[k] for k in HEAVY_FIELDS if k in res}
    light = {k: v for k, v in res.items() if k not in HEAVY_FIELDS}
    return light, heavy


class HeavyStore:
    """The heavy fields for every genome, one JSON line each, indexed by offset."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self._offsets: dict[str, int] = {}
        self._fh = None
        self._read_fh = None
        self._pos = 0

    # ── writing (Step 4)
    def open_for_write(self) -> "HeavyStore":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.path, "wb")
        self._offsets.clear()
        self._pos = 0
        return self

    def append(self, acc: str, heavy: dict) -> None:
        # The accession rides along so a lost index can be rebuilt from the file alone.
        record = {"_accession": acc, **heavy}
        line = (json.dumps(record, default=str) + "\n").encode("utf-8")
        self._fh.write(line)
        self._offsets[acc] = self._pos
        self._pos += len(line)

    def close_write(self) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None
        if self._read_fh is not None:
            self._read_fh.close()
            self._read_fh = None
        self._write_index()

    # ── reading (Step 5 and the diagnostics writers)
    def _index_path(self) -> Path:
        return self.path.with_suffix(".index.json")

    def _write_index(self) -> None:
        with open(self._index_path(), "w") as f:
            json.dump(self._offsets, f)

    def load_index(self) -> "HeavyStore":
        """Read the offset index, rebuilding it if a previous run left none."""
        index_path = self._index_path()
        if index_path.exists():
            with open(index_path) as f:
                self._offsets = json.load(f)
            return self

        self._offsets = {}
        offset = 0
        with open(self.path, "rb") as f:
            for line in f:
                acc = json.loads(line).get("_accession")
                if acc:
                    self._offsets[acc] = offset
                offset += len(line)
        return self

    def get(self, acc: str) -> dict:
        offset = self._offsets.get(acc)
        if offset is None:
            return {}
        if self._read_fh is None:
            self._read_fh = open(self.path, "rb")
        self._read_fh.seek(offset)
        record = json.loads(self._read_fh.readline())
        record.pop("_accession", None)
        return record

    def __contains__(self, acc: str) -> bool:
        return acc in self._offsets


class ResultsView(Mapping):
    """Light results in memory, heavy fields merged back in on each lookup."""

    def __init__(self, light: dict[str, dict], store: HeavyStore):
        self._light = light
        self._store = store
        self._cached_acc = None
        self._cached = None

    def __getitem__(self, acc: str) -> dict:
        if acc != self._cached_acc:
            self._cached = {**self._light[acc], **self._store.get(acc)}
            self._cached_acc = acc
        return self._cached

    def __iter__(self):
        return iter(self._light)

    def __len__(self) -> int:
        return len(self._light)


def stream_with_heavy(results: dict[str, dict], store: HeavyStore,
                      out_store: "HeavyStore | None" = None):
    """Walk *results* with heavy fields attached only while current; *out_store* keeps edits to them."""
    previous_acc = previous = None

    def release():
        if previous is None:
            return
        if out_store is not None:
            out_store.append(previous_acc, {k: previous[k]
                                            for k in HEAVY_FIELDS if k in previous})
        drop_heavy(previous)

    try:
        for acc in list(results):
            release()
            res = results[acc]
            res.update(store.get(acc))
            previous_acc, previous = acc, res
            yield acc, res
    finally:
        release()


def dump_results_json(results: Mapping, store: HeavyStore, out_path: Path) -> None:
    """Write the unchanged ruler_results.json, but a genome at a time."""
    with open(out_path, "w") as f:
        f.write("{\n")
        for i, acc in enumerate(results):
            merged = {**results[acc], **store.get(acc)}
            if i:
                f.write(",\n")
            f.write(f"{json.dumps(str(acc))}: {json.dumps(merged, indent=2, default=str)}")
        f.write("\n}\n")
