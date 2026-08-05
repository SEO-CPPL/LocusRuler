#!/usr/bin/env python3
"""Settings and locus config loading."""

import csv
import json
import sys
from pathlib import Path

from locus_status import normalize_status_role

try:
    import tomllib                    # Python 3.11+
except ImportError:
    try:
        import tomli as tomllib       # pip install tomli
    except ImportError:
        sys.exit("ERROR: Python 3.11+ or  pip install tomli  required.")


def load_settings(path: Path) -> dict:
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    root = Path(cfg["paths"]["root"])
    for tgt in cfg.get("targets", []):
        for k in ("db", "gff_dir", "fna_dir", "faa"):
            if k in tgt:
                p = Path(tgt[k])
                tgt[k] = str(p if p.is_absolute() else root / p)
    for k in ("output_root", "work_dir"):
        if k in cfg.get("paths", {}):
            p = Path(cfg["paths"][k])
            cfg["paths"][k] = str(p if p.is_absolute() else root / p)
    rules = cfg.get("classification", {}).get("rules_file", "")
    if rules:
        p = Path(rules)
        cfg["classification"]["rules_file"] = str(p if p.is_absolute() else root / p)
    warn_unread_sections(cfg, path)
    return cfg


# Thresholds are read from [cluster_blast] only.
_UNREAD_SECTIONS = {
    "piece": "cluster_blast",
    "gene_state": "cluster_blast",
}


def warn_unread_sections(cfg: dict, path: Path) -> None:
    for section, belongs_in in _UNREAD_SECTIONS.items():
        keys = cfg.get(section)
        if not isinstance(keys, dict) or not keys:
            continue
        print(
            f"WARNING: [{section}] in {path} is not read by LocusRuler.\n"
            f"         Move these keys into [{belongs_in}] or they have no "
            f"effect: {', '.join(sorted(keys))}",
            file=sys.stderr,
        )


def load_locus_cfg(path: Path) -> dict:
    cfg = json.loads(path.read_text())
    csv_path = path.with_name(path.stem + "_anchors.csv")
    if csv_path.exists():
        try:
            with open(csv_path, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                csv_anchors = list(reader)
            tag_to_data = {}
            for row in csv_anchors:
                d = {
                    "family":    row.get("family", ""),
                    "exception": row.get("exception", "").upper() == "TRUE",
                    "lenient":   row.get("lenient", "").upper() == "TRUE",
                    "status_role": normalize_status_role(
                        row.get("status_role", ""),
                        row.get("role", ""),
                        row.get("exception", "").upper() == "TRUE",
                    ),
                }
                if "pfam" in row:
                    d["pfam"] = row.get("pfam", "")
                if "pfam_split" in row:
                    d["pfam_split"] = row.get("pfam_split", "").upper() == "TRUE"
                tag_to_data[row["locus_tag"]] = d
            if "_auto" in cfg and "anchors" in cfg["_auto"]:
                count = 0
                existing_lts = {a.get("locus_tag") for a in cfg["_auto"]["anchors"]}
                for a in cfg["_auto"]["anchors"]:
                    lt = a.get("locus_tag")
                    if lt in tag_to_data:
                        d = tag_to_data[lt]
                        a["family"]    = d["family"]
                        a["exception"] = d["exception"]
                        a["lenient"]   = d["lenient"]
                        a["status_role"] = d["status_role"]
                        for col in ("pfam", "pfam_split"):
                            if col in d:
                                a[col] = d[col]
                        count += 1
                for row in csv_anchors:
                    lt = row.get("locus_tag", "")
                    if row.get("role", "").strip().lower() != "aux" or not lt:
                        continue
                    if lt in existing_lts:
                        continue
                    cfg["_auto"]["anchors"].append({
                        "locus_tag": lt,
                        "contig": row.get("contig", ""),
                        "start": int(row.get("start") or 0),
                        "end": int(row.get("end") or 0),
                        "strand": row.get("strand", ""),
                        "product": row.get("product", ""),
                        "family": row.get("family", ""),
                        "role": "aux",
                        "exception": row.get("exception", "").upper() == "TRUE",
                        "lenient": row.get("lenient", "").upper() == "TRUE",
                        "status_role": normalize_status_role(
                            row.get("status_role", ""),
                            "aux",
                            row.get("exception", "").upper() == "TRUE",
                        ),
                        "pfam": row.get("pfam", ""),
                        "pfam_split": row.get("pfam_split", "").upper() == "TRUE",
                    })
                    existing_lts.add(lt)
                    count += 1
                if count > 0:
                    print(f"[config] Applied {count} anchor CSV overrides "
                          f"(family/status_role/exception/lenient) from {csv_path.name}")
        except ValueError:
            raise
        except Exception as e:
            print(f"[WARN] Failed to load {csv_path.name}: {e}")
    return cfg


def save_locus_cfg(path: Path, cfg: dict) -> None:
    path.write_text(json.dumps(cfg, indent=2, ensure_ascii=False))


def get_target_cfg(settings: dict, name: str) -> dict:
    for t in settings.get("targets", []):
        if t["name"] == name:
            check_target_files(t)
            return t
    available = [t["name"] for t in settings.get("targets", [])]
    sys.exit(f"ERROR: target '{name}' not found in settings.  "
             f"Available: {available}")


def missing_target_files(target: dict) -> list[tuple[str, str, str]]:
    """The pieces a target is missing, without exiting."""
    missing = []
    for key, kind in (("db", "file"), ("gff_dir", "directory"),
                      ("fna_dir", "directory"), ("faa", "file")):
        value = target.get(key)
        if not value:
            missing.append((key, kind, "(not set in settings)"))
            continue
        path = Path(value)
        if not path.exists():
            missing.append((key, kind, str(path)))
    return missing


def check_target_files(target: dict) -> None:
    """Fail early, and say what to do about it."""
    name = target.get("name", "?")
    missing = missing_target_files(target)
    if not missing:
        return
    lines = [f"ERROR: target '{name}' is configured but its data is missing:"]
    for key, kind, where in missing:
        lines.append(f"         {key:8s} {kind} not found: {where}")
    lines += [
        "",
        "       Paths are resolved from [paths].root in settings.toml, so run",
        "       LocusRuler from the directory that file describes.",
        "",
        "       If you have not built this dataset yet:",
        f"         locus-ruler-setup --taxon \"Genus species\" "
        f"--outdir input/{name} --db-name {name}",
        "",
        "       That writes the .db, gff/, genomes/ and combined_proteins.faa",
        "       this target expects. Or run locus-ruler-wizard, which offers",
        "       to build it for you.",
    ]
    sys.exit("\n".join(lines))
