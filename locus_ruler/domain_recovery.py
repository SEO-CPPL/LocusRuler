#!/usr/bin/env python3
"""Optional Pfam architecture recovery."""

import csv
import re
import subprocess
from pathlib import Path
from typing import Optional


def resolve_domtblout_path(
    settings: dict,
    target_name: str,
    domain_cfg: Optional[dict] = None,
) -> Optional[Path]:
    """Return an explicitly configured domtblout path."""
    dom_cfg = domain_cfg if domain_cfg is not None else settings.get("domain_recovery", {})
    raw = dom_cfg.get("pfam_domtblout", "")
    if raw:
        p = Path(raw)
        if not p.is_absolute():
            root = Path(settings.get("paths", {}).get("root", "."))
            p = root / p
        return p if p.exists() else None
    return None


def run_hmmsearch_if_configured(
    settings: dict,
    target_faa: str,
    work_dir: Path,
    locus_id: str,
    cpu: int = 1,
    domain_cfg: Optional[dict] = None,
    pfam_accessions: Optional[list[str]] = None,
) -> Optional[Path]:
    """Run hmmsearch for configured HMM files and return domtblout path."""
    dom_cfg = domain_cfg if domain_cfg is not None else settings.get("domain_recovery", {})
    hmm_files = list(dom_cfg.get("hmm_files", []) or [])
    pfam_accessions = sorted(set(pfam_accessions or []))
    if (not hmm_files and not pfam_accessions) or not target_faa or not Path(target_faa).exists():
        return None

    root = Path(settings.get("paths", {}).get("root", "."))
    out_dir = work_dir / locus_id / "domain_recovery"
    out_dir.mkdir(parents=True, exist_ok=True)

    resolved: list[Path] = []
    for raw in hmm_files:
        p = Path(raw)
        if not p.is_absolute():
            p = root / p
        if p.exists():
            resolved.append(p)
        else:
            print(f"[domain_recovery] WARN: HMM file not found: {p}")
    pfam_a_raw = dom_cfg.get("pfam_a_hmm", "")
    if pfam_accessions and pfam_a_raw:
        pfam_a = Path(pfam_a_raw)
        if not pfam_a.is_absolute():
            pfam_a = root / pfam_a
        if pfam_a.exists():
            fetched = out_dir / "pfam_selected.hmm"
            keyfile = out_dir / "pfam_accessions.txt"
            key_map = _pfam_accession_to_name(pfam_a, pfam_accessions)
            keys = [key_map.get(acc, acc) for acc in pfam_accessions]
            keyfile.write_text("\n".join(keys) + "\n")
            hmmfetch = settings.get("tools", {}).get("hmmfetch", "hmmfetch")
            cmd = [hmmfetch, "-f", "-o", str(fetched), str(pfam_a), str(keyfile)]
            print(f"[domain_recovery] fetching {len(pfam_accessions)} Pfam profile(s) from {pfam_a}")
            try:
                subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                if fetched.exists() and fetched.stat().st_size > 0:
                    resolved.append(fetched)
            except subprocess.CalledProcessError:
                print("[domain_recovery] WARN: hmmfetch failed; check pfam_a_hmm and PF accessions")
        else:
            print(f"[domain_recovery] WARN: Pfam-A HMM not found: {pfam_a}")

    if not resolved:
        return None

    combined = out_dir / "configured_domains.hmm"
    domtbl = out_dir / "configured_domains.domtblout"
    done = out_dir / "configured_domains.done"
    if domtbl.exists() and done.exists():
        return domtbl

    with open(combined, "wb") as fout:
        for p in resolved:
            fout.write(p.read_bytes())
            fout.write(b"\n")

    hmmsearch = settings.get("tools", {}).get("hmmsearch", "hmmsearch")
    cmd = [
        hmmsearch,
        "--cpu", str(max(1, int(cpu))),
        "--domtblout", str(domtbl),
        "--noali",
    ]
    marker = domtbl.with_suffix(".cut_ga")
    if marker.exists():
        marker.unlink()
    want_ga = bool(dom_cfg.get("cut_ga", True))
    if want_ga and profiles_carry_ga(combined):
        cmd.append("--cut_ga")
        marker.write_text("ok\n")
        print("[domain_recovery] cutting at Pfam gathering thresholds")
    else:
        if want_ga:
            print("[domain_recovery] WARN: a profile carries no GA line, so "
                  "E-value thresholds are used instead")
        # HMMER sequence E-values scale with database size.
        report_evalue = str(dom_cfg.get("report_evalue", dom_cfg.get("evalue", "1e-5")))
        report_dom_evalue = str(dom_cfg.get("report_domain_evalue", report_evalue))
        cmd += ["-E", report_evalue, "--domE", report_dom_evalue]
    cmd += [str(combined), str(target_faa)]
    print(f"[domain_recovery] running hmmsearch for {len(resolved)} HMM profile file(s)")
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    done.write_text("ok\n")
    return domtbl


def profiles_carry_ga(hmm_path: Path) -> bool:
    """Whether every profile in a file has a curated gathering threshold.

    hmmsearch refuses --cut_ga outright when one profile lacks the line, so a
    hand-built HMM mixed in with Pfam's has to send the search the other way.
    """
    names = cutoffs = 0
    try:
        with open(hmm_path) as f:
            for line in f:
                if line.startswith("NAME "):
                    names += 1
                elif line.startswith("GA "):
                    cutoffs += 1
    except OSError:
        return False
    return names > 0 and cutoffs == names


def used_gathering_thresholds(domtblout: Optional[Path]) -> bool:
    """Whether the search behind `domtblout` was cut at Pfam's own thresholds.

    An externally supplied domtblout carries no marker, so it reads as False
    and keeps the E-value filtering it was written under.
    """
    return bool(domtblout) and domtblout.with_suffix(".cut_ga").exists()


def _pfam_accession_to_name(pfam_a: Path, accessions: list[str]) -> dict[str, str]:
    """Map PF accessions to HMM NAME keys for hmmfetch."""
    wanted = {acc.split(".")[0] for acc in accessions}
    found: dict[str, str] = {}
    current_name = ""
    with open(pfam_a) as f:
        for line in f:
            if line.startswith("NAME"):
                parts = line.split()
                current_name = parts[1] if len(parts) > 1 else ""
            elif line.startswith("ACC") and current_name:
                parts = line.split()
                if len(parts) > 1:
                    acc = parts[1].split(".")[0]
                    if acc in wanted:
                        found[acc] = current_name
                        if len(found) == len(wanted):
                            break
    missing = sorted(wanted - set(found))
    if missing:
        print(f"[domain_recovery] WARN: PF accession(s) not found in Pfam-A: {','.join(missing)}")
    return found


def pfam_tokens(raw: str) -> list[str]:
    """Return normalized PFxxxxx accessions from a config field."""
    out: set[str] = set()
    for token in re.split(r"[+|,;\s]+", raw or ""):
        token = token.strip()
        if re.fullmatch(r"PF\d{5}(?:\.\d+)?", token):
            out.add(token.split(".")[0])
    return sorted(out)


def pfam_accessions_from_anchors(anchors: list[dict]) -> list[str]:
    """Collect PFxxxxx accessions declared in anchor `pfam` fields."""
    out: set[str] = set()
    for a in anchors:
        out.update(pfam_tokens(a.get("pfam", "")))
    return sorted(out)


def _as_float(raw: str, default: float = float("inf")) -> float:
    try:
        return float(raw)
    except (TypeError, ValueError):
        return default


def parse_pfam_domtblout(
    path: Path,
    max_seq_evalue: Optional[float] = None,
    max_domain_cevalue: Optional[float] = None,
    max_domain_ievalue: Optional[float] = None,
) -> dict[str, set[str]]:
    """Parse HMMER domtblout into locus_tag -> set(domain names).

    Supports both:
      - hmmscan: target=domain, query=protein
      - hmmsearch: target=protein, query=domain
    """
    domains: dict[str, set[str]] = {}
    if not path or not path.exists():
        return domains
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(None, 22)
            if len(parts) < 4:
                continue
            seq_evalue = _as_float(parts[6]) if len(parts) > 6 else float("inf")
            domain_cevalue = _as_float(parts[11]) if len(parts) > 11 else float("inf")
            domain_ievalue = _as_float(parts[12]) if len(parts) > 12 else float("inf")
            passes = True
            if max_seq_evalue is not None:
                passes = seq_evalue <= max_seq_evalue
            if max_domain_cevalue is not None:
                passes = passes or domain_cevalue <= max_domain_cevalue
            if max_domain_ievalue is not None:
                passes = passes or domain_ievalue <= max_domain_ievalue
            if not passes:
                continue
            # HMMER domtblout columns: 0 target name, 1 target acc, 3 query name, 4 query acc.
            if len(parts) > 4 and parts[4].startswith("PF"):
                target = parts[0]
                domain = parts[3]
                domain_acc = parts[4].split(".")[0]
            else:
                domain = parts[0]
                target = parts[3]
                domain_acc = parts[1].split(".")[0] if len(parts) > 1 and parts[1].startswith("PF") else ""
            domains.setdefault(target, set()).add(domain)
            if domain_acc:
                domains[target].add(domain_acc)
    return domains


def domain_architecture(domains: set[str], required_pfams: set[str], label: str) -> str:
    """Classify a hit against a configured PFAM architecture."""
    present = required_pfams & domains
    if required_pfams and present == required_pfams:
        return label
    if present:
        return "+".join(sorted(present)) + "_only"
    return ""


def same_strand_neighbors(features: list[dict], target_lt: str) -> tuple[Optional[dict], Optional[dict]]:
    ordered = sorted(features, key=lambda g: (g.get("contig", ""), int(g.get("start", 0))))
    idx = next((i for i, g in enumerate(ordered) if g.get("locus_tag") == target_lt), None)
    if idx is None:
        return None, None
    contig = ordered[idx].get("contig")
    strand = ordered[idx].get("strand")
    prev_g = None
    for i in range(idx - 1, -1, -1):
        g = ordered[i]
        if g.get("contig") != contig:
            break
        if g.get("strand") == strand:
            prev_g = g
            break
    next_g = None
    for i in range(idx + 1, len(ordered)):
        g = ordered[i]
        if g.get("contig") != contig:
            break
        if g.get("strand") == strand:
            next_g = g
            break
    return prev_g, next_g


def split_domain_class(
    target_domains: set[str],
    required_pfams: set[str],
    required_pfams_ordered: list[str],
    label: str,
    prev_g: Optional[dict],
    next_g: Optional[dict],
    domains_by_lt: dict[str, set[str]],
) -> Optional[tuple[str, str, str]]:
    """Return (full_class, partner_locus_tag, kind) for split architecture calls."""
    target_present = required_pfams & target_domains
    if not target_present or target_present == required_pfams:
        return None
    missing = required_pfams - target_present
    for nb in (prev_g, next_g):
        if not nb:
            continue
        nb_lt = nb.get("locus_tag", "")
        nb_present = required_pfams & domains_by_lt.get(nb_lt, set())
        if missing and missing <= nb_present:
            return (label, nb_lt, "split_cds")
        terminal_domain = required_pfams_ordered[-1] if required_pfams_ordered else ""
        if nb.get("is_pseudo") and terminal_domain in target_present:
            return (label, nb_lt, "split_pseudo")
    return None


def write_domain_recovery_diagnostics(rows: list[dict], out_path: Path) -> None:
    if not rows:
        return
    fields = [
        "genome_acc", "locus_tag", "family", "subtype", "state", "final_class",
        "pfam_arch", "domains", "partner_locus_tag", "contig", "start", "end",
        "reason",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})
