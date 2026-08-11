#!/usr/bin/env python3
"""Interactive locus setup and run."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

_PKG_DIR = Path(__file__).resolve().parent
if str(_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(_PKG_DIR))

from build_config import check_anchors_csv
from config_utils import get_target_cfg, load_settings, missing_target_files
from domain_recovery import pfam_tokens
from discover import (
    CONTEXT_GENES,
    find_genomes,
    paint,
    search_any_genome,
    format_candidates,
    format_neighborhood,
    genome_label,
    group,
    neighborhood,
    search,
)


# Returned by any prompt when the answer is "take me back one step".
BACK = object()
_BACK_WORDS = ("b", "back")

# Answers that mean "show me more genes", and how many each adds
# (before, after). The gene the user wants is often just outside the
# window the search hit produced.
_WIDEN_STEP = 10
_WIDEN = {"<": (_WIDEN_STEP, 0),
          ">": (0, _WIDEN_STEP),
          "+": (_WIDEN_STEP, _WIDEN_STEP)}


def _ask(prompt: str, valid: range, default: int | None = None,
         widen: bool = False):
    """Read a number, re-asking until it is one of the offered ones.

    With `widen`, also accept <, > and + and hand them back as sentinels,
    so a caller listing genes can show more of them rather than making the
    user go back and search again.
    """
    while True:
        try:
            raw = input(prompt).strip()
        except EOFError:
            sys.exit("\nAborted (no input available).")
        if raw.lower() in _BACK_WORDS:
            return BACK
        if widen and raw in _WIDEN:
            return raw
        if not raw:
            if default is not None:
                return default
            sys.exit("Aborted.")
        if raw.isdigit() and int(raw) in valid:
            return int(raw)
        print(f"  Enter a number between {valid.start} and {valid.stop - 1}, "
              + ("<, > or + to list more genes, " if widen else "")
              + "or b to go back.")


def _confirm(prompt: str, default: bool = False) -> bool:
    try:
        answer = input(prompt).strip().lower()
    except EOFError:
        return False
    if not answer:
        return default
    return answer in ("y", "yes")


def _editor() -> list[str] | None:
    """The user's editor, found the way git and crontab find it."""
    for variable in ("VISUAL", "EDITOR"):
        value = os.environ.get(variable)
        if value:
            return shlex.split(value)
    for name in ("nano", "vim", "vi"):
        found = shutil.which(name)
        if found:
            return [found]
    return None


def _prompt(prompt: str, allow_back: bool = False):
    """Read a line."""
    try:
        raw = input(prompt).strip()
    except EOFError:
        sys.exit("\nAborted (no input available).")
    if allow_back and raw.lower() in _BACK_WORDS:
        return BACK
    return raw


def _runner(script: str = "locus-ruler", module: str = "locus_ruler.run") -> list[str]:
    """How to invoke a sibling command from here."""
    sibling = Path(sys.executable).parent / script
    if sibling.exists():
        return [str(sibling)]
    return [sys.executable, "-m", module]


# The four NCBI levels, most finished first.
_LEVELS = [
    ("complete", "one closed sequence per replicon"),
    ("chromosome", "essentially whole, may still have gaps"),
    ("scaffold", "gaps, but the order of the pieces is known"),
    ("contig", "unordered draft assembly"),
]


def _parse_levels(raw: str, names: list[str]) -> list[str]:
    """Read '1 3', '1,3' or 'complete scaffold' as a set of levels."""
    picked = set()
    for token in raw.replace(",", " ").split():
        if token.isdigit() and 1 <= int(token) <= len(names):
            picked.add(names[int(token) - 1])
            continue
        matches = [n for n in names if n.startswith(token.lower())]
        if len(matches) != 1:
            return []
        picked.add(matches[0])
    return [n for n in names if n in picked]


def _pick_levels() -> list[str]:
    """Which assembly levels to download."""
    names = [name for name, _ in _LEVELS]
    print(paint("\nWhich assembly levels do you want?", "bold"))
    print("A locus split across two contigs can read as a deletion.\n")
    for index, (name, why) in enumerate(_LEVELS, start=1):
        print(f"  {paint(f'[{index}]', 'bold', 'cyan')} {name:<12s}"
              + paint(why, "dim"))
    print("\nPick as many as you like: " + paint("1 2", "bold")
          + ", or the names, or Enter for all.")
    while True:
        raw = _prompt("  levels [all]: ")
        if not raw:
            return names
        picked = _parse_levels(raw, names)
        if picked:
            return picked
        print(f"  Use numbers 1-{len(names)} or the level names.")


def _pick_contig_cap(default: int = 100) -> int:
    """The contig-count ceiling, asked without using the word seq_count."""
    print(paint("\nHow fragmented an assembly will you accept?", "bold"))
    print("This is the number of contigs a genome is broken into: a finished")
    print("one is a handful, a rough draft can be hundreds. 0 = no limit.")
    while True:
        raw = _prompt(f"  most contigs per genome [{default}]: ")
        if not raw:
            return default
        if raw.isdigit():
            return int(raw)
        print("  Enter a whole number, or press Enter for "
              f"{default}.")


def _dataset_dir(name: str) -> Path:
    """Where a downloaded dataset's db/gff/genomes/faa land: input/<name>/."""
    return Path("input") / name


def _build_dataset(settings_path: Path) -> str | None:
    """Offer to download a dataset, rather than stopping at a command to copy."""
    print("LocusRuler reads genomes from a local database, which it can")
    print("download from NCBI for you.")
    print("\nWhich organisms? A genus, a species, or several separated by spaces.")
    print("  e.g. Escherichia")
    print("       Escherichia Salmonella Klebsiella")
    taxa = _prompt("  taxon: ").split()
    if not taxa:
        return None

    default_name = taxa[0].lower() if len(taxa) == 1 else "dataset"
    name = _prompt(f"  name for this dataset [{default_name}]: ") or default_name
    outdir = _dataset_dir(name)

    setup = _runner("locus-ruler-setup", "locus_ruler.setup")
    where = ["--taxon", *taxa, "--outdir", str(outdir),
             "--assembly-level", *_pick_levels(),
             "--seq-count", str(_pick_contig_cap())]

    print(paint("\nChecking how many genomes that is...", "dim"))
    if subprocess.call(setup + where + ["--dry-run"]) != 0:
        print(paint("  Survey failed; check the taxon spelling.", "yellow"))
        return None

    build = setup + where + [
        "--db-name", name, "--add-target", str(settings_path)]
    print(paint("\nDownloading can take a long time", "yellow")
          + " and needs disk space for every genome.")
    if not _confirm("Download now? [y/N] "):
        print("\nRun this when you are ready, then start the wizard again:")
        print("  " + paint(" ".join(build), "bold"))
        return None

    if subprocess.call(build) != 0:
        print(paint("\nDataset build failed.", "yellow"))
        return None
    print(paint(f"\nDataset '{name}' is ready.", "green"))
    return name


def _pick_settings(given: str | None) -> Path:
    """Find settings.toml, or offer to create it from the example."""
    if given:
        path = Path(given)
        if not path.exists():
            sys.exit(f"ERROR: settings file not found: {path}")
        return path
    for candidate in (Path("settings.toml"), Path("locus_ruler/settings.toml")):
        if candidate.exists():
            print(f"Using {candidate}")
            return candidate

    example = Path("settings.example.toml")
    if example.exists():
        print(paint("\nNo settings.toml yet.", "bold")
              + f" This is usually a one-time copy of {example.name}.")
        if _confirm(f"  Create settings.toml from it now? [Y/n] ", default=True):
            dest = Path("settings.toml")
            dest.write_text(example.read_text(encoding="utf-8"),
                            encoding="utf-8")
            print(paint(f"  Wrote {dest}.", "green"))
            return dest

    sys.exit("ERROR: no settings.toml here. Copy the example first:\n"
             "         cp settings.example.toml settings.toml\n"
             "       or pass --settings <path>.")


def _pick_target(settings: dict, given: str | None, settings_path: Path):
    """Offer datasets whose files actually exist, not just ones declared."""
    declared = settings.get("targets", [])
    ready = [t["name"] for t in declared if not missing_target_files(t)]
    unready = [t["name"] for t in declared if t["name"] not in ready]
    if given:
        return given

    if not ready:
        if unready:
            print(paint("\nA dataset is declared in settings.toml but its "
                        "files aren't downloaded yet: ", "bold")
                  + ", ".join(unready))
        else:
            print(paint("\nNo dataset yet.", "bold"))
        built = _build_dataset(settings_path)
        if not built:
            sys.exit("\nNothing to analyze yet. Start the wizard again once "
                     "you have a dataset.")
        return built

    # Building a dataset is on this menu even when one already exists.
    print(paint("\nWhich dataset?", "bold") + "\n")
    for index, name in enumerate(ready, start=1):
        print(f"  {paint(f'[{index}]', 'bold', 'cyan')} {name}")
    fresh = len(ready) + 1
    print(f"  {paint(f'[{fresh}]', 'bold', 'cyan')} "
          + paint("download a new one from NCBI", "dim"))

    choice = _ask(f"\nWhich one? [1-{fresh}, Enter = 1] ",
                  range(1, fresh + 1), default=1)
    if choice is BACK:
        return BACK
    if choice < fresh:
        return ready[choice - 1]
    built = _build_dataset(settings_path)
    if not built:
        sys.exit("\nNo new dataset. Start the wizard again when you have one.")
    return built


def _pick_accession(db: Path, given: str | None):
    """Choose the reference genome the same way as everything else: by number."""
    if given:
        return given
    print(paint("\nWhich genome should be the reference?", "bold"))
    print("The cluster is defined in this genome, and every other genome is")
    print("compared against it.")
    print("Search by strain name, species, or accession.")
    while True:
        term = _prompt("  search (b = back): ", allow_back=True)
        if term is BACK:
            return BACK
        if not term:
            sys.exit("Aborted.")
        matches = find_genomes(db, term)
        if not matches:
            print(f"  Nothing matches {term!r}. Try a shorter word.")
            continue
        if len(matches) > 25:
            print(f"  {len(matches)} genomes match; narrow it down.")
            continue
        print()
        for index, g in enumerate(matches, start=1):
            print(f"  {paint(f'[{index}]', 'bold', 'cyan')} "
                  f"{g['accession']:18s} {g['species']} {g['strain']}"
                  + paint(f"  [{g['assembly_level']}]", "dim"))
        choice = _ask(f"\nWhich one? [1-{len(matches)}, b = back] ",
                      range(1, len(matches) + 1))
        if choice is BACK:
            continue
        return matches[choice - 1]["accession"]


def _db_for(ctx: dict, target: str) -> Path:
    return Path(get_target_cfg(ctx["settings"], target)["db"])


def _given(ctx: dict, key: str):
    """A value from the command line, usable once."""
    return ctx["given"].pop(key, None)


def _step_target(state: dict, ctx: dict):
    target = _pick_target(ctx["settings"], _given(ctx, "target"),
                          ctx["settings_path"])
    if target is BACK:
        return BACK
    if target not in [t["name"] for t in ctx["settings"].get("targets", [])]:
        # _build_dataset appended a new [[targets]] block; reread it
        ctx["settings"] = load_settings(ctx["settings_path"])
    return target


def _step_accession(state: dict, ctx: dict):
    return _pick_accession(_db_for(ctx, state["target"]),
                           _given(ctx, "accession"))


def _step_term(state: dict, ctx: dict):
    db = _db_for(ctx, state["target"])
    accession = state["accession"]
    print("\n" + paint("Reference: ", "bold")
          + paint(accession, "cyan") + f"  {genome_label(db, accession)}")
    term = _given(ctx, "term")
    while True:
        if not term:
            print(paint("\nWhat is the cluster about?", "bold"))
            print("A word from a gene product (transporter, secretion, "
                  "transposase), or a locus_tag.")
            term = _prompt("  search (b = back): ", allow_back=True)
            if term is BACK:
                return BACK
            if not term:
                sys.exit("Aborted.")
        if search(db, accession, term):
            return term
        print(f"  Nothing in {accession} matches {term!r}. Try another word.")
        term = None


def _step_candidate(state: dict, ctx: dict):
    db = _db_for(ctx, state["target"])
    hits = search(db, state["accession"], state["term"])
    groups = group(hits)
    found = "gene matches" if len(hits) == 1 else "genes match"
    where = "neighborhood" if len(groups) == 1 else "neighborhoods"
    print(f"\n{len(hits)} {found} {state['term']!r}, "
          f"in {len(groups)} {where}:\n")
    for line in format_candidates(groups):
        print(line)
    choice = _ask(f"\nWhich one? [1-{len(groups)}, b = back] ",
                  range(1, len(groups) + 1))
    if choice is BACK:
        return BACK
    return groups[choice - 1]


def _step_bounds(state: dict, ctx: dict):
    db = _db_for(ctx, state["target"])
    picked = state["candidate"]
    pad_before = pad_after = CONTEXT_GENES

    def _show():
        genes = neighborhood(db, state["accession"], picked["contig"],
                              picked["start"], picked["end"],
                              pad_before=pad_before, pad_after=pad_after)
        # The legend sits at the prompt, where it is still on screen.
        print(f"\n{picked['contig']}, genes in order:\n")
        for line in format_neighborhood(genes):
            print(line)
        return genes

    genes = _show()

    # Ask for the locus itself, not for what surrounds it.
    while True:
        print(paint("\nWhich genes are the two ends of the cluster?", "bold"))
        print(f"(order does not matter; * = matched {state['term']!r}; "
              "b = back)")
        print(paint("(gene not listed? < shows more above, > more below, "
                    "+ both)", "dim"))
        a = _ask(f"  From [1-{len(genes)}] ", range(1, len(genes) + 1),
                 widen=True)
        if a is BACK:
            return BACK
        if a in _WIDEN:
            before, after = _WIDEN[a]
            pad_before += before
            pad_after += after
            genes = _show()
            continue
        b = _ask(f"  To   [1-{len(genes)}] ", range(1, len(genes) + 1))
        if b is BACK:
            continue
        first, last = sorted((a, b))
        if first == last:
            print(paint("\n  A locus needs at least two genes; "
                        "pick two different ends.", "yellow"))
            continue
        if first == 1 or last == len(genes):
            print(paint("\n  That end has no neighbor in view to use as a "
                        "flank.", "yellow"))
            print("  Pick ends further inside the listing.")
            continue

        flank_l, flank_r = genes[first - 2], genes[last]
        inside = genes[first - 1:last]
        width = 58
        print("\n  " + paint(f"{'flank':>7}  {flank_l['locus_tag']:<18s} "
                             f"{flank_l['product'][:44]}", "dim"))
        print(f"  {'-' * width}")
        print("  " + paint(f"{len(inside)} genes, "
                           f"{flank_r['start'] - flank_l['end']:,} bp",
                           "bold", "green"))
        for index, gene in enumerate(inside, start=1):
            edge = ("first" if index == 1
                    else ("last" if index == len(inside) else ""))
            tag = paint(f"{gene['locus_tag']:<18s}", "cyan")
            print(f"  {paint(f'{edge:>7}', 'green')}  {tag} "
                  f"{gene['product'][:44]}")
        print(f"  {'-' * width}")
        print("  " + paint(f"{'flank':>7}  {flank_r['locus_tag']:<18s} "
                           f"{flank_r['product'][:44]}", "dim"))

        if _confirm("\nUse this as the locus? [y/N] "):
            return {"flank_L": flank_l["locus_tag"],
                    "flank_R": flank_r["locus_tag"]}
        # Declining means the ends were wrong, so ask again rather than quit.


def _default_outdir(ctx: dict, target: str, locus_id: str) -> Path:
    """Where the config (and everything step 1 writes beside it) goes
    when `--outdir` is not given: next to where the results will land.

    Without this, the config, the anchor table, anchors.faa and the cluster
    FASTA default to the current directory -- wherever the wizard happened
    to be started from -- while run.py writes tables/, diagnostics/ and the
    workbook to `output_root/<target>/<locus_id>`. Two trees for one locus,
    connected only by memory of which shell you were in. `output_root` is
    read from the same settings run.py itself resolves the path from
    (config_utils.load_settings makes it absolute at load time), so this
    stays correct wherever the wizard is launched from.
    """
    output_root = Path(ctx["settings"]["paths"]["output_root"])
    return output_root / target / locus_id


def _write_config(state: dict, ctx: dict) -> Path:
    args = ctx["args"]
    default_locus_id = f"{state['term']}_{state['accession'].replace('.', '_')}"
    if args.locus_id or not sys.stdin.isatty():
        locus_id = args.locus_id or default_locus_id
    else:
        # Left to itself this defaults to term_accession, which is unique but
        # long -- a real barrier to re-running (--config <that>.json) by hand.
        locus_id = _prompt(
            f"\n  name for this locus [{default_locus_id}]: ") or default_locus_id
    outdir = (Path(args.outdir) if args.outdir
             else _default_outdir(ctx, state["target"], locus_id))
    out = outdir / f"{locus_id}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "_description": "Written by locus-ruler-wizard.",
        "locus_id": locus_id,
        "reference": {
            "accession": state["accession"],
            "target": state["target"],
            **state["bounds"],
        },
    }, indent=2) + "\n", encoding="utf-8")
    return out


# What can happen once the config exists.
_ACTIONS = [
    ("run the whole pipeline", "BLAST every genome; minutes to hours"),
    ("set the gene roles here", "one prompt per gene, then run"),
    ("stop at the anchor table", "edit it yourself, then resume"),
    ("print the command only", "nothing is run"),
    ("go back", "change the ends, the reference, the dataset"),
]


# The roles worth asking about for a cluster gene, in the order offered.
_ROLE_HELP = [
    ("CORE", "defines the locus. ABSENT, DECAYED and DIVERGENT are decided",
             "on these genes alone"),
    ("SHARED", "used here but also elsewhere in the genome. Missing, it turns",
               "an otherwise complete locus CONDITIONAL"),
    ("ASSOCIATED", "belongs to the locus but is not required for it to work,",
                   "such as an uptake receptor beside a biosynthesis core"),
    ("IGNORE", "not scored. Use it for a gene that sits inside the locus but",
               "should not count toward it (flanks are CONTEXT automatically)"),
]
# Only the genes between the ends.
_CLUSTER_ROLES = ("cluster_L", "inner", "cluster_R")

# How many aux search hits fit on screen before the list stops being readable.
_AUX_SHOWN = 30


# Printed once above the list, then echoed in short form beside every prompt.
_ROLE_HINT = "l/x set a flag and stay on this gene; a role number moves on"


def _print_role_legend() -> None:
    print(paint("\nWhat part does each gene play?", "bold"))
    for index, (name, first, second) in enumerate(_ROLE_HELP, start=1):
        print(f"  {paint(str(index), 'bold', 'cyan')} {name:<11s} {first}")
        if second:
            print(f"    {'':<11s} {second}")
    print("\n  " + paint("Enter", "bold") + " keeps what is shown  "
          + paint("l", "bold") + " toggle lenient  "
          + paint("x", "bold") + " toggle exception  "
          + paint("0", "bold") + " clear both")
    print("  " + paint("?", "bold") + " this list  "
          + paint("p", "bold") + " previous gene  "
          + paint("b", "bold") + " back to the menu")
    print("  " + paint("lenient", "dim")
          + paint(" = accept a weaker match for this gene; ", "dim")
          + paint("exception", "dim")
          + paint(" = annotated but inert", "dim"))
    print("  " + paint(_ROLE_HINT, "dim"))


def _bracket_value(row: dict, current: str) -> str:
    """The `[current]` a prompt shows, with any active flags folded in."""
    marks = [name for name in ("lenient", "exception")
             if (row.get(name) or "").strip().upper() == "TRUE"]
    if not marks:
        return paint(current, "bold", "green")
    return paint(current, "bold", "green") + paint(f", {' '.join(marks)}", "yellow")


def _save_anchors(anchors: Path, fields: list[str], rows: list[dict]) -> None:
    with anchors.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)
    print(paint(f"\n  Saved {anchors.name}.", "green"))


# Previous-gene words, kept distinct from `_BACK_WORDS`, which leaves the step entirely.
_PREV_WORDS = ("p", "prev", "previous")


def _curate_roles(anchors: Path) -> bool:
    """Ask for each cluster gene's role, instead of handing over a CSV."""
    with anchors.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    editable = [i for i, row in enumerate(rows)
                if (row.get("role") or "") in _CLUSTER_ROLES]
    if not editable:
        print(paint("  No cluster genes in the anchor table.", "yellow"))
        return False

    _print_role_legend()
    names = [name for name, _, _ in _ROLE_HELP]
    position = 0
    while position < len(editable):
        row = rows[editable[position]]
        current = (row.get("status_role") or "").strip().upper() or "CORE"
        print(f"\n  {paint(f'{position + 1}/{len(editable)}', 'dim')} "
              + paint(f"{row.get('locus_tag', ''):<18s}", "cyan")
              + f"{(row.get('product') or '')[:46]}")
        print("       " + paint(_ROLE_HINT, "dim"))
        answer = _prompt(f"       role [{_bracket_value(row, current)}]: "
                         ).strip().lower()

        if answer in _BACK_WORDS:
            _save_anchors(anchors, fields, rows)
            return True
        if answer in _PREV_WORDS:
            position = max(position - 1, 0)
            continue
        if answer == "?":
            _print_role_legend()
            continue
        if answer == "0":
            row["lenient"] = "FALSE"
            row["exception"] = "FALSE"
            continue
        if answer in ("l", "x"):
            column = "lenient" if answer == "l" else "exception"
            on = (row.get(column) or "").strip().upper() == "TRUE"
            row[column] = "FALSE" if on else "TRUE"
            continue          # same gene, so the change is visible before moving on
        if not answer:
            position += 1
            continue
        if answer.isdigit() and 1 <= int(answer) <= len(names):
            row["status_role"] = names[int(answer) - 1]
            position += 1
            continue
        match = [name for name in names if name.lower().startswith(answer)]
        if len(match) == 1:
            row["status_role"] = match[0]
            position += 1
            continue
        print(f"       Enter a number 1-{len(names)}, a role name, "
              "l, x, 0, ?, p or b.")

    _save_anchors(anchors, fields, rows)
    return True


def _curate_families(anchors: Path) -> None:
    """Ask for each cluster gene's family: the paralog-count cap."""
    with anchors.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    editable = [i for i, row in enumerate(rows)
                if (row.get("role") or "") in _CLUSTER_ROLES]
    if not editable:
        print(paint("  No cluster genes in the anchor table.", "yellow"))
        return

    print(paint("\nFamily labels cap paralog counts:", "bold")
          + " two genes with the same")
    print("label are counted as one system. Enter keeps the automatic guess.")
    print(paint("(p previous gene, b back to the menu)\n", "dim"))
    position = 0
    while position < len(editable):
        row = rows[editable[position]]
        current = row.get("family") or "(none)"
        print(f"  {paint(f'{position + 1}/{len(editable)}', 'dim')} "
              + paint(f"{row.get('locus_tag', ''):<18s}", "cyan")
              + f"{(row.get('product') or '')[:46]}")
        answer = _prompt(f"       family [{paint(current, 'bold', 'green')}]: ")
        stripped = answer.strip()
        if stripped.lower() in _BACK_WORDS:
            _save_anchors(anchors, fields, rows)
            return
        if stripped.lower() in _PREV_WORDS:
            position = max(position - 1, 0)
            continue
        if stripped:
            row["family"] = stripped
        position += 1

    _save_anchors(anchors, fields, rows)


def _pick_aux_gene(db: Path) -> dict | None:
    """Search until something is chosen, rather than out at the first miss."""
    while True:
        term = _prompt("\n  search (Enter to cancel): ").strip()
        if not term:
            return None
        hits = search_any_genome(db, term, limit=_AUX_SHOWN)
        if not hits:
            print(f"  Nothing matches {term!r}.")
            continue
        shown = hits[:_AUX_SHOWN]
        print()
        for index, h in enumerate(shown, start=1):
            where = f"{h['species']} {h['strain']}".strip() or h["genome_acc"]
            tag = paint(f"{h['locus_tag']:<18s}", "cyan")
            number = paint(f"[{index}]".rjust(4), "bold", "cyan")
            print(f"  {number} {tag} "
                  f"{where[:28]:<28s} {h['product'][:40]}")
        if len(hits) > len(shown):
            print(paint(f"  ... and more beyond the first {_AUX_SHOWN}. "
                        "A longer search word shortens this.", "dim"))
        choice = _ask(f"\nWhich one? [1-{len(shown)}, b = search again] ",
                      range(1, len(shown) + 1))
        if choice is BACK:
            continue
        return shown[choice - 1]


def _curate_add_aux(anchors: Path, db: Path) -> None:
    """Add a gene left out of the cluster that may still sit beside it."""
    print(paint("\nAn aux anchor is a gene outside the cluster you defined,",
                "bold"))
    print("which may still turn up next to it in some genomes. LocusRuler")
    print("then looks for it and reports where it appears, without letting it")
    print("move the cluster boundaries.")
    print(paint("  Search matches a gene product, a gene name or a locus_tag,",
                "dim"))
    print(paint("  across every genome in this dataset.", "dim"))
    picked = _pick_aux_gene(db)
    if picked is None:
        return

    names = [name for name, _, _ in _ROLE_HELP]
    print("\nRole for this gene (blank = diagnostic only, not scored):")
    for index, (name, first, _) in enumerate(_ROLE_HELP, start=1):
        print(f"  {paint(str(index), 'bold', 'cyan')} {name:<11s} {first}")
    role_answer = _prompt(f"  role [none]: ").strip().lower()
    status_role = ""
    if role_answer.isdigit() and 1 <= int(role_answer) <= len(names):
        status_role = names[int(role_answer) - 1]
    else:
        match = [name for name in names if name.lower().startswith(role_answer)]
        if len(match) == 1:
            status_role = match[0]

    with anchors.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    rows.append({
        "role": "aux", "locus_tag": picked["locus_tag"], "family": "",
        "status_role": status_role, "exception": "FALSE", "lenient": "FALSE",
        "pfam": "", "pfam_split": "FALSE",
    })
    _save_anchors(anchors, fields, rows)
    print(f"  Added {picked['locus_tag']} from {picked['genome_acc']} "
         f"as aux ({status_role or 'diagnostic only'}).")


def _print_pfam_legend() -> None:
    print(paint("\nWhich genes are too diverged for BLAST to find?", "bold"))
    print("  Name the Pfam domains such a gene must carry and LocusRuler")
    print("  looks for them with HMMER once BLAST has come up empty.")
    print("  Write one accession, or several joined by + where the protein")
    print("  carries all of them. " + paint("PF00005", "bold")
          + " and " + paint("PF00005+PF00664", "bold") + " are both valid.")
    print("\n  " + paint("Enter", "bold") + " keeps what is shown  "
          + paint("s", "bold") + " toggle split  "
          + paint("0", "bold") + " clear this gene")
    print("  " + paint("?", "bold") + " this list  "
          + paint("p", "bold") + " previous gene  "
          + paint("b", "bold") + " back to the menu")
    print("  " + paint("split", "dim")
          + paint(" = also accept a protein carrying only some of the domains,",
                  "dim"))
    print("  " + paint("        when the neighbor beside it carries the rest",
                       "dim"))


def _curate_pfam(anchors: Path) -> None:
    """Ask which genes get a Pfam profile to fall back on."""
    with anchors.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    editable = [i for i, row in enumerate(rows)
                if (row.get("role") or "") in _CLUSTER_ROLES]
    if not editable:
        print(paint("  No cluster genes in the anchor table.", "yellow"))
        return

    _print_pfam_legend()
    position = 0
    while position < len(editable):
        row = rows[editable[position]]
        current = row.get("pfam") or "(none)"
        split = (row.get("pfam_split") or "").strip().upper() == "TRUE"
        shown = paint(current, "bold", "green")
        if split:
            shown += paint(", split", "yellow")
        print(f"\n  {paint(f'{position + 1}/{len(editable)}', 'dim')} "
              + paint(f"{row.get('locus_tag', ''):<18s}", "cyan")
              + f"{(row.get('product') or '')[:46]}")
        answer = _prompt(f"       pfam [{shown}]: ").strip()

        if answer.lower() in _BACK_WORDS:
            break
        if answer.lower() in _PREV_WORDS:
            position = max(position - 1, 0)
            continue
        if answer == "?":
            _print_pfam_legend()
            continue
        if answer.lower() == "s":
            row["pfam_split"] = "FALSE" if split else "TRUE"
            continue          # same gene, so the change is visible before moving on
        if answer == "0":
            row["pfam"] = ""
            row["pfam_split"] = "FALSE"
            position += 1
            continue
        if not answer:
            position += 1
            continue
        accessions = pfam_tokens(answer)
        if not accessions:
            print(paint("       Not a Pfam accession. They read PF followed "
                        "by five digits.", "yellow"))
            continue
        row["pfam"] = "+".join(accessions)
        position += 1

    _save_anchors(anchors, fields, rows)


def _build_anchors(config: Path, base: list[str]) -> tuple[int, Path | None]:
    """Run step 1 only, so the anchor table exists and nothing has BLASTed."""
    code = subprocess.call(base + ["--to", "1"])
    if code != 0:
        return code, None
    anchors = config.with_name(config.stem + "_anchors.csv")
    if not anchors.exists():
        print(paint(f"\nNo anchor table at {anchors}.", "yellow"))
        return 1, None
    return 0, anchors


def _resume_hint(resume: list[str]) -> None:
    print("\nWhen the anchor table looks right, continue with:")
    print("  " + paint(" ".join(resume), "bold"))


def _stop_at_anchors(config: Path, base: list[str]) -> int:
    """Write the anchor table and stop, for editing outside the wizard."""
    code, anchors = _build_anchors(config, base)
    if anchors is None:
        return code
    # run.py already printed the file and the resume command; this adds which columns
    # are the user's to change.
    print("\nYours to edit: " + paint("status_role", "cyan") + ", "
          + paint("family", "cyan") + ", " + paint("exception", "cyan")
          + ", " + paint("lenient", "cyan") + ", the HMM")
    print("columns, and aux rows. The rest comes from the genome.")
    editor = _editor()
    if editor:
        print("\n  " + paint(" ".join(editor + [str(anchors)]), "bold"))
    return 0


# A menu, not a forced sequence: each choice returns here, so family, aux, and HMM can
# be set in any combination before running.
_CURATE_OPTIONS = [
    ("gene roles", "CORE / SHARED / ASSOCIATED / IGNORE, plus lenient / exception"),
    ("family labels", "group genes into paralog-count families"),
    ("add an anchor gene (aux)", "pull in an extra gene from any genome"),
    ("Pfam rescue", "find a gene too diverged for BLAST, by its domains"),
    ("run the pipeline now", "BLAST every genome in the target"),
    ("stop here", "edit the file yourself, resume later"),
]


def _curate_then_run(config: Path, base: list[str], db: Path, settings: dict) -> int:
    """Set as much or as little as wanted, then run -- looping, not linear."""
    code, anchors = _build_anchors(config, base)
    if anchors is None:
        return code
    resume = base + ["--from", "2"]
    run_choice = len(_CURATE_OPTIONS) - 1

    while True:
        print(paint("\nAnchor table ready.", "bold") + " What do you want to set?\n")
        for index, (label, why) in enumerate(_CURATE_OPTIONS, start=1):
            print(f"  {paint(f'[{index}]', 'bold', 'cyan')} {label:<26s}"
                  + paint(why, "dim"))
        choice = _ask(f"\nWhich one? [1-{len(_CURATE_OPTIONS)}, "
                      f"Enter = {run_choice}] ",
                      range(1, len(_CURATE_OPTIONS) + 1), default=run_choice)
        if choice is BACK:
            continue

        if choice == 1:
            _curate_roles(anchors)
        elif choice == 2:
            _curate_families(anchors)
        elif choice == 3:
            _curate_add_aux(anchors, db)
        elif choice == 4:
            _curate_pfam(anchors)
        elif choice == run_choice:
            problems = check_anchors_csv(anchors, settings)
            if problems:
                print(paint("\nThat would not do what you meant:", "yellow"))
                for problem in problems[:8]:
                    print(f"  - {problem}")
                if len(problems) > 8:
                    print(paint(f"  ... and {len(problems) - 8} more", "dim"))
                continue
            return subprocess.call(resume)
        else:
            _resume_hint(resume)
            return 0


def _step_action(state: dict, ctx: dict):
    out = _write_config(state, ctx)
    print("\n" + paint("Wrote ", "green") + str(out))

    base = _runner() + ["--config", str(out),
                        "--settings", str(ctx["settings_path"]),
                        "--target", state["target"]]
    print(paint("\nWhat now?", "bold") + "\n")
    for index, (label, why) in enumerate(_ACTIONS, start=1):
        print(f"  {paint(f'[{index}]', 'bold', 'cyan')} {label:<28s}"
              + paint(why, "dim"))
    choice = _ask(f"\nWhich one? [1-{len(_ACTIONS)}, Enter = 1] ",
                  range(1, len(_ACTIONS) + 1), default=1)
    if choice is BACK or choice == len(_ACTIONS):
        return BACK

    if choice == 1:
        print(paint("\nRunning. This BLASTs every genome in the target.",
                    "yellow"))
        return subprocess.call(base)
    if choice == 2:
        return _curate_then_run(out, base, _db_for(ctx, state["target"]), ctx["settings"])
    if choice == 3:
        return _stop_at_anchors(out, base)

    print("\nNot run. Copy this when you are ready:")
    print("  " + paint(" ".join(base), "bold"))
    return 0


# Order matters: each step may read earlier decisions, and going back discards
# everything downstream.
_STEPS = [
    ("target", _step_target),
    ("accession", _step_accession),
    ("term", _step_term),
    ("candidate", _step_candidate),
    ("bounds", _step_bounds),
    ("action", _step_action),
]


def _walk(ctx: dict) -> int:
    """Run the steps in order, letting `b` at any prompt undo the last one."""
    state: dict = {}
    index = 0
    names = [name for name, _ in _STEPS]
    while index < len(_STEPS):
        name, handler = _STEPS[index]
        result = handler(state, ctx)
        if result is BACK:
            if index == 0:
                print(paint("  Already at the first question.", "dim"))
                continue
            index -= 1
            for later in names[index:]:
                state.pop(later, None)
            continue
        if name == "action":
            return result
        state[name] = result
        index += 1
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Find a locus by browsing, then set it up and run it.",
        epilog="Everything is optional: whatever you leave out is asked for.")
    ap.add_argument("find", nargs="?", default=None,
                    help="A word from the product name (transporter), or a "
                         "locus_tag. Plain text, not a regular expression")
    ap.add_argument("--find", dest="find_opt", default=None,
                    help=argparse.SUPPRESS)
    ap.add_argument("--settings", default=None,
                    help="Default: settings.toml in the working directory")
    ap.add_argument("--target", default=None,
                    help="Default: the only target, or you are asked")
    ap.add_argument("--accession", default=None,
                    help="Reference genome, the one the cluster is defined in")
    ap.add_argument("--locus-id", default=None,
                    help="Name for this locus. Skips the interactive prompt "
                         "for it; otherwise you are asked, with "
                         "<find>_<accession> offered as the default")
    ap.add_argument("--outdir", default=None,
                    help="Where to write the config and anchor table. "
                         "Default: output_root/<target>/<locus_id>, the same "
                         "place run.py puts the results")
    args = ap.parse_args(argv)

    settings_path = _pick_settings(args.settings)
    return _walk({
        "args": args,
        "settings_path": settings_path,
        "settings": load_settings(settings_path),
        "given": {"target": args.target,
                  "accession": args.accession,
                  "term": args.find or args.find_opt},
    })


if __name__ == "__main__":
    sys.exit(main())
