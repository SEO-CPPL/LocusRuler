#!/usr/bin/env python3
"""Genome database builder."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import sqlite3
import subprocess
import sys
import time
import urllib.request
from collections import Counter
from datetime import date
from pathlib import Path

_PKG_DIR = Path(__file__).resolve().parent
if str(_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(_PKG_DIR))

# FAA building lives in db_utils so other modules can import it too
from db_utils import (
    append_to_faa,
    build_genome_faa,
    build_pid_to_lt_map,
)

TODAY = date.today().isoformat()
_ATTR_RE = re.compile(r'(\w+)=([^;]+)')

_EPILOG = """\
Several taxa can go into one database: --taxon GenusA GenusB

Output layout (<outdir>/)
  <db_name>.db          SQLite genome + protein database
  gff/                  one *.gff3 per genome
  genomes/              *_genomic.fna.gz and *_protein.faa.gz
  combined_proteins.faa all proteins, locus_tag headers
  _download_tmp/        staging, cleaned after each run
"""


# ── Helpers ──────────────────────────────────────────────────────
def _slug(name: str) -> str:
    """'Gram-negative Enterobacteriaceae' → 'gram_negative_enterobacteriaceae'."""
    return re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_')


def _base_id(acc: str) -> str:
    """'GCF_012345678.1' → '012345678'  (strip prefix and version)."""
    return acc.split('_', 1)[1].split('.')[0]


# ── NCBI query ──────────────────────────────────────────────────────
# What NCBI calls each level, keyed by the shorter word people actually say.
ASSEMBLY_LEVELS = {
    'complete':   'complete_genome',
    'chromosome': 'chromosome',
    'scaffold':   'scaffold',
    'contig':     'contig',
}


def query_ncbi(
    taxon: str,
    seq_max: int,
    existing_accs: set[str],
    fa_dir: Path,
    levels: list[str] | None = None,
) -> list[dict]:
    """Return assemblies for *taxon* that are missing from the local DB/files."""
    existing_fna = {f.name for f in fa_dir.glob('*_genomic.fna.gz')}
    existing_faa = {f.name for f in fa_dir.glob('*_protein.faa.gz')}

    wanted = levels or list(ASSEMBLY_LEVELS)
    level_filter = ''.join(
        f'&filters.assembly_level={ASSEMBLY_LEVELS[lv]}' for lv in wanted)

    results: list[dict] = []
    page_token: str | None = None

    while True:
        url = (
            f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxon}'
            '/dataset_report'
            f'?page_size=500{level_filter}'
        )
        if page_token:
            url += f'&page_token={page_token}'
        req = urllib.request.Request(
            url, headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'}
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read())

        for rpt in data.get('reports', []):
            acc = rpt.get('accession', '')
            if not acc:
                continue
            asm   = rpt.get('assembly_info', {})
            stats = rpt.get('assembly_stats', {})
            org   = rpt.get('organism', {})

            if asm.get('assembly_status', '') == 'suppressed':
                continue

            seq_count = int(
                stats.get('number_of_contigs') or
                stats.get('number_of_component_sequences') or
                9999
            )
            if seq_max and seq_count > seq_max:
                continue

            has_fna = any(f.startswith(acc) for f in existing_fna)
            has_faa = any(f.startswith(acc) for f in existing_faa)
            if acc in existing_accs and has_fna and has_faa:
                continue  # fully present, skip

            results.append({
                'accession':      acc,
                'assembly_name':  asm.get('assembly_name', ''),
                'organism_name':  org.get('organism_name', ''),
                'species':        org.get('organism_name', '').split(' var.')[0],
                'strain':         org.get('infraspecific_names', {}).get('strain', ''),
                'assembly_level': asm.get('assembly_level', ''),
                'seq_count':      seq_count,
                'total_length':   int(stats.get('total_sequence_length') or 0),
            })

        page_token = data.get('next_page_token')
        if not page_token:
            break
        time.sleep(0.2)

    return results


# ── Download helpers ──────────────────────────────────────────────────────
def _download_batch_gcf(
    accessions: list[str], tmp_dir: Path, batch_size: int
) -> None:
    """Download a list of GCF accessions via ncbi-genome-download in batches."""
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        n     = i // batch_size + 1
        total = (len(accessions) - 1) // batch_size + 1
        print(f'  Batch {n}/{total}: {len(batch)} assemblies')
        cmd = [
            # sys.executable, since the venv interpreter may be the only one with the package.
            sys.executable, '-m', 'ncbi_genome_download',
            '--section', 'refseq',
            '--formats', 'fasta,gff,protein-fasta',
            '--assembly-accessions', ','.join(batch),
            '--output-folder', str(tmp_dir),
            '--flat-output',
            '--retries', '3',
            '--parallel', '4',
            'bacteria',
        ]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            if 'No module named' in r.stderr:
                sys.exit(
                    'ERROR: ncbi-genome-download is not installed for '
                    f'{sys.executable}.\n'
                    '       Install it with:  python -m pip install ncbi-genome-download'
                )
            print(f'  [WARN] ncbi-genome-download error: {r.stderr[:300]}')
        time.sleep(1)


def _download_gca_ftp(acc: str, tmp_dir: Path) -> bool:
    """Download a GCA assembly directly from NCBI FTP.  Returns True on success."""
    num  = acc.split('_')[1].split('.')[0]
    path = f'{num[:3]}/{num[3:6]}/{num[6:9]}'
    base = f'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{path}/'

    try:
        req = urllib.request.Request(base, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15) as r:
            html = r.read().decode()
    except Exception as exc:
        print(f'  [WARN] {acc}: FTP directory listing failed ({exc})')
        return False

    dirs = re.findall(r'href="(GCA_[^"]+/)"', html)
    target_dir = next((d for d in dirs if acc.split('.')[0] in d), None)
    if not target_dir and dirs:
        target_dir = dirs[0]
    if not target_dir:
        print(f'  [WARN] {acc}: no FTP directory found')
        return False

    stem      = target_dir.rstrip('/')
    inner_url = base + target_dir

    try:
        req2 = urllib.request.Request(inner_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req2, timeout=15) as r:
            inner = r.read().decode()
    except Exception:
        inner = ''

    available = set(re.findall(r'href="([^"]+)"', inner))
    if not any('_genomic.gff' in f for f in available):
        print(f'  [SKIP] {acc}: no GFF annotation (unannotated draft)')
        return False

    files_to_get = [
        (f'{inner_url}{stem}_genomic.gff.gz', f'{stem}_genomic.gff.gz'),
        (f'{inner_url}{stem}_genomic.fna.gz', f'{stem}_genomic.fna.gz'),
    ]
    if any('_protein.faa' in f for f in available):
        files_to_get.append(
            (f'{inner_url}{stem}_protein.faa.gz', f'{stem}_protein.faa.gz')
        )

    ok = True
    for url, fname in files_to_get:
        dst = tmp_dir / fname
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=120) as r:
                dst.write_bytes(r.read())
        except Exception as exc:
            print(f'  [WARN] {acc} {fname}: {exc}')
            ok = False
    return ok


# ── File organization ──────────────────────────────────────────────────────
def _organize_files(
    tmp_dir: Path, fa_dir: Path, gff_dir: Path
) -> dict[str, dict[str, Path]]:
    """Move downloaded files from *tmp_dir* to their final locations."""
    downloaded: dict[str, dict[str, Path]] = {}

    for fna in list(tmp_dir.glob('*_genomic.fna.gz')):
        m = re.search(r'(GC[FA]_\d+\.\d+)', fna.name)
        if not m:
            continue
        acc = m.group(1)
        dst = fa_dir / fna.name
        shutil.move(str(fna), str(dst))
        downloaded.setdefault(acc, {})['fna'] = dst

    for gff_gz in list(tmp_dir.glob('*_genomic.gff.gz')):
        m = re.search(r'(GC[FA]_\d+\.\d+)', gff_gz.name)
        if not m:
            continue
        acc = m.group(1)
        dst = gff_dir / f'{acc}.gff3'
        with gzip.open(gff_gz, 'rt') as fin, open(dst, 'w') as fout:
            shutil.copyfileobj(fin, fout)
        gff_gz.unlink(missing_ok=True)
        downloaded.setdefault(acc, {})['gff'] = dst

    for faa_gz in list(tmp_dir.glob('*_protein.faa.gz')):
        m = re.search(r'(GC[FA]_\d+\.\d+)', faa_gz.name)
        if not m:
            continue
        acc = m.group(1)
        dst = fa_dir / faa_gz.name
        shutil.move(str(faa_gz), str(dst))
        downloaded.setdefault(acc, {})['faa'] = dst

    return downloaded


# ── DB helpers ──────────────────────────────────────────────────────
DB_SCHEMA = """
CREATE TABLE IF NOT EXISTS genomes (
    accession      TEXT PRIMARY KEY,
    assembly_name  TEXT,
    organism_name  TEXT,
    species        TEXT,
    strain         TEXT,
    assembly_level TEXT,
    seq_count      INTEGER,
    total_length   INTEGER,
    protein_count  INTEGER,
    download_date  TEXT
);
CREATE TABLE IF NOT EXISTS proteins (
    locus_tag      TEXT PRIMARY KEY,
    genome_acc     TEXT,
    protein_acc    TEXT,
    gene_name      TEXT,
    product        TEXT,
    contig         TEXT,
    start          INTEGER,
    end            INTEGER,
    strand         TEXT,
    protein_len    INTEGER
);
CREATE INDEX IF NOT EXISTS idx_prot_genome  ON proteins (genome_acc);
CREATE INDEX IF NOT EXISTS idx_prot_product ON proteins (product);
CREATE INDEX IF NOT EXISTS idx_prot_contig  ON proteins (genome_acc, contig, start);
"""


def _init_db(db_path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(str(db_path))
    for stmt in DB_SCHEMA.strip().split(';'):
        stmt = stmt.strip()
        if stmt:
            con.execute(stmt)
    # Migrate older DBs that might be missing columns
    for col, coltype in [('seq_count', 'INTEGER')]:
        try:
            con.execute(f'ALTER TABLE genomes ADD COLUMN {col} {coltype}')
        except Exception:
            pass
    con.commit()
    return con


def _parse_gff3(gff_path: Path, genome_acc: str) -> list[dict]:
    rows: list[dict] = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            attrs = dict(_ATTR_RE.findall(parts[8]))
            lt = attrs.get('locus_tag', '')
            if not lt:
                continue
            if attrs.get('pseudo', '').lower() == 'true':
                continue
            rows.append({
                'locus_tag':      lt,
                'genome_acc':     genome_acc,
                'protein_acc':    attrs.get('protein_id', ''),
                'gene_name':      attrs.get('gene', ''),
                'product':        attrs.get('product', 'hypothetical protein'),
                'contig':         parts[0],
                'start':          int(parts[3]),
                'end':            int(parts[4]),
                'strand':         parts[6],
                'protein_len':    (int(parts[4]) - int(parts[3]) + 1) // 3 - 1,
            })
    return rows


def _register_genomes(
    con: sqlite3.Connection,
    downloaded: dict[str, dict[str, Path]],
    meta_map: dict[str, dict],
) -> tuple[int, int, list[tuple[str, str]]]:
    """Insert downloaded genomes + proteins into the DB."""
    cur = con.cursor()
    n_genomes = n_proteins = 0
    skipped: list[tuple[str, str]] = []

    for acc, files in downloaded.items():
        if 'gff' not in files:
            skipped.append((acc, 'no GFF'))
            continue

        meta = meta_map.get(acc, {})
        cur.execute(
            '''INSERT OR IGNORE INTO genomes
               (accession, assembly_name, organism_name, species, strain,
                assembly_level, seq_count, total_length, protein_count, download_date)
               VALUES (?,?,?,?,?,?,?,?,?,?)''',
            (acc,
             meta.get('assembly_name', ''),
             meta.get('organism_name', ''),
             meta.get('species', ''),
             meta.get('strain', ''),
             meta.get('assembly_level', ''),
             meta.get('seq_count'),
             meta.get('total_length', 0),
             0,
             TODAY),
        )
        if cur.rowcount == 0:
            cur.execute(
                'UPDATE genomes SET seq_count=? WHERE accession=? AND seq_count IS NULL',
                (meta.get('seq_count'), acc),
            )
            skipped.append((acc, 'already in DB (files repaired)'))
            continue
        n_genomes += 1

        proteins = _parse_gff3(files['gff'], acc)
        if not proteins:
            skipped.append((acc, 'zero proteins after GFF parse'))
            continue

        cur.executemany(
            '''INSERT OR IGNORE INTO proteins
               (locus_tag, genome_acc, protein_acc, gene_name, product,
                contig, start, end, strand, protein_len)
               VALUES (:locus_tag,:genome_acc,:protein_acc,:gene_name,:product,
                       :contig,:start,:end,:strand,:protein_len)''',
            proteins,
        )
        n_proteins += len(proteins)
        cur.execute(
            'UPDATE genomes SET protein_count=? WHERE accession=?',
            (len(proteins), acc),
        )

        if n_genomes % 50 == 0:
            con.commit()
            print(f'  {n_genomes:,} genomes registered, {n_proteins:,} proteins...')

    con.commit()
    return n_genomes, n_proteins, skipped


# ── combined_proteins.faa ──────────────────────────────────────────────────────
def _build_combined_faa(
    faa_path: Path,
    con: sqlite3.Connection,
    gff_dir: Path,
    fa_dir: Path,
    new_downloaded: dict[str, dict[str, Path]],
) -> int:
    """Build or append to *faa_path*."""
    if not faa_path.exists():
        print(f'  Full rebuild from all {faa_path.name} DB entries...')
        faa_path.parent.mkdir(parents=True, exist_ok=True)
        all_accs = {r[0] for r in con.execute('SELECT accession FROM genomes')}
        total = 0
        for acc in sorted(all_accs):
            gff_p    = gff_dir / f'{acc}.gff3'
            fna_cands = sorted(fa_dir.glob(f'{acc}*_genomic.fna.gz'))
            faa_cands = sorted(fa_dir.glob(f'{acc}*_protein.faa.gz'))
            if not gff_p.exists() or not fna_cands:
                continue
            pid_map = build_pid_to_lt_map(gff_p) if faa_cands else None
            seqs    = build_genome_faa(gff_p, fna_cands[0],
                                       faa_cands[0] if faa_cands else None,
                                       pid_map)
            total  += append_to_faa(seqs, faa_path)
        return total

    # Incremental append
    total = 0
    for acc, files in new_downloaded.items():
        gff_p = files.get('gff')
        fna_p = files.get('fna')
        faa_p = files.get('faa')
        if not gff_p or not fna_p:
            continue
        pid_map = build_pid_to_lt_map(gff_p) if faa_p else None
        seqs    = build_genome_faa(gff_p, fna_p, faa_p, pid_map)
        total  += append_to_faa(seqs, faa_path)
    return total


# ── settings.toml integration ──────────────────────────────────────────────────────
def _add_target_to_settings(
    settings_path: Path,
    target_name: str,
    db_path: Path,
    gff_dir: Path,
    fa_dir: Path,
    faa_path: Path,
) -> None:
    """Append a [[targets]] block to *settings_path* if the target is absent."""
    text = settings_path.read_text()

    if f'name    = "{target_name}"' in text or f"name = '{target_name}'" in text:
        print(f'  [skip] target "{target_name}" already in {settings_path.name}')
        return

    # Try to derive paths relative to [paths].root
    root_m = re.search(r'root\s*=\s*"([^"]+)"', text)
    root   = Path(root_m.group(1)) if root_m else None

    def _rel(p: Path) -> str:
        if root:
            try:
                return str(p.relative_to(root))
            except ValueError:
                pass
        return str(p)

    block = (
        f'\n[[targets]]\n'
        f'name    = "{target_name}"\n'
        f'db      = "{_rel(db_path)}"\n'
        f'gff_dir = "{_rel(gff_dir)}"\n'
        f'fna_dir = "{_rel(fa_dir)}"\n'
        f'faa     = "{_rel(faa_path)}"\n'
    )

    # Insert before the first non-targets section after the targets blocks, else append.
    insert_re = re.compile(
        r'^(\[(?!targets\b|\[))',  # first top-level section that is not [[targets]]
        re.MULTILINE,
    )
    m = insert_re.search(text)
    if m:
        pos  = m.start()
        text = text[:pos] + block + '\n' + text[pos:]
    else:
        text = text.rstrip('\n') + '\n' + block

    settings_path.write_text(text)
    print(f'  Added target "{target_name}" to {settings_path}')


# ── CLI ──────────────────────────────────────────────────────
def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description='Download & index genomes for LocusRuler.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EPILOG,
    )
    ap.add_argument(
        '--taxon', nargs='+', required=True, metavar='NAME',
        help='NCBI taxon name(s), e.g. "Genus species" or "GenusA GenusB"',
    )
    ap.add_argument(
        '--outdir', required=True, metavar='PATH',
        help='Output directory; sub-dirs gff/ and genomes/ are created automatically.',
    )
    ap.add_argument(
        '--db-name', default=None, metavar='NAME',
        help='Database / target name (no extension).  Default: slug of first taxon.',
    )
    ap.add_argument(
        '--seq-count', type=int, default=100, metavar='N',
        help='Max contigs per genome; 0 = no filter.  Default: 100.',
    )
    ap.add_argument(
        '--assembly-level', nargs='+', default=list(ASSEMBLY_LEVELS),
        choices=list(ASSEMBLY_LEVELS), metavar='LEVEL',
        help='Assembly levels to accept: '
             + ', '.join(ASSEMBLY_LEVELS) + '.  Default: all four.',
    )
    ap.add_argument(
        '--init-from', default=None, metavar='PATH',
        help='Copy an existing DB to --outdir before downloading (incremental update).',
    )
    ap.add_argument(
        '--add-target', default=None, metavar='SETTINGS_TOML',
        help='settings.toml to update with a new [[targets]] block after setup.',
    )
    ap.add_argument(
        '--batch-size', type=int, default=50, metavar='N',
        help='Assemblies per ncbi-genome-download call.  Default: 50.',
    )
    ap.add_argument(
        '--dry-run', action='store_true',
        help='Query NCBI and print stats; do not download anything.',
    )
    return ap.parse_args()


def main() -> None:
    args = _parse_args()

    outdir   = Path(args.outdir).resolve()
    db_name  = args.db_name or _slug(args.taxon[0])
    db_path  = outdir / f'{db_name}.db'
    gff_dir  = outdir / 'gff'
    fa_dir   = outdir / 'genomes'
    tmp_dir  = outdir / '_download_tmp'
    faa_path = outdir / 'combined_proteins.faa'
    seq_max  = args.seq_count

    # ── --init-from 
    if args.init_from:
        src = Path(args.init_from)
        if not src.exists():
            sys.exit(f'[ERROR] --init-from path not found: {src}')
        if db_path.exists():
            print(f'[INFO] {db_path.name} already exists — skipping --init-from copy')
        else:
            db_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, db_path)
            print(f'[INFO] {src} → {db_path}')

    # ── Create directories 
    for d in (outdir, gff_dir, fa_dir, tmp_dir):
        d.mkdir(parents=True, exist_ok=True)

    # ── DB init 
    con = _init_db(db_path)
    try:
        existing_accs = {r[0] for r in con.execute('SELECT accession FROM genomes')}
    except Exception:
        existing_accs = set()
    print(f'DB: {db_path.name}  existing genomes: {len(existing_accs):,}')

    # ── NCBI query 
    levels = [lv for lv in ASSEMBLY_LEVELS if lv in args.assembly_level]
    cap = f'≤ {seq_max} contigs per genome' if seq_max else 'any number of contigs'
    print(f'\nQuerying NCBI ({", ".join(levels)}; {cap})...')
    all_targets: list[dict] = []
    for taxon in args.taxon:
        hits = query_ncbi(taxon, seq_max, existing_accs, fa_dir, levels)
        print(f'  {taxon}: {len(hits)} new / incomplete')
        all_targets.extend(hits)

    # ── GCF/GCA deduplication 
    existing_gcf_bases = {_base_id(a) for a in existing_accs if a.startswith('GCF_')}
    query_gcf_bases    = {_base_id(t['accession']) for t in all_targets
                          if t['accession'].startswith('GCF_')}
    all_gcf_bases = existing_gcf_bases | query_gcf_bases

    deduped: list[dict] = []
    for t in all_targets:
        if t['accession'].startswith('GCA_') and _base_id(t['accession']) in all_gcf_bases:
            continue  # a GCF version exists — prefer that
        deduped.append(t)
    all_targets = deduped

    gcf_targets = [t for t in all_targets if t['accession'].startswith('GCF_')]
    gca_targets = [t for t in all_targets if t['accession'].startswith('GCA_')]

    print(f'\nNew / incomplete (after dedup): {len(all_targets)}')
    print(f'  GCF (RefSeq): {len(gcf_targets)}')
    print(f'  GCA (FTP):    {len(gca_targets)}')

    # ── Dry-run summary 
    if args.dry_run:
        if all_targets:
            try:
                import pandas as pd  # optional dependency for pretty tables
                df     = pd.DataFrame(all_targets)
                bins   = [0, 1, 5, 10, 20, 50, 100, 9999]
                labels = ['1', '2-5', '6-10', '11-20', '21-50', '51-100', '101+']
                df['contigs'] = pd.cut(df['seq_count'], bins=bins, labels=labels,
                                       right=True)
                print('\nAssembly level × contigs per genome:')
                print(df.groupby(['assembly_level', 'contigs'], observed=True)
                         .size().unstack(fill_value=0).to_string())
                print('\nTop species:')
                print(df['species'].value_counts().head(15).to_string())
            except ImportError:
                lv_cnt = Counter(t['assembly_level'] for t in all_targets)
                for lv, n in sorted(lv_cnt.items()):
                    print(f'  {lv}: {n}')
        print('\n[dry-run] No files downloaded.')
        con.close()
        return

    if not all_targets:
        print('Nothing to download.  Done.')
        con.close()
        return

    # ── Download 
    print('\nDownloading...')
    if gcf_targets:
        print(f'GCF via ncbi-genome-download ({len(gcf_targets)} assemblies)...')
        _download_batch_gcf([t['accession'] for t in gcf_targets],
                            tmp_dir, args.batch_size)

    if gca_targets:
        print(f'GCA via NCBI FTP ({len(gca_targets)} assemblies)...')
        for i, t in enumerate(gca_targets):
            _download_gca_ftp(t['accession'], tmp_dir)
            if (i + 1) % 5 == 0 or i + 1 == len(gca_targets):
                print(f'  {i + 1}/{len(gca_targets)} processed')
            time.sleep(0.5)

    # ── Organize files 
    print('\nOrganising files...')
    downloaded = _organize_files(tmp_dir, fa_dir, gff_dir)
    n_fna = sum(1 for v in downloaded.values() if 'fna' in v)
    n_gff = sum(1 for v in downloaded.values() if 'gff' in v)
    n_faa = sum(1 for v in downloaded.values() if 'faa' in v)
    print(f'Collected: {len(downloaded)} genomes  FNA {n_fna}  GFF {n_gff}  FAA {n_faa}')

    # ── Register in DB 
    print('\nRegistering in DB...')
    meta_map = {t['accession']: t for t in all_targets}
    n_g, n_p, skipped = _register_genomes(con, downloaded, meta_map)
    print(f'Registered: {n_g:,} genomes, {n_p:,} proteins')
    if skipped:
        reason_cnt = Counter(r for _, r in skipped)
        print('Skipped: ' + ', '.join(f'{n} ({r})' for r, n in reason_cnt.items()))

    # ── combined_proteins.faa 
    print(f'\nBuilding {faa_path.name}...')
    total_seqs = _build_combined_faa(faa_path, con, gff_dir, fa_dir, downloaded)
    if total_seqs:
        print(f'  → {total_seqs:,} sequences written')
    else:
        print('  → no new sequences (FAA unchanged)')

    # ── Summary 
    rows = con.execute(
        'SELECT assembly_level, COUNT(*) FROM genomes GROUP BY assembly_level'
    ).fetchall()
    total_db = sum(r[1] for r in rows)
    print(f'\n=== {db_path.name} ({total_db:,} total genomes) ===')
    for level, n in sorted(rows):
        print(f'  {level}: {n}')
    con.close()

    # ── settings.toml update 
    if args.add_target:
        settings_path = Path(args.add_target).resolve()
        if not settings_path.exists():
            print(f'[WARN] --add-target path not found: {settings_path}')
        else:
            print(f'\nUpdating {settings_path.name}...')
            _add_target_to_settings(
                settings_path, db_name, db_path, gff_dir, fa_dir, faa_path
            )

    print(f'\nDone.  Next step:')
    print(f'  python locus_ruler/run.py --config <locus.json> '
          f'--settings <settings.toml> --target {db_name}')


if __name__ == '__main__':
    main()
