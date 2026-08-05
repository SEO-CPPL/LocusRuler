"""Genome FAA building utilities."""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Iterator

# NCBI translation table 11, from ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.
# The codon mapping is the standard genetic code; what makes it table 11 is the
# initiator set below, applied by translate_dna.
CODON_TABLE: dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
# Table 11 initiators. Any of these translates to M in the first position,
# which is how NCBI writes the protein FASTA this fallback has to match.
START_CODONS = frozenset({'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'})

_RC_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')
_ATTR_RE  = re.compile(r'(\w+)=([^;]+)')


def revcomp(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def translate_dna(dna: str, cds: bool = True) -> str:
    """Translate a CDS the way NCBI table 11 does, initiator codon included."""
    aa: list[str] = []
    for i in range(0, len(dna) - 2, 3):
        aa.append(CODON_TABLE.get(dna[i:i + 3].upper(), 'X'))
    if cds and aa and dna[:3].upper() in START_CODONS:
        aa[0] = 'M'
    return ''.join(aa).rstrip('*').replace('*', 'X')


def load_fna(fna_path: Path) -> dict[str, str]:
    """Return {contig_id: sequence} from a (possibly gzipped) FASTA."""
    contigs: dict[str, str] = {}
    open_fn = gzip.open if str(fna_path).endswith('.gz') else open
    header: str | None = None
    parts: list[str] = []
    with open_fn(fna_path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    contigs[header] = ''.join(parts)
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        contigs[header] = ''.join(parts)
    return contigs


def build_pid_to_lt_map(gff_path: Path) -> dict[str, str]:
    """Return {protein_id: locus_tag} from a GFF3 file."""
    mapping: dict[str, str] = {}
    open_fn = gzip.open if str(gff_path).endswith('.gz') else open
    with open_fn(gff_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            attrs = dict(_ATTR_RE.findall(parts[8]))
            lt = attrs.get('locus_tag', '') or attrs.get('ID', '')
            if not lt:
                continue
            pid = attrs.get('protein_id', '')
            if pid:
                mapping[pid] = lt
            feat_id = attrs.get('ID', '')
            if feat_id and feat_id.startswith('fig|'):
                mapping[feat_id] = lt
    return mapping


def iter_faa_renamed(
    faa_path: Path, pid_to_lt: dict[str, str]
) -> Iterator[tuple[str, str]]:
    """Yield (locus_tag, sequence) from a protein FAA, renaming headers."""
    open_fn = gzip.open if str(faa_path).endswith('.gz') else open
    current_lt: str | None = None
    seq_lines: list[str] = []
    with open_fn(faa_path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_lt and seq_lines:
                    yield current_lt, ''.join(seq_lines)
                pid = line[1:].split()[0]
                current_lt = pid_to_lt.get(pid)
                seq_lines = []
            elif current_lt:
                seq_lines.append(line)
    if current_lt and seq_lines:
        yield current_lt, ''.join(seq_lines)


def iter_faa_translated(
    gff_path: Path, fna_path: Path
) -> Iterator[tuple[str, str]]:
    """Translate CDS records for an assembly that ships no protein FASTA.

    The result joins the same combined FAA as NCBI-supplied proteins, so it
    has to match how NCBI would have translated the record.
    """
    raw = load_fna(fna_path)
    if not raw:
        return
    ctg_cache: dict[str, str] = {}
    open_fn = gzip.open if str(gff_path).endswith('.gz') else open
    with open_fn(gff_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            attrs = dict(_ATTR_RE.findall(parts[8]))
            lt = attrs.get('locus_tag', '') or attrs.get('ID', '')
            if not lt:
                continue
            ctg = parts[0]
            if ctg not in ctg_cache:
                clean = ctg.split('|')[-1]
                found = ''
                for hdr, seq in raw.items():
                    if hdr.split()[0] in (ctg, clean) or ctg in hdr or clean in hdr:
                        found = seq
                        break
                ctg_cache[ctg] = found
            dna = ctg_cache[ctg]
            if not dna:
                continue
            start, end, strand = int(parts[3]) - 1, int(parts[4]), parts[6]
            cds = dna[start:end]
            if strand == '-':
                cds = revcomp(cds)
            pep = translate_dna(cds)
            if len(pep) >= 10:
                yield lt, pep


def build_genome_faa(
    gff_path: Path,
    fna_path: Path,
    faa_path: Path | None,
    pid_to_lt: dict[str, str] | None = None,
) -> Iterator[tuple[str, str]]:
    """Yield (locus_tag, sequence) per CDS, preferring NCBI's protein FASTA."""
    if faa_path and faa_path.exists():
        if pid_to_lt is None:
            pid_to_lt = build_pid_to_lt_map(gff_path)
        yield from iter_faa_renamed(faa_path, pid_to_lt)
    else:
        yield from iter_faa_translated(gff_path, fna_path)


def append_to_faa(
    sequences: Iterator[tuple[str, str]],
    out_faa: Path,
    line_width: int = 60,
) -> int:
    """Append *sequences* to *out_faa* in FASTA format.  Returns count written."""
    n = 0
    with open(out_faa, 'a') as fh:
        for lt, seq in sequences:
            fh.write(f'>{lt}\n')
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + '\n')
            n += 1
    return n
