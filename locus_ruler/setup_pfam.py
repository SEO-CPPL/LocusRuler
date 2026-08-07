#!/usr/bin/env python3
"""Download and index Pfam-A.hmm for optional domain recovery."""

import argparse
import gzip
import shutil
import subprocess
import urllib.request
from pathlib import Path


PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Download and index Pfam-A.hmm for [domain_recovery]")
    ap.add_argument("--out-dir", default="input/pfam",
                    help="Where to write Pfam-A.hmm and its hmmpress index "
                         "(default: input/pfam)")
    ap.add_argument("--url", default=PFAM_URL,
                    help="Source URL for Pfam-A.hmm.gz (default: the "
                         "current EBI release)")
    ap.add_argument("--force", action="store_true",
                    help="Re-download and re-index even if Pfam-A.hmm "
                         "already exists in --out-dir")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gz_path = out_dir / "Pfam-A.hmm.gz"
    hmm_path = out_dir / "Pfam-A.hmm"

    if not gz_path.exists() or args.force:
        print(f"[setup_pfam] downloading {args.url}")
        urllib.request.urlretrieve(args.url, gz_path)
    else:
        print(f"[setup_pfam] using existing {gz_path}")

    if not hmm_path.exists() or args.force:
        print(f"[setup_pfam] decompressing → {hmm_path}")
        with gzip.open(gz_path, "rb") as fin, open(hmm_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        print(f"[setup_pfam] using existing {hmm_path}")

    hmmpress = shutil.which("hmmpress")
    if hmmpress:
        print("[setup_pfam] hmmpress indexing")
        subprocess.run([hmmpress, "-f", str(hmm_path)], check=True)
    else:
        print("[setup_pfam] hmmpress not found; hmmfetch still works but may be slower")

    print("\nAdd this to settings.toml:")
    print("[domain_recovery]")
    print(f'pfam_a_hmm = "{hmm_path}"')


if __name__ == "__main__":
    main()
