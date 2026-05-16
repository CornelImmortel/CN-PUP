#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip, shutil
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    files = sorted(Path(args.input_dir).glob("*.long.tsv.gz"))
    if not files:
        raise SystemExit(f"No *.long.tsv.gz files found in {args.input_dir}")
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = False
    total = 0
    with gzip.open(out, "wt", encoding="utf-8", newline="") as oh:
        for f in files:
            with gzip.open(f, "rt", encoding="utf-8") as ih:
                header = ih.readline()
                if not wrote_header:
                    oh.write(header); wrote_header = True
                for line in ih:
                    oh.write(line); total += 1
    print(f"Merged {len(files)} files into {out} ({total} rows)")
if __name__ == "__main__":
    main()
