#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, itertools
from collections import defaultdict
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")


def read_long(path: Path):
    with open_text(path) as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def main():
    ap = argparse.ArgumentParser(description="Summarize pairwise and all-way overlap in a merged long mutation table.")
    ap.add_argument("--input", required=True, help="Merged long table TSV/TSV.GZ")
    ap.add_argument("--output", required=True, help="Output overlap summary TSV")
    args = ap.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sets = defaultdict(set)
    rows = 0
    for r in read_long(input_path):
        sample_id = r.get("sample_id", "")
        caller = r.get("caller", "")
        var_id = r.get("var_id", "")
        if not sample_id or not caller or not var_id:
            continue
        sets[f"{sample_id}__{caller}"].add(var_id)
        rows += 1

    names = sorted(sets)
    with output_path.open("w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["comparison", "n_a", "n_b", "n_intersection", "jaccard"])
        for name in names:
            w.writerow([name, len(sets[name]), "", len(sets[name]), ""])
        for a, b in itertools.combinations(names, 2):
            inter = len(sets[a] & sets[b])
            union = len(sets[a] | sets[b])
            jac = inter / union if union else 0
            w.writerow([f"{a} vs {b}", len(sets[a]), len(sets[b]), inter, f"{jac:.6f}"])
        if names:
            all_inter = len(set.intersection(*(sets[n] for n in names)))
            all_union = len(set.union(*(sets[n] for n in names)))
            w.writerow(["ALL", "", "", all_inter, f"{(all_inter / all_union if all_union else 0):.6f}"])

    print(f"Read {rows} rows from {input_path}")
    print(f"Wrote overlap QC: {output_path}")
    if names:
        print(f"All-way shared variants: {all_inter}")


if __name__ == "__main__":
    main()