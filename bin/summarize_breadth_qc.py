#!/usr/bin/env python3
"""Summarize breadth-of-coverage (fraction of target bases at >=Nx) per cell from mosdepth dist files."""

import argparse
import csv
import re
from pathlib import Path

CELL_ID_RE = re.compile(r"^(.*)\.mosdepth\.(region|global)\.dist\.txt$")


def parse_dist_file(path):
    depth_to_frac = {}
    with open(path, newline="") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            chrom, depth, frac = parts
            if chrom != "total":
                continue
            try:
                depth_to_frac[int(depth)] = float(frac)
            except ValueError:
                continue
    return depth_to_frac


def breadth_at(depth_to_frac, threshold):
    if not depth_to_frac:
        return ""
    if threshold in depth_to_frac:
        return depth_to_frac[threshold]
    if threshold > max(depth_to_frac):
        return 0.0
    return ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--thresholds", default="1,10,20", help="Comma-separated depth thresholds")
    parser.add_argument("mosdepth_files", nargs="+")
    args = parser.parse_args()

    thresholds = [int(t) for t in args.thresholds.split(",") if t != ""]

    by_cell = {}
    for f in args.mosdepth_files:
        name = Path(f).name
        m = CELL_ID_RE.match(name)
        if not m:
            continue
        cell_id, kind = m.group(1), m.group(2)
        # Prefer the target-restricted "region" distribution over "global" if both exist.
        if cell_id in by_cell and by_cell[cell_id][1] == "region":
            continue
        by_cell[cell_id] = (f, kind)

    rows = []
    for cell_id in sorted(by_cell):
        path, kind = by_cell[cell_id]
        depth_to_frac = parse_dist_file(path)
        row = {"cell_id": cell_id, "coverage_scope": "on_target_panel" if kind == "region" else "whole_genome"}
        for t in thresholds:
            row[f"breadth_{t}x"] = breadth_at(depth_to_frac, t)
        rows.append(row)

    out_path = Path(args.out_prefix).with_suffix(".breadth_summary.tsv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["cell_id", "coverage_scope"] + [f"breadth_{t}x" for t in thresholds]
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
