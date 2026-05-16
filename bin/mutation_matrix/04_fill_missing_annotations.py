#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip
from collections import defaultdict
from pathlib import Path

ANNOT_FIELDS = ["gene_symbol", "gene_id", "effect", "impact", "transcript", "hgvs_c", "hgvs_p"]
MISSING = {"", ".", "NA", "N/A", "None", "Unknown", "unknown", "UNANNOTATED", "unannotated"}
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3, "UNANNOTATED": 4, "": 9}


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")


def open_write(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(path, "wt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open("w", encoding="utf-8", newline="")


def missing(value) -> bool:
    return value is None or str(value).strip() in MISSING


def ann_score(row: dict) -> tuple:
    gene = row.get("gene_symbol", "")
    impact = row.get("impact", "")
    effect = row.get("effect", "")
    return (
        0 if not missing(gene) else 1,
        IMPACT_RANK.get(str(impact).upper(), 8),
        0 if not missing(effect) else 1,
        str(gene),
    )


def combine_unique(values):
    out = []
    seen = set()
    for value in values:
        if missing(value):
            continue
        for part in str(value).split(";"):
            part = part.strip()
            if not part or part in MISSING or part in seen:
                continue
            seen.add(part)
            out.append(part)
    return ";".join(out)


def main():
    ap = argparse.ArgumentParser(description="Fill missing annotations by var_id from rows where another caller has ANN-derived fields.")
    ap.add_argument("--input", required=True, help="Merged long table TSV/TSV.GZ")
    ap.add_argument("--output", required=True, help="Annotated merged long table TSV/TSV.GZ")
    args = ap.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    with open_text(input_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    for field in ANNOT_FIELDS:
        if field not in fieldnames:
            fieldnames.append(field)
            for row in rows:
                row[field] = ""

    by_var = defaultdict(list)
    for row in rows:
        var_id = row.get("var_id", "")
        if var_id:
            by_var[var_id].append(row)

    annotation = {}
    for var_id, var_rows in by_var.items():
        candidates = [r for r in var_rows if any(not missing(r.get(f, "")) for f in ANNOT_FIELDS)]
        if not candidates:
            continue
        best = sorted(candidates, key=ann_score)[0]
        ann = {f: best.get(f, "") for f in ANNOT_FIELDS}
        ann["gene_symbol"] = combine_unique(r.get("gene_symbol", "") for r in candidates) or ann.get("gene_symbol", "")
        ann["gene_id"] = combine_unique(r.get("gene_id", "") for r in candidates) or ann.get("gene_id", "")
        annotation[var_id] = ann

    filled_cells = 0
    filled_rows = set()
    for i, row in enumerate(rows):
        ann = annotation.get(row.get("var_id", ""))
        if not ann:
            continue
        for field in ANNOT_FIELDS:
            if missing(row.get(field, "")) and not missing(ann.get(field, "")):
                row[field] = ann[field]
                filled_cells += 1
                filled_rows.add(i)

    with open_write(output_path) as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Read rows: {len(rows)}")
    print(f"Variants with annotation source: {len(annotation)}")
    print(f"Rows updated: {len(filled_rows)}")
    print(f"Annotation cells filled: {filled_cells}")
    print(f"Wrote: {output_path}")


if __name__ == "__main__":
    main()