#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip
from pathlib import Path

def to_float(x):
    try:
        if x in (None, "", "."):
            return None
        return float(x)
    except ValueError:
        return None

def read_rows(path):
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        yield from csv.DictReader(fh, delimiter="\t")

def write_matrix(path, variants, columns, values, fill):
    with open(path, "w", encoding="utf-8", newline="") as out:
        out.write("var_id\t" + "\t".join(columns) + "\n")
        for var in variants:
            out.write(var + "\t" + "\t".join(str(values.get((var, col), fill)) for col in columns) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--matrix-dir", required=True)
    ap.add_argument("--ctc-scite-dir", required=True)
    args = ap.parse_args()
    matrix_dir, scite_dir = Path(args.matrix_dir), Path(args.ctc_scite_dir)
    matrix_dir.mkdir(parents=True, exist_ok=True); scite_dir.mkdir(parents=True, exist_ok=True)
    records = {}
    meta = {}
    samples = set()
    for r in read_rows(args.input):
        sample_col = f"{r['sample_id']}__{r['caller']}"
        samples.add(sample_col)
        key = (r["var_id"], sample_col)
        alt = to_float(r.get("altread"))
        total = to_float(r.get("total_depth"))
        vaf = to_float(r.get("vaf"))
        current = records.get(key)
        score = (-1 if alt is None else alt, -1 if total is None else total)
        if current is None or score > current["score"]:
            records[key] = {"alt": alt, "total": total, "vaf": vaf, "filter": r.get("filter"), "score": score}
        meta.setdefault(r["var_id"], r)
    variants = sorted(meta, key=lambda v: (meta[v].get("chrom", ""), int(float(meta[v].get("pos", 0))), v))
    columns = sorted(samples)
    binary, vafs, alts, totals = {}, {}, {}, {}
    for var in variants:
        for col in columns:
            rec = records.get((var, col))
            if rec is None:
                binary[(var, col)] = -1
                continue
            alt = rec["alt"] or 0
            total = rec["total"] or 0
            vf = rec["vaf"] if rec["vaf"] is not None else (alt / total if total > 0 else None)
            # Depth/alt/VAF no longer decide keep-vs-drop upstream (see
            # FILTER_MONOVAR_AND_SUBTRACT_GERMLINE in main.nf) -- the VCF
            # FILTER column (PASS vs LOWCONF) is the single source of truth
            # for confidence, carried straight through by 03_vcf_to_long_table.py.
            binary[(var, col)] = 1 if rec.get("filter") == "PASS" else 0
            vafs[(var, col)] = "" if vf is None else round(vf, 6)
            alts[(var, col)] = int(alt)
            totals[(var, col)] = int(total)
    write_matrix(matrix_dir / "mutation_binary.tsv", variants, columns, binary, -1)
    write_matrix(matrix_dir / "vaf.tsv", variants, columns, vafs, "")
    write_matrix(matrix_dir / "alt_reads.tsv", variants, columns, alts, 0)
    write_matrix(matrix_dir / "total_depth.tsv", variants, columns, totals, 0)

    # CTC-SCITE total/alt format without header: chr pos ref alt sample_total sample_alt ...
    with open(scite_dir / "ctc_scite_input_snv_total.tsv", "w", encoding="utf-8", newline="") as out:
        for var in variants:
            m = meta[var]
            row = [m["chrom"], m["pos"], m["ref"], m["alt"]]
            for col in columns:
                row += [str(totals.get((var, col), 0)), str(alts.get((var, col), 0))]
            out.write("\t".join(row) + "\n")
    with open(scite_dir / "sample_description_graphviz.tsv", "w", encoding="utf-8", newline="") as out:
        for col in columns:
            out.write(f"{col}\t1\t1\n")
    print(f"Wrote matrices for {len(variants)} variants and {len(columns)} sample/caller columns")
if __name__ == "__main__":
    main()
