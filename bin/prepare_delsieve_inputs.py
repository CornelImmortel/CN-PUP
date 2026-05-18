#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def to_float(value):
    try:
        if value in (None, "", "."):
            return None
        return float(value)
    except ValueError:
        return None


def load_ctc_bams(cell_metadata: Path):
    rows = []
    with open(cell_metadata, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cell_type = row.get("cell_type", "").lower()
            cell_id = row.get("cell_id", "")
            if cell_type == "ctc" or cell_id.upper().startswith("CTC"):
                rows.append((cell_id, row["bam_path"]))
    rows.sort(key=lambda x: (re.sub(r"\d+", "", x[0]), int(re.search(r"\d+", x[0]).group(0)) if re.search(r"\d+", x[0]) else 0, x[0]))
    return rows


def select_candidates(args):
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    allowed_callers = None if args.callers.lower() == "all" else {x.strip() for x in args.callers.split(",") if x.strip()}
    ctc_ids = {cell_id for cell_id, _ in load_ctc_bams(Path(args.cell_metadata))}

    by_var = {}
    support = defaultdict(set)
    callers = defaultdict(set)
    genes = defaultdict(set)
    proteins = defaultdict(set)
    consequences = defaultdict(set)

    with open_text(args.long_table) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sample_id = row.get("sample_id", "")
            caller = row.get("caller", "")
            var_id = row.get("var_id", "")
            if not var_id or sample_id not in ctc_ids:
                continue
            if allowed_callers is not None and caller not in allowed_callers:
                continue
            total = to_float(row.get("total_depth") or row.get("DP"))
            alt = to_float(row.get("altread"))
            vaf = to_float(row.get("vaf"))
            if vaf is None and total and alt is not None:
                vaf = alt / total if total > 0 else None
            if total is None or alt is None:
                continue
            if total < args.min_total_depth or alt < args.min_alt_reads or (vaf or 0.0) < args.min_vaf:
                continue
            support[var_id].add(sample_id)
            callers[var_id].add(caller)
            for field, target in (("gene_symbol", genes), ("hgvs_p", proteins), ("effect", consequences)):
                value = row.get(field, "")
                if value:
                    target[var_id].add(value)
            by_var.setdefault(var_id, row)

    kept = []
    for var_id, row in by_var.items():
        if len(support[var_id]) >= args.min_ctc_support:
            kept.append(row)

    def sort_key(row):
        chrom = row.get("chrom", "")
        pos = int(float(row.get("pos", 0)))
        return (chrom, pos, row.get("ref", ""), row.get("alt", ""))

    kept.sort(key=sort_key)
    with open(outdir / "candidate_sites.tsv", "w", encoding="utf-8", newline="") as out:
        cols = [
            "var_id", "chrom", "pos", "ref", "alt", "ctc_support", "supporting_ctcs",
            "callers", "gene_symbol", "protein_change", "consequence",
        ]
        writer = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for row in kept:
            var_id = row["var_id"]
            writer.writerow({
                "var_id": var_id,
                "chrom": row.get("chrom", ""),
                "pos": row.get("pos", ""),
                "ref": row.get("ref", ""),
                "alt": row.get("alt", ""),
                "ctc_support": len(support[var_id]),
                "supporting_ctcs": ",".join(sorted(support[var_id])),
                "callers": ",".join(sorted(callers[var_id])),
                "gene_symbol": ",".join(sorted(genes[var_id])),
                "protein_change": ",".join(sorted(proteins[var_id])),
                "consequence": ",".join(sorted(consequences[var_id])),
            })

    with open(outdir / "candidate_sites.bed", "w", encoding="utf-8") as bed:
        for row in kept:
            chrom, pos = row.get("chrom", ""), int(float(row.get("pos", 0)))
            bed.write(f"{chrom}\t{pos - 1}\t{pos}\t{row['var_id']}\n")

    cells = load_ctc_bams(Path(args.cell_metadata))
    with open(outdir / "cell_names", "w", encoding="utf-8") as out:
        out.write(" ".join(cell for cell, _ in cells) + "\n")
    with open(outdir / "cell_names.tsv", "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["cell_id", "bam_path"])
        writer.writerows(cells)
    with open(outdir / "bams.txt", "w", encoding="utf-8") as out:
        for _, bam in cells:
            out.write(f"{bam}\n")
    with open(outdir / "candidate_summary.tsv", "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["candidate_sites", len(kept)])
        writer.writerow(["ctc_cells", len(cells)])
        writer.writerow(["min_ctc_support", args.min_ctc_support])
        writer.writerow(["callers", args.callers])


def parse_bases(ref_base: str, bases: str):
    counts = Counter()
    i = 0
    ref = ref_base.upper()
    while i < len(bases):
        c = bases[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            i += 1
            j = i
            while j < len(bases) and bases[j].isdigit():
                j += 1
            n = int(bases[i:j] or 0)
            i = j + n
            continue
        if c in ".,": 
            counts[ref] += 1
        else:
            b = c.upper()
            if b in {"A", "C", "G", "T"}:
                counts[b] += 1
        i += 1
    return counts


def parse_mpileup(args):
    outdir = Path(args.output_dir)
    cells = []
    with open(outdir / "cell_names.tsv", "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cells.append(row["cell_id"])

    site_by_key = {}
    with open(outdir / "candidate_sites.tsv", "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            site_by_key[(row["chrom"], str(int(float(row["pos"]))))] = row

    with open(args.mpileup, "r", encoding="utf-8") as pile, \
            open(outdir / "read_counts.full_support_coverage.tsv", "w", encoding="utf-8", newline="") as counts_out, \
            open(outdir / "read_counts.human_readable.tsv", "w", encoding="utf-8", newline="") as readable_out:
        counts_writer = csv.writer(counts_out, delimiter="\t")
        readable_cols = ["var_id", "chrom", "pos", "ref", "candidate_alt", "cell_id", "A", "C", "G", "T", "alt1", "alt2", "alt3", "ref_count", "coverage"]
        readable_writer = csv.DictWriter(readable_out, fieldnames=readable_cols, delimiter="\t")
        readable_writer.writeheader()
        header = ["chrom", "pos", "ref", "candidate_alt", "var_id"]
        for cell in cells:
            header.extend([f"{cell}.alt1", f"{cell}.alt2", f"{cell}.alt3", f"{cell}.ref", f"{cell}.coverage"])
        counts_writer.writerow(header)

        for line in pile:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref = fields[0], fields[1], fields[2].upper()
            site = site_by_key.get((chrom, pos), {})
            row = [chrom, pos, ref, site.get("alt", ""), site.get("var_id", f"{chrom}:{pos}")]
            sample_fields = fields[3:]
            for idx, cell in enumerate(cells):
                base_idx = idx * 3
                if base_idx + 2 >= len(sample_fields):
                    depth, bases = 0, ""
                else:
                    depth = int(sample_fields[base_idx])
                    bases = sample_fields[base_idx + 1]
                base_counts = parse_bases(ref, bases)
                alt_counts = sorted([base_counts[b] for b in "ACGT" if b != ref], reverse=True)
                ref_count = base_counts[ref]
                row.extend([alt_counts[0], alt_counts[1], alt_counts[2], ref_count, depth])
                readable_writer.writerow({
                    "var_id": site.get("var_id", f"{chrom}:{pos}"),
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "candidate_alt": site.get("alt", ""),
                    "cell_id": cell,
                    "A": base_counts["A"],
                    "C": base_counts["C"],
                    "G": base_counts["G"],
                    "T": base_counts["T"],
                    "alt1": alt_counts[0],
                    "alt2": alt_counts[1],
                    "alt3": alt_counts[2],
                    "ref_count": ref_count,
                    "coverage": depth,
                })
            counts_writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Prepare candidate-site read counts for DelSIEVE from CN-PUP outputs.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("select-candidates")
    p.add_argument("--long-table", required=True)
    p.add_argument("--cell-metadata", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--min-ctc-support", type=int, default=2)
    p.add_argument("--min-total-depth", type=float, default=20)
    p.add_argument("--min-alt-reads", type=float, default=4)
    p.add_argument("--min-vaf", type=float, default=0.0)
    p.add_argument("--callers", default="all")
    p.set_defaults(func=select_candidates)

    p = sub.add_parser("parse-mpileup")
    p.add_argument("--mpileup", required=True)
    p.add_argument("--output-dir", required=True)
    p.set_defaults(func=parse_mpileup)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
