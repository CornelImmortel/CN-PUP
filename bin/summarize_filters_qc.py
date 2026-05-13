#!/usr/bin/env python3
"""Write QC summaries for MonoVar split VCFs and final matrices."""

import argparse
import csv
import gzip
from pathlib import Path
from statistics import median


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", newline="")


def parse_sample(format_field, sample_field):
    keys = format_field.split(":") if format_field else []
    values = sample_field.split(":") if sample_field else []
    return {k: values[i] if i < len(values) else "" for i, k in enumerate(keys)}


def parse_ad(ad):
    if not ad or ad == ".":
        return None, None
    parts = ad.split(",")
    try:
        ref = int(parts[0]) if len(parts) >= 1 and parts[0] not in ("", ".") else None
        alt = int(parts[1]) if len(parts) >= 2 and parts[1] not in ("", ".") else None
        return ref, alt
    except ValueError:
        return None, None


def quantile(values, q):
    if not values:
        return ""
    values = sorted(values)
    idx = (len(values) - 1) * q
    lo = int(idx)
    hi = min(lo + 1, len(values) - 1)
    frac = idx - lo
    return values[lo] * (1 - frac) + values[hi] * frac


def summarize_split_vcf(path):
    cell_id = Path(path).name.replace(".monovar.split.vcf", "")
    total_sites = 0
    nonref_sites = 0
    dp_values = []
    alt_values = []
    vaf_values = []
    nonref_dp_values = []
    nonref_alt_values = []
    nonref_vaf_values = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            total_sites += 1
            sample = parse_sample(fields[8], fields[9])
            gt = sample.get("GT", "")
            dp = sample.get("DP", "")
            try:
                dp_int = int(dp) if dp not in ("", ".") else None
            except ValueError:
                dp_int = None
            ref, alt = parse_ad(sample.get("AD", ""))
            if dp_int is not None:
                dp_values.append(dp_int)
            if alt is not None:
                alt_values.append(alt)
            if ref is not None and alt is not None and (ref + alt) > 0:
                vaf_values.append(alt / (ref + alt))
            nonref = gt not in ("0/0", "0|0", "./.", ".|.", "")
            if nonref:
                nonref_sites += 1
                if dp_int is not None:
                    nonref_dp_values.append(dp_int)
                if alt is not None:
                    nonref_alt_values.append(alt)
                if ref is not None and alt is not None and (ref + alt) > 0:
                    nonref_vaf_values.append(alt / (ref + alt))

    def stats(prefix, values):
        return {
            f"{prefix}_min": min(values) if values else "",
            f"{prefix}_q25": quantile(values, 0.25) if values else "",
            f"{prefix}_median": median(values) if values else "",
            f"{prefix}_q75": quantile(values, 0.75) if values else "",
            f"{prefix}_max": max(values) if values else "",
        }

    row = {
        "cell_id": cell_id,
        "total_sites": total_sites,
        "nonref_sites": nonref_sites,
    }
    row.update(stats("all_dp", dp_values))
    row.update(stats("all_alt", alt_values))
    row.update(stats("all_vaf", vaf_values))
    row.update(stats("nonref_dp", nonref_dp_values))
    row.update(stats("nonref_alt", nonref_alt_values))
    row.update(stats("nonref_vaf", nonref_vaf_values))
    return row


def read_final_summary(path):
    if not path or not Path(path).exists():
        return {}
    out = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[row["ctc_id"]] = row
    return out


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--caller", default="monovar")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--final-summary")
    parser.add_argument("--min-total-depth", required=True)
    parser.add_argument("--min-alt-reads", required=True)
    parser.add_argument("--min-vaf", required=True)
    parser.add_argument("--max-population-af", required=True)
    parser.add_argument("--keep-cosmic", required=True)
    parser.add_argument("--germline-mode", required=True)
    parser.add_argument("--vep-sites-per-chunk", required=True)
    parser.add_argument("split_vcfs", nargs="+")
    args = parser.parse_args()

    split_rows = [summarize_split_vcf(p) for p in args.split_vcfs]
    split_rows.sort(key=lambda r: r["cell_id"])
    final = read_final_summary(args.final_summary)
    for row in split_rows:
        cell_id = row["cell_id"]
        row["final_variant_count"] = final.get(cell_id, {}).get("variant_count", 0)
        row["final_gene_count"] = final.get(cell_id, {}).get("gene_count", 0)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fields = list(split_rows[0].keys()) if split_rows else ["cell_id"]
    write_tsv(out_prefix.with_suffix(".prefilter_depth_alt_qc.tsv"), split_rows, fields)

    settings = [
        {"setting": "patient_id", "value": args.patient_id},
        {"setting": "caller", "value": args.caller},
        {"setting": "germline_mode", "value": args.germline_mode},
        {"setting": "min_total_depth", "value": args.min_total_depth},
        {"setting": "min_alt_reads", "value": args.min_alt_reads},
        {"setting": "min_vaf", "value": args.min_vaf},
        {"setting": "max_population_af", "value": args.max_population_af},
        {"setting": "keep_cosmic", "value": args.keep_cosmic},
        {"setting": "vep_sites_per_chunk", "value": args.vep_sites_per_chunk},
    ]
    write_tsv(out_prefix.with_suffix(".filter_settings.tsv"), settings, ["setting", "value"])

    total_prefilter_nonref = sum(int(r["nonref_sites"]) for r in split_rows)
    total_final = sum(int(r["final_variant_count"]) for r in split_rows)
    impact = [
        {"metric": "split_vcfs", "value": len(split_rows)},
        {"metric": "prefilter_nonref_calls", "value": total_prefilter_nonref},
        {"metric": "final_retained_calls", "value": total_final},
        {"metric": "removed_calls", "value": total_prefilter_nonref - total_final},
    ]
    write_tsv(out_prefix.with_suffix(".filter_impact.tsv"), impact, ["metric", "value"])


if __name__ == "__main__":
    main()
