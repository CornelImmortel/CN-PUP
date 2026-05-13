#!/usr/bin/env python3
"""Filter a VCF using VEP population AF annotations, preserving COSMIC hits."""

import argparse
import csv
import gzip
import re
from pathlib import Path

COSMIC_RE = re.compile(r"\bCOS[MV]\d+\b", re.IGNORECASE)


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, newline="")


def parse_float(value):
    if value in (None, "", "-"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def is_cosmic(existing_variation):
    return bool(COSMIC_RE.search(existing_variation or ""))


def read_vep_keep_set(vep_tsv, max_af_threshold, keep_cosmic):
    keep = set()
    annotations = {}
    with open(vep_tsv, newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("##")), delimiter="\t")
        required = {"Uploaded_variation", "Existing_variation", "MAX_AF"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"VEP TSV missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            var_id = row["Uploaded_variation"]
            max_af = parse_float(row.get("MAX_AF"))
            cosmic = is_cosmic(row.get("Existing_variation", ""))
            low_af = max_af is None or max_af <= max_af_threshold
            row_keep = low_af or (keep_cosmic and cosmic)
            previous = annotations.get(var_id)
            if previous is None or row_keep:
                annotations[var_id] = {**row, "keep_population_cosmic": str(row_keep), "cosmic_hit": str(cosmic)}
            if row_keep:
                keep.add(var_id)
    return keep, annotations


def vcf_var_id(fields):
    chrom, pos, _id, ref, alt = fields[:5]
    return f"{chrom}:{pos}_{ref}/{alt}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--vep", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--output-annotations", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--max-af", type=float, default=0.001)
    parser.add_argument("--keep-cosmic", action="store_true")
    args = parser.parse_args()

    keep, annotations = read_vep_keep_set(args.vep, args.max_af, args.keep_cosmic)
    total = 0
    kept = 0

    with open_text(args.vcf, "rt") as src, open(args.output_vcf, "w", newline="") as out_vcf:
        for line in src:
            if line.startswith("#"):
                out_vcf.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            total += 1
            if vcf_var_id(fields) in keep:
                kept += 1
                out_vcf.write(line)

    fieldnames = [
        "Uploaded_variation", "Location", "Allele", "Gene", "SYMBOL", "Feature",
        "Consequence", "IMPACT", "HGVSc", "HGVSp", "Existing_variation", "MAX_AF",
        "MAX_AF_POPS", "AF", "gnomADe_AF", "gnomADg_AF", "keep_population_cosmic", "cosmic_hit",
    ]
    with open(args.output_annotations, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for var_id in sorted(annotations):
            writer.writerow({k: annotations[var_id].get(k, "") for k in fieldnames})

    with open(args.summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["input_variants", total])
        writer.writerow(["kept_variants", kept])
        writer.writerow(["removed_variants", total - kept])
        writer.writerow(["max_population_af", args.max_af])
        writer.writerow(["keep_cosmic", args.keep_cosmic])


if __name__ == "__main__":
    main()
