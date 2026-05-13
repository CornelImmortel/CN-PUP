#!/usr/bin/env python3
"""Build long table and mutation matrices from final per-cell VCFs."""

import argparse
import csv
import gzip
import re
from pathlib import Path


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", newline="")


def parse_info(info):
    out = {}
    if not info or info == ".":
        return out
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = True
    return out


def parse_sample(format_field, sample_field):
    keys = format_field.split(":") if format_field else []
    values = sample_field.split(":") if sample_field else []
    return {k: values[i] if i < len(values) else "" for i, k in enumerate(keys)}


def parse_ad(ad):
    if not ad or ad == ".":
        return "", ""
    parts = ad.split(",")
    ref = parts[0] if len(parts) >= 1 else ""
    alt = parts[1] if len(parts) >= 2 else ""
    return ref, alt


def gene_from_info(info_dict):
    ann = info_dict.get("ANN", "")
    if ann:
        genes = []
        ensg = []
        consequences = []
        impacts = []
        for entry in ann.split(","):
            fields = entry.split("|")
            if len(fields) > 1 and fields[1]:
                consequences.append(fields[1])
            if len(fields) > 2 and fields[2]:
                impacts.append(fields[2])
            if len(fields) > 3 and fields[3]:
                genes.append(fields[3])
            if len(fields) > 4 and fields[4]:
                ensg.append(fields[4])
        return ";".join(sorted(set(genes))), ";".join(sorted(set(ensg))), ";".join(sorted(set(consequences))), ";".join(sorted(set(impacts)))
    return "", "", "", ""


def infer_cell_id(path):
    name = Path(path).name
    for suffix in [
        ".monovar.population_cosmic.vcf.gz",
        ".monovar.population_cosmic.vcf",
        ".population_cosmic.vcf.gz",
        ".population_cosmic.vcf",
    ]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return re.sub(r"\.vcf(\.gz)?$", "", name)


def read_vcf(path, patient_id, caller):
    cell_id = infer_cell_id(path)
    rows = []
    with open_text(path) as handle:
        sample_name = cell_id
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#CHROM"):
                parts = line.split("\t")
                if len(parts) > 9:
                    sample_name = parts[9]
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info, fmt, sample = fields[:10]
            info_dict = parse_info(info)
            sample_dict = parse_sample(fmt, sample)
            refread, altread = parse_ad(sample_dict.get("AD", ""))
            gene_symbol, gene_ensg, consequence, impact = gene_from_info(info_dict)
            var_id = f"{chrom}:{pos}_{ref}/{alt}"
            rows.append(
                {
                    "patient_id": patient_id,
                    "ctc_id": cell_id,
                    "cell_id": cell_id,
                    "sample_name": sample_name,
                    "caller": caller,
                    "var_id": var_id,
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": qual,
                    "FILTER": filt,
                    "GT": sample_dict.get("GT", ""),
                    "DP": sample_dict.get("DP", ""),
                    "AD": sample_dict.get("AD", ""),
                    "refread": refread,
                    "altread": altread,
                    "GQ": sample_dict.get("GQ", ""),
                    "gene_symbol": gene_symbol,
                    "gene_ensg": gene_ensg,
                    "consequence": consequence,
                    "impact": impact,
                }
            )
    return rows


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_matrix(path, var_ids, cell_ids, values):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["var_id"] + cell_ids)
        for var_id in var_ids:
            writer.writerow([var_id] + [values.get((var_id, cell_id), 0) for cell_id in cell_ids])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--caller", default="monovar")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("vcfs", nargs="+")
    args = parser.parse_args()

    rows = []
    for vcf in args.vcfs:
        rows.extend(read_vcf(vcf, args.patient_id, args.caller))

    rows.sort(key=lambda r: (r["ctc_id"], r["CHROM"], int(r["POS"]), r["REF"], r["ALT"]))
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    long_fields = [
        "patient_id", "ctc_id", "cell_id", "sample_name", "caller", "var_id", "CHROM", "POS", "REF", "ALT",
        "QUAL", "FILTER", "GT", "DP", "AD", "refread", "altread", "GQ", "gene_symbol", "gene_ensg", "consequence", "impact",
    ]
    write_tsv(out_prefix.with_suffix(".long.tsv"), rows, long_fields)

    var_ids = sorted({r["var_id"] for r in rows})
    cell_ids = sorted({r["ctc_id"] for r in rows})
    binary = {(r["var_id"], r["ctc_id"]): 1 for r in rows}
    alt = {(r["var_id"], r["ctc_id"]): r["altread"] for r in rows}
    ref = {(r["var_id"], r["ctc_id"]): r["refread"] for r in rows}

    write_matrix(out_prefix.with_suffix(".binary_matrix.tsv"), var_ids, cell_ids, binary)
    write_matrix(out_prefix.with_suffix(".altread_matrix.tsv"), var_ids, cell_ids, alt)
    write_matrix(out_prefix.with_suffix(".refread_matrix.tsv"), var_ids, cell_ids, ref)

    summary_rows = []
    for cell_id in cell_ids:
        cell_rows = [r for r in rows if r["ctc_id"] == cell_id]
        summary_rows.append(
            {
                "patient_id": args.patient_id,
                "caller": args.caller,
                "ctc_id": cell_id,
                "variant_count": len(cell_rows),
                "gene_count": len({g for r in cell_rows for g in r["gene_symbol"].split(";") if g}),
            }
        )
    write_tsv(out_prefix.with_suffix(".summary.tsv"), summary_rows, ["patient_id", "caller", "ctc_id", "variant_count", "gene_count"])


if __name__ == "__main__":
    main()
