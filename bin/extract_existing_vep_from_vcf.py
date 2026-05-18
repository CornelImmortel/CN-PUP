#!/usr/bin/env python3
"""Extract existing VEP-style VCF annotations into the pipeline VEP TSV format."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


DEFAULT_FIELDS = [
    "Uploaded_variation", "Location", "Allele", "Gene", "SYMBOL", "Feature",
    "Consequence", "IMPACT", "HGVSc", "HGVSp", "Existing_variation", "MAX_AF",
    "MAX_AF_POPS", "AF", "gnomADe_AF", "gnomADg_AF",
]


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else path.open("r", encoding="utf-8", errors="replace")


def parse_info(info: str) -> dict[str, str | bool]:
    out: dict[str, str | bool] = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = True
    return out


def parse_csq_fields(line: str) -> list[str] | None:
    match = re.search(r"Format: ([^\">]+)", line)
    if not match:
        return None
    return [item.strip() for item in match.group(1).split("|")]


def empty_row(uploaded: str, location: str, allele: str) -> dict[str, str]:
    row = {field: "" for field in DEFAULT_FIELDS}
    row["Uploaded_variation"] = uploaded
    row["Location"] = location
    row["Allele"] = allele
    return row


def pick(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = row.get(key, "")
        if value not in ("", "."):
            return value
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    vcf = Path(args.vcf)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    csq_fields: list[str] | None = None
    rows: list[dict[str, str]] = []
    saw_csq = False

    with open_text(vcf) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = parse_csq_fields(line)
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _vid, ref, alts, _qual, _filter, info_s = fields[:8]
            info = parse_info(info_s)
            csq = str(info.get("CSQ", ""))
            if not csq or not csq_fields:
                continue
            saw_csq = True
            for alt in alts.split(","):
                uploaded = f"{chrom}:{pos}_{ref}/{alt}"
                location = f"{chrom}:{pos}"
                matching = []
                for entry in csq.split(","):
                    values = entry.split("|")
                    values += [""] * max(0, len(csq_fields) - len(values))
                    ann = dict(zip(csq_fields, values))
                    if ann.get("Allele") == alt:
                        matching.append(ann)
                if not matching:
                    matching = [dict(zip(csq_fields, entry.split("|") + [""] * len(csq_fields))) for entry in csq.split(",")]
                for ann in matching:
                    row = empty_row(uploaded, location, alt)
                    row["Gene"] = pick(ann, "Gene")
                    row["SYMBOL"] = pick(ann, "SYMBOL", "Gene")
                    row["Feature"] = pick(ann, "Feature")
                    row["Consequence"] = pick(ann, "Consequence")
                    row["IMPACT"] = pick(ann, "IMPACT")
                    row["HGVSc"] = pick(ann, "HGVSc")
                    row["HGVSp"] = pick(ann, "HGVSp")
                    row["Existing_variation"] = pick(ann, "Existing_variation", "Existing_variation_dbSNP")
                    row["MAX_AF"] = pick(ann, "MAX_AF")
                    row["MAX_AF_POPS"] = pick(ann, "MAX_AF_POPS")
                    row["AF"] = pick(ann, "AF")
                    row["gnomADe_AF"] = pick(ann, "gnomADe_AF", "gnomAD_exomes_AF")
                    row["gnomADg_AF"] = pick(ann, "gnomADg_AF", "gnomAD_genomes_AF")
                    rows.append(row)

    with output.open("w", encoding="utf-8", newline="") as out:
        out.write("\t".join(DEFAULT_FIELDS) + "\n")
        for row in rows:
            out.write("\t".join(row.get(field, "") for field in DEFAULT_FIELDS) + "\n")

    if not saw_csq:
        raise SystemExit(f"No VEP CSQ annotation found in {vcf}")


if __name__ == "__main__":
    main()
