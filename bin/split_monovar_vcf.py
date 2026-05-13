#!/usr/bin/env python3
"""Split a multi-sample MonoVar VCF into one VCF per cell.

MonoVar can add one extra trailing passcode/barcode column after the sample
columns. This script keeps that value as INFO/MONOVAR_PASSCODE in each split
record so no information is silently lost.
"""

import argparse
import csv
import re
from pathlib import Path


def srr_token(path_or_name):
    match = re.search(r"(SRR\d+)", path_or_name or "")
    return match.group(1) if match else None


def read_cells(path, patient_id):
    cells = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"patient_id", "cell_id", "bam_path", "cell_type"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"cell metadata missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            if not row.get("patient_id"):
                continue
            if row["patient_id"] != patient_id:
                continue
            cells.append(
                {
                    "patient_id": row["patient_id"],
                    "cell_id": row["cell_id"],
                    "bam_path": row["bam_path"],
                    "cell_type": row["cell_type"],
                    "srr": srr_token(row.get("bam_path", "")),
                }
            )
    if not cells:
        raise SystemExit(f"no cells for patient_id={patient_id} in {path}")
    return cells


def build_sample_map(vcf_samples, cells):
    used_cells = set()
    mapping = []

    for idx, sample in enumerate(vcf_samples):
        chosen = None
        sample_token = srr_token(sample)
        for cell_index, cell in enumerate(cells):
            if cell_index in used_cells:
                continue
            if cell["srr"] and (cell["srr"] in sample or cell["srr"] == sample_token):
                chosen = (cell_index, cell)
                break
        if chosen is None and idx < len(cells) and idx not in used_cells:
            chosen = (idx, cells[idx])
        if chosen is None:
            raise SystemExit(f"could not map MonoVar VCF sample '{sample}' to a cell metadata row")
        used_cells.add(chosen[0])
        mapping.append({"sample_index": idx, "vcf_sample": sample, **chosen[1]})

    return mapping


def add_info_value(info, key, value):
    if value in (None, "", "."):
        return info
    addition = f"{key}={value}"
    if info in ("", "."):
        return addition
    return f"{info};{addition}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Multi-sample MonoVar VCF")
    parser.add_argument("--cell-metadata", required=True, help="TSV with patient_id, cell_id, bam_path, cell_type")
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    in_vcf = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cells = read_cells(args.cell_metadata, args.patient_id)

    meta_headers = []
    chrom_header = None
    records_seen = 0
    records_written = 0
    out_handles = []
    mapping = None

    with open(in_vcf, "r", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("##"):
                meta_headers.append(line)
                continue
            if line.startswith("#CHROM"):
                chrom_header = line.split("\t")
                if len(chrom_header) < 10:
                    raise SystemExit("MonoVar VCF has no sample columns")
                fixed_header = chrom_header[:9]
                samples = chrom_header[9:]
                mapping = build_sample_map(samples, cells)

                if not any(h.startswith("##INFO=<ID=MONOVAR_PASSCODE,") for h in meta_headers):
                    meta_headers.append('##INFO=<ID=MONOVAR_PASSCODE,Number=1,Type=String,Description="MonoVar trailing per-site passcode/barcode from the multi-sample VCF">')

                for item in mapping:
                    path = outdir / f"{item['cell_id']}.monovar.split.vcf"
                    fh = open(path, "w", newline="")
                    for header in meta_headers:
                        print(header, file=fh)
                    print(f"##SAMPLE=<ID={item['cell_id']},OriginalMonoVarSample={item['vcf_sample']},CellType={item['cell_type']}>", file=fh)
                    print("\t".join(fixed_header + [item["cell_id"]]), file=fh)
                    out_handles.append(fh)
                continue
            if line.startswith("#"):
                continue

            if chrom_header is None or mapping is None:
                raise SystemExit("VCF body found before #CHROM header")

            fields = line.split("\t")
            n_samples = len(mapping)
            expected = 9 + n_samples
            if len(fields) < expected:
                raise SystemExit(f"record has too few columns at line starting: {line[:120]}")

            passcode = None
            if len(fields) > expected:
                passcode = fields[expected]

            fixed = fields[:9]
            if passcode is not None:
                fixed[7] = add_info_value(fixed[7], "MONOVAR_PASSCODE", passcode)

            for idx, fh in enumerate(out_handles):
                sample_field = fields[9 + idx]
                print("\t".join(fixed + [sample_field]), file=fh)
                records_written += 1
            records_seen += 1

    for fh in out_handles:
        fh.close()

    if chrom_header is None:
        raise SystemExit("VCF missing #CHROM header")

    map_path = outdir / f"{args.patient_id}.monovar.split_sample_map.tsv"
    with open(map_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["sample_index", "cell_id", "cell_type", "vcf_sample", "srr", "bam_path"],
            lineterminator="\n",
        )
        writer.writeheader()
        for item in mapping:
            writer.writerow({k: item.get(k, "") for k in writer.fieldnames})

    print(f"Input records: {records_seen}")
    print(f"Split VCFs: {len(mapping)}")
    print(f"Written per-cell records: {records_written}")
    print(f"Sample map: {map_path}")


if __name__ == "__main__":
    main()
