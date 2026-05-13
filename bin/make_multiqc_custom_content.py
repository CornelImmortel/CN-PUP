#!/usr/bin/env python3
"""Convert CN-PUP QC TSVs to MultiQC custom-content JSON files."""

import argparse
import csv
import json
from pathlib import Path


def read_key_value_tsv(path, key_col, value_col):
    out = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[row[key_col]] = row[value_col]
    return out


def read_rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def coerce(value):
    if value is None or value == "":
        return None
    try:
        if any(x in str(value) for x in [".", "e", "E"]):
            return float(value)
        return int(value)
    except ValueError:
        return value


def write_json(path, payload):
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--settings", required=True)
    parser.add_argument("--prefilter-qc", required=True)
    parser.add_argument("--filter-impact", required=True)
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    settings = read_key_value_tsv(args.settings, "setting", "value")
    impact = read_key_value_tsv(args.filter_impact, "metric", "value")
    prefilter_rows = read_rows(args.prefilter_qc)

    impact_data = {
        args.patient_id: {
            "prefilter_nonref_calls": coerce(impact.get("prefilter_nonref_calls", "")),
            "final_retained_calls": coerce(impact.get("final_retained_calls", "")),
            "removed_calls": coerce(impact.get("removed_calls", "")),
        }
    }
    write_json(
        outdir / f"{args.patient_id}.filter_impact_mqc.json",
        {
            "id": "cnpup_filter_impact",
            "section_name": "CN-PUP filter impact",
            "description": "Number of MonoVar non-reference calls before filtering and retained after germline, quality, population AF, and COSMIC filtering.",
            "plot_type": "bargraph",
            "pconfig": {
                "id": "cnpup_filter_impact_plot",
                "title": "Filter impact",
                "ylab": "Variant calls",
            },
            "data": impact_data,
        },
    )

    prefilter_data = {}
    for row in prefilter_rows:
        sample = row.get("cell_id", "unknown")
        prefilter_data[sample] = {
            "total_sites": coerce(row.get("total_sites", "")),
            "nonref_sites": coerce(row.get("nonref_sites", "")),
            "nonref_dp_median": coerce(row.get("nonref_dp_median", "")),
            "nonref_alt_median": coerce(row.get("nonref_alt_median", "")),
            "nonref_vaf_median": coerce(row.get("nonref_vaf_median", "")),
            "final_variant_count": coerce(row.get("final_variant_count", "")),
        }
    write_json(
        outdir / f"{args.patient_id}.prefilter_depth_alt_qc_mqc.json",
        {
            "id": "cnpup_prefilter_depth_alt_qc",
            "section_name": "CN-PUP pre-filter depth / alt-read QC",
            "description": "Per-cell MonoVar candidate-site depth, alternative-read and VAF summaries before final filtering.",
            "plot_type": "table",
            "pconfig": {
                "id": "cnpup_prefilter_depth_alt_qc_table",
                "title": "Pre-filter depth and alt-read QC",
            },
            "data": prefilter_data,
        },
    )

    write_json(
        outdir / f"{args.patient_id}.filter_settings_mqc.json",
        {
            "id": "cnpup_filter_settings",
            "section_name": "CN-PUP filter settings",
            "description": "Thresholds and modes used for this run.",
            "plot_type": "table",
            "pconfig": {
                "id": "cnpup_filter_settings_table",
                "title": "Filter settings",
            },
            "data": {args.patient_id: {k: coerce(v) for k, v in settings.items()}},
        },
    )


if __name__ == "__main__":
    main()
