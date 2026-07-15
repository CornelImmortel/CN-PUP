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


def write_bargraph(path, section_id, section_name, description, title, ylab, data):
    write_json(
        path,
        {
            "id": section_id,
            "section_name": section_name,
            "description": description,
            "plot_type": "bargraph",
            "pconfig": {
                "id": f"{section_id}_plot",
                "title": title,
                "ylab": ylab,
            },
            "data": data,
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--settings", required=True)
    parser.add_argument("--prefilter-qc", required=True)
    parser.add_argument("--filter-impact", required=True)
    parser.add_argument("--breadth-summary", default=None, help="Optional reports/<patient>.breadth_summary.tsv from SUMMARIZE_BREADTH_QC")
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
            "final_retained_calls_high_confidence": coerce(impact.get("final_retained_calls_high_confidence", "")),
            "flagged_low_depth_alt_calls": coerce(impact.get("flagged_low_depth_alt_calls", "")),
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
            "nonref_het_sites": coerce(row.get("nonref_het_sites", "")),
            "nonref_hom_alt_sites": coerce(row.get("nonref_hom_alt_sites", "")),
            "nonref_other_gt_sites": coerce(row.get("nonref_other_gt_sites", "")),
            "nonref_pass_sites": coerce(row.get("nonref_pass_sites", "")),
            "nonref_dot_filter_sites": coerce(row.get("nonref_dot_filter_sites", "")),
            "nonref_filtered_sites": coerce(row.get("nonref_filtered_sites", "")),
            "final_retention_fraction": coerce(row.get("final_retention_fraction", "")),
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

    final_counts = {
        sample: {"final_variant_count": values.get("final_variant_count")}
        for sample, values in prefilter_data.items()
    }
    write_bargraph(
        outdir / f"{args.patient_id}.final_variant_counts_mqc.json",
        "cnpup_final_variant_counts",
        "CN-PUP final variant counts",
        "Retained variants per cell after germline/background, quality, population AF and COSMIC filtering.",
        "Final retained variants per cell",
        "Variants",
        final_counts,
    )

    nonref_counts = {
        sample: {"prefilter_nonref_sites": values.get("nonref_sites")}
        for sample, values in prefilter_data.items()
    }
    write_bargraph(
        outdir / f"{args.patient_id}.prefilter_nonref_counts_mqc.json",
        "cnpup_prefilter_nonref_counts",
        "CN-PUP pre-filter non-reference calls",
        "MonoVar non-reference calls per cell before downstream filtering.",
        "Pre-filter non-reference calls per cell",
        "Non-reference calls",
        nonref_counts,
    )

    depth_medians = {
        sample: {"median_depth": values.get("nonref_dp_median")}
        for sample, values in prefilter_data.items()
    }
    write_bargraph(
        outdir / f"{args.patient_id}.prefilter_depth_median_mqc.json",
        "cnpup_prefilter_depth_median",
        "CN-PUP median depth before filtering",
        "Median depth among MonoVar non-reference candidate sites before downstream filtering.",
        "Median depth at non-reference sites",
        "Depth",
        depth_medians,
    )

    alt_medians = {
        sample: {"median_alt_reads": values.get("nonref_alt_median")}
        for sample, values in prefilter_data.items()
    }
    write_bargraph(
        outdir / f"{args.patient_id}.prefilter_alt_median_mqc.json",
        "cnpup_prefilter_alt_median",
        "CN-PUP median alt reads before filtering",
        "Median alternative-read count among MonoVar non-reference candidate sites before downstream filtering.",
        "Median alternative reads at non-reference sites",
        "Alternative reads",
        alt_medians,
    )

    vaf_medians = {
        sample: {"median_vaf": values.get("nonref_vaf_median")}
        for sample, values in prefilter_data.items()
    }
    write_bargraph(
        outdir / f"{args.patient_id}.prefilter_vaf_median_mqc.json",
        "cnpup_prefilter_vaf_median",
        "CN-PUP median VAF before filtering",
        "Median VAF among MonoVar non-reference candidate sites before downstream filtering.",
        "Median VAF at non-reference sites",
        "VAF",
        vaf_medians,
    )

    # MonoVar QUAL = -log10(P(no variant present anywhere in the cohort at
    # this site)) -- site-level, not per-cell, and NOT true Phred scale (no
    # x10 factor). Meaningful range is roughly 0-15.
    qual_medians = {
        sample: {"median_qual": values.get("nonref_qual_median")}
        for sample, values in prefilter_data.items()
    }
    write_bargraph(
        outdir / f"{args.patient_id}.prefilter_qual_median_mqc.json",
        "cnpup_prefilter_qual_median",
        "CN-PUP median QUAL before filtering",
        "Median MonoVar QUAL among non-reference candidate sites before downstream filtering "
        "(-log10 of the site's joint probability of having no variant; not Phred-scaled).",
        "Median QUAL at non-reference sites",
        "QUAL",
        qual_medians,
    )

    per_cell_drop = {}
    for sample, values in prefilter_data.items():
        nonref = values.get("nonref_sites") or 0
        final = values.get("final_variant_count") or 0
        per_cell_drop[sample] = {
            "removed_calls": max(nonref - final, 0),
            "final_retained_calls": final,
        }
    write_bargraph(
        outdir / f"{args.patient_id}.per_cell_filter_drop_mqc.json",
        "cnpup_per_cell_filter_drop",
        "CN-PUP per-cell filter drop",
        "Per-cell calls retained and removed by downstream filters after MonoVar non-reference candidate calling.",
        "Per-cell filter impact",
        "Variant calls",
        per_cell_drop,
    )

    genotype_balance = {
        sample: {
            "heterozygous_like": values.get("nonref_het_sites"),
            "homozygous_alt_like": values.get("nonref_hom_alt_sites"),
            "other_nonref_gt": values.get("nonref_other_gt_sites"),
        }
        for sample, values in prefilter_data.items()
    }
    write_bargraph(
        outdir / f"{args.patient_id}.single_cell_gt_balance_mqc.json",
        "cnpup_single_cell_gt_balance",
        "CN-PUP single-cell genotype balance",
        "MonoVar non-reference genotype classes before downstream filtering. Strong shifts can indicate allele dropout, amplification bias, or caller-specific behavior.",
        "Pre-filter genotype classes",
        "Non-reference calls",
        genotype_balance,
    )

    filter_status = {
        sample: {
            "PASS": values.get("nonref_pass_sites"),
            "dot_unfiltered": values.get("nonref_dot_filter_sites"),
            "filtered": values.get("nonref_filtered_sites"),
        }
        for sample, values in prefilter_data.items()
    }
    write_bargraph(
        outdir / f"{args.patient_id}.single_cell_filter_status_mqc.json",
        "cnpup_single_cell_filter_status",
        "CN-PUP pre-filter VCF FILTER status",
        "FILTER values among MonoVar non-reference candidate calls before CN-PUP downstream filtering.",
        "Pre-filter VCF FILTER status",
        "Non-reference calls",
        filter_status,
    )

    retention_fraction = {
        sample: {"retention_fraction": values.get("final_retention_fraction")}
        for sample, values in prefilter_data.items()
    }
    write_bargraph(
        outdir / f"{args.patient_id}.single_cell_retention_fraction_mqc.json",
        "cnpup_single_cell_retention_fraction",
        "CN-PUP final retention fraction",
        "Fraction of MonoVar non-reference candidate calls retained after background, quality, population AF and COSMIC filtering.",
        "Final retention fraction",
        "Fraction retained",
        retention_fraction,
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

    if args.breadth_summary and Path(args.breadth_summary).is_file():
        breadth_rows = read_rows(args.breadth_summary)
        breadth_thresholds = sorted(
            {key[len("breadth_"):] for row in breadth_rows for key in row if key.startswith("breadth_")}
        )
        breadth_table = {
            row.get("cell_id", "unknown"): {
                "coverage_scope": row.get("coverage_scope", ""),
                **{f"breadth_{t}": coerce(row.get(f"breadth_{t}", "")) for t in breadth_thresholds},
            }
            for row in breadth_rows
        }
        write_json(
            outdir / f"{args.patient_id}.breadth_summary_mqc.json",
            {
                "id": "cnpup_breadth_summary",
                "section_name": "CN-PUP breadth of coverage",
                "description": "Fraction of the target panel covered at or above each depth threshold, per cell (from mosdepth region.dist.txt).",
                "plot_type": "table",
                "pconfig": {
                    "id": "cnpup_breadth_summary_table",
                    "title": "Breadth of coverage",
                },
                "data": breadth_table,
            },
        )
        if breadth_thresholds:
            breadth_bargraph = {
                cell_id: {f"breadth_{t}": values.get(f"breadth_{t}") for t in breadth_thresholds}
                for cell_id, values in breadth_table.items()
            }
            write_bargraph(
                outdir / f"{args.patient_id}.breadth_of_coverage_mqc.json",
                "cnpup_breadth_of_coverage",
                "CN-PUP breadth of coverage (bar)",
                "Fraction of the target panel covered at or above each depth threshold, per cell.",
                "Breadth of coverage",
                "Fraction of panel covered",
                breadth_bargraph,
            )


if __name__ == "__main__":
    main()
