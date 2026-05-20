#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".gz":
        return pd.read_csv(path, sep="\t", dtype=str, compression="gzip")
    return pd.read_csv(path, sep="\t", dtype=str)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def num_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def read_cell_metadata(path: str | Path, patient_id: str) -> pd.DataFrame:
    meta = read_table(path)
    if "patient_id" in meta.columns:
        meta = meta[meta["patient_id"].astype(str) == patient_id].copy()
    expected = ["patient_id", "cell_id", "bam_path", "cell_type"]
    for col in expected:
        if col not in meta.columns:
            meta[col] = ""
    meta["cell_id"] = meta["cell_id"].astype(str)
    meta["cell_type"] = meta["cell_type"].fillna("").astype(str)
    return meta[expected].drop_duplicates("cell_id")


def is_ctc_type(x: str) -> bool:
    x = str(x).lower()
    return "ctc" in x or x in {"single_cell", "single-cell", "cell"}


def is_leukocyte_type(x: str) -> bool:
    x = str(x).lower()
    return any(k in x for k in ["leuk", "wbc", "white"])


def normalize_long(long_df: pd.DataFrame) -> pd.DataFrame:
    rename = {
        "CHROM": "chrom",
        "POS": "pos",
        "REF": "ref",
        "ALT": "alt",
        "DP": "total_depth",
        "VAF": "vaf",
        "FILTER": "filter",
    }
    long_df = long_df.rename(columns={k: v for k, v in rename.items() if k in long_df.columns})
    for col in [
        "sample_id",
        "caller",
        "var_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene_symbol",
        "effect",
        "impact",
        "total_depth",
        "altread",
        "vaf",
    ]:
        if col not in long_df.columns:
            long_df[col] = ""
    long_df["sample_id"] = long_df["sample_id"].astype(str)
    long_df["caller"] = long_df["caller"].astype(str)
    long_df["var_id"] = long_df["var_id"].astype(str)
    long_df["total_depth_num"] = num_series(long_df["total_depth"])
    long_df["altread_num"] = num_series(long_df["altread"])
    long_df["vaf_num"] = num_series(long_df["vaf"])
    return long_df


def read_binary_matrix(path: Path, callers: set[str] | None = None) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    usecols = None
    if callers:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
        usecols = ["var_id"] + [
            col for col in header[1:]
            if "__" in col and col.rsplit("__", 1)[1] in callers
        ]
    mat = pd.read_csv(path, sep="\t", dtype=str, usecols=usecols)
    if "var_id" not in mat.columns:
        return pd.DataFrame()
    mat = mat.set_index("var_id")
    return mat.apply(pd.to_numeric, errors="coerce")


def parse_samtools_stats(path: Path) -> dict[str, float | str]:
    metrics = {}
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("SN\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            key = parts[1].rstrip(":")
            value = parts[2]
            try:
                metrics[key] = float(value)
            except ValueError:
                metrics[key] = value
    raw = metrics.get("raw total sequences")
    mapped = metrics.get("reads mapped")
    if isinstance(raw, float) and raw > 0 and isinstance(mapped, float):
        metrics["mapped_fraction"] = mapped / raw
    return metrics


def parse_mosdepth_summary(path: Path) -> dict[str, float]:
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return {}
    if df.empty:
        return {}
    name_col = "chrom" if "chrom" in df.columns else df.columns[0]
    row = df[df[name_col].astype(str).str.lower() == "total"]
    if row.empty:
        row = df.tail(1)
    out = {}
    for col in ["mean", "min", "max"]:
        if col in row.columns:
            val = pd.to_numeric(row[col], errors="coerce").iloc[0]
            if pd.notna(val):
                out[f"mosdepth_{col}"] = float(val)
    return out


def read_alignment_qc(path: str | Path | None) -> dict[str, dict[str, float | str]]:
    if not path:
        return {}
    root = Path(path)
    if not root.exists():
        return {}
    qc = defaultdict(dict)
    for stats in root.rglob("*.samtools.stats.out"):
        cell_id = stats.name.replace(".samtools.stats.out", "")
        qc[cell_id].update(parse_samtools_stats(stats))
    for summary in root.rglob("*.mosdepth.summary.txt"):
        cell_id = summary.name.replace(".mosdepth.summary.txt", "")
        qc[cell_id].update(parse_mosdepth_summary(summary))
    return dict(qc)


def sample_from_matrix_col(col: str) -> str:
    return col.rsplit("__", 1)[0] if "__" in col else col


def classify_samples(long_df: pd.DataFrame, meta: pd.DataFrame) -> tuple[set[str], set[str], set[str], set[str]]:
    ctc_from_meta = set(meta.loc[meta["cell_type"].map(is_ctc_type), "cell_id"].astype(str))
    leuk_from_meta = set(meta.loc[meta["cell_type"].map(is_leukocyte_type), "cell_id"].astype(str))
    all_samples = set(long_df["sample_id"].dropna().astype(str))
    tissue_samples = {s for s in all_samples if s not in ctc_from_meta and s not in leuk_from_meta}
    metastasis_samples = {s for s in tissue_samples if any(k in s.lower() for k in ["met", "meta", "metast"])}
    return ctc_from_meta, tissue_samples, metastasis_samples, leuk_from_meta


def ctc_support_sets(binary: pd.DataFrame, ctc_ids: set[str]) -> dict[str, set[str]]:
    support = defaultdict(set)
    if binary.empty:
        return support
    for col in binary.columns:
        sample = sample_from_matrix_col(col)
        if sample not in ctc_ids:
            continue
        present_vars = binary.index[binary[col] == 1]
        for var_id in present_vars:
            support[str(var_id)].add(sample)
    return support


def build_ctc_qc(long_df: pd.DataFrame, meta: pd.DataFrame, support: dict[str, set[str]], ctc_ids: set[str], alignment_qc: dict[str, dict[str, float | str]]) -> pd.DataFrame:
    columns = [
        "cell_id",
        "cell_type",
        "bam_path",
        "qc_status",
        "qc_warnings",
        "raw_reads",
        "mapped_reads",
        "mapped_fraction",
        "mean_bam_coverage",
        "variant_count",
        "recurrent_variant_count",
        "private_variant_count",
        "median_variant_site_depth",
        "median_variant_site_vaf",
        "callers",
    ]
    rows = []
    recurrent = {v for v, cells in support.items() if len(cells) >= 2}
    for _, m in meta.iterrows():
        cell_id = str(m["cell_id"])
        if cell_id not in ctc_ids:
            continue
        sub = long_df[long_df["sample_id"] == cell_id].copy()
        var_set = set(sub["var_id"])
        warnings = []
        median_depth = sub["total_depth_num"].median() if not sub.empty else np.nan
        median_vaf = sub["vaf_num"].median() if not sub.empty else np.nan
        if sub.empty:
            warnings.append("No retained variant rows")
        if pd.isna(median_depth):
            warnings.append("No variant-site depth available")
        elif median_depth < 10:
            warnings.append("Low median depth at retained variant sites")
        aqc = alignment_qc.get(cell_id, {})
        mean_bam_coverage = aqc.get("mosdepth_mean", "")
        mapped_fraction = aqc.get("mapped_fraction", "")
        if mean_bam_coverage != "":
            try:
                if float(mean_bam_coverage) < 10:
                    warnings.append("Low mean BAM coverage")
            except (TypeError, ValueError):
                pass
        if mapped_fraction != "":
            try:
                if float(mapped_fraction) < 0.8:
                    warnings.append("Low mapped-read fraction")
            except (TypeError, ValueError):
                pass
        status = "PASS"
        if warnings:
            status = "WARN"
        rows.append(
            {
                "cell_id": cell_id,
                "cell_type": m.get("cell_type", ""),
                "bam_path": m.get("bam_path", ""),
                "qc_status": status,
                "qc_warnings": "; ".join(warnings),
                "raw_reads": aqc.get("raw total sequences", ""),
                "mapped_reads": aqc.get("reads mapped", ""),
                "mapped_fraction": "" if mapped_fraction == "" else round(float(mapped_fraction), 5),
                "mean_bam_coverage": "" if mean_bam_coverage == "" else round(float(mean_bam_coverage), 3),
                "variant_count": len(var_set),
                "recurrent_variant_count": len(var_set & recurrent),
                "private_variant_count": sum(1 for v in var_set if len(support.get(v, set())) == 1),
                "median_variant_site_depth": "" if pd.isna(median_depth) else round(float(median_depth), 3),
                "median_variant_site_vaf": "" if pd.isna(median_vaf) else round(float(median_vaf), 5),
                "callers": ", ".join(sorted(set(sub["caller"].dropna().astype(str)))) if not sub.empty else "",
            }
        )
    return pd.DataFrame(rows, columns=columns)


def build_sharing_spectrum(support: dict[str, set[str]]) -> pd.DataFrame:
    counts = defaultdict(int)
    for cells in support.values():
        counts[len(cells)] += 1
    return pd.DataFrame(
        [{"ctc_support_count": k, "variant_count": v} for k, v in sorted(counts.items())],
        columns=["ctc_support_count", "variant_count"],
    )


def build_top_variants(long_df: pd.DataFrame, support: dict[str, set[str]], limit: int) -> pd.DataFrame:
    columns = [
        "support_ctc_count",
        "ctc_samples",
        "var_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene_symbol",
        "effect",
        "impact",
    ]
    rows = []
    by_var = long_df.drop_duplicates("var_id").set_index("var_id", drop=False)
    for var_id, cells in support.items():
        if var_id not in by_var.index:
            continue
        r = by_var.loc[var_id]
        rows.append(
            {
                "support_ctc_count": len(cells),
                "ctc_samples": ", ".join(sorted(cells)),
                "var_id": var_id,
                "chrom": r.get("chrom", ""),
                "pos": r.get("pos", ""),
                "ref": r.get("ref", ""),
                "alt": r.get("alt", ""),
                "gene_symbol": r.get("gene_symbol", ""),
                "effect": r.get("effect", ""),
                "impact": r.get("impact", ""),
            }
        )
    out = pd.DataFrame(rows, columns=columns)
    if out.empty:
        return out
    return out.sort_values(["support_ctc_count", "var_id"], ascending=[False, True]).head(limit)


def build_tissue_overlap(long_df: pd.DataFrame, support: dict[str, set[str]], ctc_ids: set[str], tissue_samples: set[str]) -> pd.DataFrame:
    if not tissue_samples:
        return pd.DataFrame(columns=["cell_id", "tissue_variant_count", "ctc_variant_count", "shared_count", "jaccard", "ctc_fraction_shared_with_tissue"])
    tissue_vars = set(long_df.loc[long_df["sample_id"].isin(tissue_samples), "var_id"])
    rows = []
    for cell_id in sorted(ctc_ids):
        ctc_vars = {v for v, cells in support.items() if cell_id in cells}
        inter = ctc_vars & tissue_vars
        union = ctc_vars | tissue_vars
        rows.append(
            {
                "cell_id": cell_id,
                "tissue_variant_count": len(tissue_vars),
                "ctc_variant_count": len(ctc_vars),
                "shared_count": len(inter),
                "jaccard": round(len(inter) / len(union), 6) if union else 0,
                "ctc_fraction_shared_with_tissue": round(len(inter) / len(ctc_vars), 6) if ctc_vars else 0,
            }
        )
    return pd.DataFrame(rows)


def build_rarefaction(support: dict[str, set[str]], ctc_ids: set[str], iterations: int) -> pd.DataFrame:
    ctc_list = sorted(ctc_ids)
    if not ctc_list:
        return pd.DataFrame(columns=["n_ctc", "mean_variants_observed", "mean_recurrent_variants_observed"])
    rng = random.Random(1)
    rows = []
    for n in range(1, len(ctc_list) + 1):
        totals = []
        recurrent_totals = []
        reps = 1 if n == len(ctc_list) else iterations
        for _ in range(reps):
            selected = set(ctc_list if n == len(ctc_list) else rng.sample(ctc_list, n))
            observed = {v for v, cells in support.items() if cells & selected}
            recurrent = {v for v, cells in support.items() if len(cells & selected) >= 2}
            totals.append(len(observed))
            recurrent_totals.append(len(recurrent))
        rows.append(
            {
                "n_ctc": n,
                "mean_variants_observed": round(float(np.mean(totals)), 3),
                "mean_recurrent_variants_observed": round(float(np.mean(recurrent_totals)), 3),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build compact data tables for the CN-PUP final patient report.")
    ap.add_argument("--patient-id", required=True)
    ap.add_argument("--long-table", required=True)
    ap.add_argument("--matrices-dir", required=True)
    ap.add_argument("--cell-metadata", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--alignment-qc-dir", default="")
    ap.add_argument("--callers", default="", help="Comma-separated callers to include, for example monovar.")
    ap.add_argument("--top-variant-limit", type=int, default=200)
    ap.add_argument("--rarefaction-iterations", type=int, default=200)
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    long_df = normalize_long(read_table(args.long_table))
    callers_filter = {x.strip() for x in args.callers.split(",") if x.strip()}
    if callers_filter:
        long_df = long_df[long_df["caller"].isin(callers_filter)].copy()
    meta = read_cell_metadata(args.cell_metadata, args.patient_id)
    ctc_ids, tissue_samples, metastasis_samples, leukocyte_samples = classify_samples(long_df, meta)
    binary = read_binary_matrix(Path(args.matrices_dir) / "mutation_binary.tsv", callers_filter or None)
    support = ctc_support_sets(binary, ctc_ids)
    alignment_qc = read_alignment_qc(args.alignment_qc_dir)

    warnings = []
    if binary.empty:
        warnings.append({"level": "WARN", "message": "mutation_binary.tsv was missing or empty; CTC support metrics are unavailable."})
    if not tissue_samples:
        warnings.append({"level": "INFO", "message": "No non-CTC tissue samples were detected in the comparison table; tissue representativity metrics are not available."})
    if alignment_qc:
        warnings.append({"level": "INFO", "message": "BAM-level coverage was read from samtools/mosdepth outputs; variant-site depth is also shown for called SNVs."})
        coverage_source = "bam_mosdepth_and_variant_site_depth_from_call_tables"
    else:
        warnings.append({"level": "INFO", "message": "No samtools/mosdepth files were provided to the final report; coverage metrics are variant-site depths from caller tables only."})
        coverage_source = "variant_site_depth_from_call_tables"

    unique_variants = long_df["var_id"].nunique()
    recurrent_ge2 = sum(1 for cells in support.values() if len(cells) >= 2)
    recurrent_ge3 = sum(1 for cells in support.values() if len(cells) >= 3)
    max_support = max([len(cells) for cells in support.values()] or [0])
    summary = pd.DataFrame(
        [
            {"metric": "patient_id", "value": args.patient_id},
            {"metric": "ctc_count", "value": len(ctc_ids)},
            {"metric": "tissue_sample_count", "value": len(tissue_samples)},
            {"metric": "metastasis_sample_count", "value": len(metastasis_samples)},
            {"metric": "leukocyte_sample_count", "value": len(leukocyte_samples)},
            {"metric": "caller_count", "value": long_df["caller"].nunique()},
            {"metric": "callers", "value": ", ".join(sorted(long_df["caller"].dropna().unique()))},
            {"metric": "retained_variant_rows", "value": len(long_df)},
            {"metric": "unique_variants_all_samples", "value": unique_variants},
            {"metric": "ctc_variants_with_support_ge2", "value": recurrent_ge2},
            {"metric": "ctc_variants_with_support_ge3", "value": recurrent_ge3},
            {"metric": "max_ctc_support_for_one_variant", "value": max_support},
            {"metric": "coverage_source", "value": coverage_source},
        ]
    )

    ctc_qc = build_ctc_qc(long_df, meta, support, ctc_ids, alignment_qc)
    write_tsv(summary, out_dir / "patient_summary.tsv")
    write_tsv(ctc_qc, out_dir / "ctc_qc_summary.tsv")
    write_tsv(build_sharing_spectrum(support), out_dir / "sharing_spectrum.tsv")
    write_tsv(build_top_variants(long_df, support, args.top_variant_limit), out_dir / "top_recurrent_variants.tsv")
    write_tsv(build_tissue_overlap(long_df, support, ctc_ids, tissue_samples), out_dir / "ctc_tissue_overlap.tsv")
    write_tsv(build_rarefaction(support, ctc_ids, args.rarefaction_iterations), out_dir / "rarefaction.tsv")
    write_tsv(pd.DataFrame(warnings), out_dir / "warnings.tsv")


if __name__ == "__main__":
    main()
