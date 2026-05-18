#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else path.open("r", encoding="utf-8", errors="replace")


def split_info(value: str) -> dict[str, str | bool]:
    out: dict[str, str | bool] = {}
    if not value or value == ".":
        return out
    for item in value.split(";"):
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            out[key] = val
        else:
            out[item] = True
    return out


def parse_sample(format_field: str, sample_field: str) -> dict[str, str]:
    if not format_field or not sample_field:
        return {}
    keys = format_field.split(":")
    values = sample_field.split(":")
    return {key: values[i] if i < len(values) else "" for i, key in enumerate(keys)}


def ann_categories(info: dict[str, str | bool], alt: str) -> tuple[bool, bool, bool]:
    raw_ann = str(info.get("ANN", ""))
    raw_csq = str(info.get("CSQ", ""))
    annotation = raw_ann or raw_csq
    if not annotation:
        return False, False, False
    entries = [entry for entry in annotation.split(",") if entry]
    matched = [entry for entry in entries if entry.split("|", 1)[0] == alt] or entries
    text = ",".join(matched).lower()
    non_synonymous = any(token in text for token in [
        "missense", "stop_gained", "stop_lost", "start_lost", "frameshift",
        "inframe", "splice", "protein_altering", "transcript_ablation",
    ])
    predictive_impact = any(token in text for token in ["high", "moderate", "deleterious", "damaging"])
    listed_db = any(key in info for key in ["DB", "COSMIC", "COSMIC_CNT", "CNT"]) or "cosmic" in text or "dbsnp" in text
    return non_synonymous, listed_db, predictive_impact


def read_vcf_variants(path: Path, annotated: bool = False) -> dict[str, dict[str, bool]]:
    variants: dict[str, dict[str, bool]] = {}
    if not path or not path.exists():
        return variants
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _vid, ref, alts, _qual, _filt, info_s = fields[:8]
            info = split_info(info_s)
            for alt in alts.split(","):
                if len(ref) != 1 or len(alt) != 1 or alt in {".", "*"}:
                    continue
                var_id = f"{chrom}:{pos}_{ref}/{alt}"
                non_syn, listed, predictive = ann_categories(info, alt) if annotated else (False, False, False)
                variants[var_id] = {
                    "non_synonymous": non_syn,
                    "listed_db": listed,
                    "predictive_impact": predictive,
                }
    return variants


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = []
        for row in csv.DictReader(handle, delimiter="\t"):
            sample_id = row.get("sample_id", "").strip()
            caller = row.get("caller", "").strip().lower()
            raw_vcf = row.get("raw_vcf", "").strip()
            if sample_id and caller and raw_vcf and not sample_id.startswith("#"):
                rows.append({"sample_id": sample_id, "caller": caller, "raw_vcf": raw_vcf})
        return rows


def first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def find_stage_vcf(project: Path, stage: str, caller: str, sample_id: str) -> Path | None:
    root = project / stage / caller
    if not root.exists():
        return None
    patterns = [
        f"{sample_id}.{caller}.comparable.vcf.gz",
        f"{sample_id}.{caller}.filtered.snvs.norm.vcf.gz",
        f"{sample_id}.{caller}.filtered.snvs.vcf.gz",
        f"{sample_id}.{caller}.filtered.*.snvs.norm.vcf.gz",
        f"{sample_id}.{caller}.filtered.*.snvs.vcf.gz",
    ]
    for pattern in patterns:
        matches = sorted(root.glob(pattern))
        if matches:
            return matches[0]
    return None


def read_tissue_sets(path: Path | None) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = defaultdict(set)
    if not path or not path.exists():
        return sets
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_type = row.get("sample_type", "").strip().lower()
            vcf_path = row.get("vcf_path", "").strip()
            if not sample_type or not vcf_path:
                continue
            variants = read_vcf_variants(Path(vcf_path))
            sets[sample_type].update(variants)
    return sets


def safe_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", value)
    return cleaned.strip("_") or "sample"


def label(text: str, value: int | str) -> str:
    return f"{text}<br/>{value}"


def write_mermaid(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    title = f"{row['sample_id']} / {row['caller']}"
    metastasis_value = row["shared_with_metastasis"] if row["shared_with_metastasis"] != "" else "NA"
    lines = [
        "flowchart TD",
        f"  T[\"{title}\"]",
        f"  A[\"{label('Raw input SNVs', row['raw_input_snvs'])}\"]",
        f"  B[\"{label('After QC filters and normalization', row['normalized_snvs'])}\"]",
        f"  C[\"{label('Seen in bulk blood exclusion', row['seen_in_bulk_blood'])}\"]",
        f"  D[\"{label('Seen in leukocyte exclusion', row['seen_in_leukocyte'])}\"]",
        f"  E[\"{label('After background subtraction', row['comparable_snvs'])}\"]",
        f"  F[\"{label('Shared with metastasis', metastasis_value)}\"]",
        f"  G[\"{label('Shared with at least one other CTC', row['shared_with_at_least_one_other_ctc'])}\"]",
        f"  H[\"{label('Shared with all CTCs for this caller', row['shared_with_all_ctcs'])}\"]",
        f"  I[\"{label('Non-synonymous', row['non_synonymous'])}\"]",
        f"  J[\"{label('Listed in dbSNP/COSMIC annotation', row['listed_db_or_cosmic'])}\"]",
        f"  K[\"{label('Predicted protein impact', row['predictive_protein_impact'])}\"]",
        "  T --> A",
        "  A --> B",
        "  B --> C",
        "  B --> D",
        "  B --> E",
        "  E --> F",
        "  E --> G",
        "  E --> H",
        "  E --> I",
        "  I --> J",
        "  J --> K",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate per-sample variant attrition flowchart counts and Mermaid diagrams.")
    parser.add_argument("--project", required=True, help="Comparison project directory containing processed_calls, normalized_calls, comparable_calls, and bulk_exclusion.")
    parser.add_argument("--samples-tsv", required=True, help="samples.tsv used by the comparison branch.")
    parser.add_argument("--output-dir", required=True, help="Output directory for summary TSV and Mermaid flowcharts.")
    parser.add_argument("--tissue-vcfs", help="Optional TSV with sample_id, sample_type, vcf_path columns for primary/metastasis overlap.")
    args = parser.parse_args()

    project = Path(args.project)
    output_dir = Path(args.output_dir)
    chart_dir = output_dir / "mermaid"
    output_dir.mkdir(parents=True, exist_ok=True)
    chart_dir.mkdir(parents=True, exist_ok=True)

    samples = read_samples(Path(args.samples_tsv))
    tissue_sets = read_tissue_sets(Path(args.tissue_vcfs) if args.tissue_vcfs else None)
    metastasis_set = set().union(*(vars for key, vars in tissue_sets.items() if "metast" in key)) if tissue_sets else set()

    bulk_vcf = first_existing(sorted((project / "bulk_exclusion").glob("*.bulk_blood.norm.vcf.gz")))
    leukocyte_vcf = first_existing(sorted((project / "bulk_exclusion").glob("*.leukocyte.norm.vcf.gz")))
    bulk_set = set(read_vcf_variants(bulk_vcf).keys()) if bulk_vcf else set()
    leukocyte_set = set(read_vcf_variants(leukocyte_vcf).keys()) if leukocyte_vcf else set()

    comparable_sets: dict[tuple[str, str], set[str]] = {}
    comparable_meta: dict[tuple[str, str], dict[str, dict[str, bool]]] = {}
    raw_sets: dict[tuple[str, str], set[str]] = {}
    norm_sets: dict[tuple[str, str], set[str]] = {}

    for sample in samples:
        sample_id, caller = sample["sample_id"], sample["caller"]
        key = (sample_id, caller)
        raw_sets[key] = set(read_vcf_variants(Path(sample["raw_vcf"])).keys())
        norm_vcf = find_stage_vcf(project, "normalized_calls", caller, sample_id)
        comp_vcf = find_stage_vcf(project, "comparable_calls", caller, sample_id)
        norm_sets[key] = set(read_vcf_variants(norm_vcf).keys()) if norm_vcf else set()
        comp_variants = read_vcf_variants(comp_vcf, annotated=True) if comp_vcf else {}
        comparable_sets[key] = set(comp_variants.keys())
        comparable_meta[key] = comp_variants

    rows: list[dict[str, str]] = []
    fieldnames = [
        "sample_id", "caller", "raw_input_snvs", "normalized_snvs",
        "seen_in_bulk_blood", "seen_in_leukocyte", "removed_by_any_background",
        "comparable_snvs", "shared_with_metastasis",
        "shared_with_at_least_one_other_ctc", "shared_with_all_ctcs",
        "non_synonymous", "listed_db_or_cosmic", "predictive_protein_impact",
        "mermaid_path",
    ]

    callers = sorted({caller for _sample, caller in comparable_sets})
    all_by_caller: dict[str, set[str]] = {}
    for caller in callers:
        caller_sets = [vars for (sample_id, c), vars in comparable_sets.items() if c == caller and sample_id.lower().startswith("ctc")]
        all_by_caller[caller] = set.intersection(*caller_sets) if caller_sets else set()

    for sample in samples:
        sample_id, caller = sample["sample_id"], sample["caller"]
        key = (sample_id, caller)
        raw = raw_sets.get(key, set())
        norm = norm_sets.get(key, set())
        comp = comparable_sets.get(key, set())
        other_ctc_sets = [
            vars for (other_sample, other_caller), vars in comparable_sets.items()
            if other_caller == caller and other_sample != sample_id and other_sample.lower().startswith("ctc")
        ]
        shared_other = comp & set().union(*other_ctc_sets) if other_ctc_sets else set()
        shared_all = comp & all_by_caller.get(caller, set())
        meta = comparable_meta.get(key, {})
        non_syn = {var for var, ann in meta.items() if ann.get("non_synonymous")}
        listed = {var for var, ann in meta.items() if ann.get("listed_db")}
        predictive = {var for var, ann in meta.items() if ann.get("predictive_impact")}

        mermaid_rel = Path("mermaid") / f"{safe_id(sample_id)}.{safe_id(caller)}.variant_flow.mmd"
        row = {
            "sample_id": sample_id,
            "caller": caller,
            "raw_input_snvs": str(len(raw)),
            "normalized_snvs": str(len(norm)),
            "seen_in_bulk_blood": str(len(norm & bulk_set)) if bulk_set else "0",
            "seen_in_leukocyte": str(len(norm & leukocyte_set)) if leukocyte_set else "0",
            "removed_by_any_background": str(len(norm - comp)) if norm else "0",
            "comparable_snvs": str(len(comp)),
            "shared_with_metastasis": str(len(comp & metastasis_set)) if metastasis_set else "",
            "shared_with_at_least_one_other_ctc": str(len(shared_other)),
            "shared_with_all_ctcs": str(len(shared_all)),
            "non_synonymous": str(len(non_syn)),
            "listed_db_or_cosmic": str(len(listed)),
            "predictive_protein_impact": str(len(predictive)),
            "mermaid_path": str(mermaid_rel),
        }
        write_mermaid(output_dir / mermaid_rel, row)
        rows.append(row)

    summary_path = output_dir / "variant_flow_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    index_path = output_dir / "variant_flowcharts.md"
    with index_path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# Variant Flowcharts\n\n")
        out.write("Each Mermaid diagram shows variant attrition and sharing for one sample/caller pair.\n\n")
        for row in rows:
            out.write(f"## {row['sample_id']} / {row['caller']}\n\n")
            out.write("```{mermaid}\n")
            out.write((output_dir / row["mermaid_path"]).read_text(encoding="utf-8"))
            out.write("```\n\n")

    print(f"Wrote {summary_path}")
    print(f"Wrote {len(rows)} Mermaid flowcharts under {chart_dir}")


if __name__ == "__main__":
    main()
