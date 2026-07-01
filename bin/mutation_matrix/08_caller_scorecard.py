#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, itertools, statistics
from pathlib import Path
from scipy import stats as scipy_stats


def read_matrix(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        columns = header[1:]
        rows = {}
        for row in reader:
            var_id = row[0]
            rows[var_id] = dict(zip(columns, row[1:]))
    return columns, rows


def is_ctc_sample(sample_id: str) -> bool:
    return "ctc" in sample_id.lower()


def split_column(col: str):
    sample_id, _, caller = col.rpartition("__")
    return sample_id, caller


def to_float(x):
    try:
        if x in (None, "", "."):
            return None
        return float(x)
    except ValueError:
        return None


def median(values):
    values = [v for v in values if v is not None]
    return statistics.median(values) if values else ""


def jaccard(a: set, b: set):
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def overlap_coefficient(a: set, b: set):
    smaller = min(len(a), len(b))
    return len(a & b) / smaller if smaller else 0.0


def build_caller_data(columns, binary_rows, alt_rows, depth_rows, vaf_rows):
    """For each caller: {sample_id: set(var_id called), ...} plus per-caller flat call lists."""
    by_caller_ctc_sets = {}
    by_caller_calls = {}  # caller -> list of (sample_id, var_id)
    all_ctc_samples = set()
    background_samples_by_caller = {}

    caller_columns = {}
    for col in columns:
        sample_id, caller = split_column(col)
        caller_columns.setdefault(caller, []).append((sample_id, col))
        if is_ctc_sample(sample_id):
            all_ctc_samples.add(sample_id)

    for caller, sample_cols in caller_columns.items():
        ctc_sets = {sample_id: set() for sample_id, _ in sample_cols if is_ctc_sample(sample_id)}
        background_samples_by_caller[caller] = [col for sample_id, col in sample_cols if not is_ctc_sample(sample_id)]
        by_caller_ctc_sets[caller] = ctc_sets
        by_caller_calls[caller] = []

    # Single pass over the matrix: for every "1" in a CTC column, record it against that column's caller.
    for var_id, row in binary_rows.items():
        for col, val in row.items():
            if val != "1":
                continue
            sample_id, caller = split_column(col)
            if not is_ctc_sample(sample_id):
                continue
            by_caller_ctc_sets[caller][sample_id].add(var_id)
            by_caller_calls[caller].append((sample_id, var_id))

    return by_caller_ctc_sets, by_caller_calls, background_samples_by_caller, sorted(all_ctc_samples)


def compute_scorecard_row(caller, ctc_sets, calls, background_cols, binary_rows,
                           alt_rows, depth_rows, vaf_rows, caller_columns_by_var_sample):
    ctc_samples = sorted(ctc_sets)
    n_ctcs = len(ctc_samples)
    total_calls = len(calls)
    union_variants = set()
    for s in ctc_samples:
        union_variants |= ctc_sets[s]

    support_count = {}
    for s in ctc_samples:
        for v in ctc_sets[s]:
            support_count[v] = support_count.get(v, 0) + 1

    shared_all_n_ctcs = 0
    phylo_informative = 0
    if n_ctcs >= 2:
        shared_all_n_ctcs = sum(1 for n in support_count.values() if n == n_ctcs)
        phylo_informative = sum(1 for n in support_count.values() if 2 <= n < n_ctcs)

    # A "private call" is a variant seen in exactly one CTC sample for this caller.
    private_calls = sum(1 for n in support_count.values() if n == 1)

    pairwise_jaccards = []
    pairwise_overlaps = []
    if n_ctcs >= 2:
        for s1, s2 in itertools.combinations(ctc_samples, 2):
            pairwise_jaccards.append(jaccard(ctc_sets[s1], ctc_sets[s2]))
            pairwise_overlaps.append(overlap_coefficient(ctc_sets[s1], ctc_sets[s2]))

    cross_caller_supported = 0
    for sample_id, var_id in calls:
        other_callers = caller_columns_by_var_sample.get((var_id, sample_id), set()) - {caller}
        if other_callers:
            cross_caller_supported += 1
    caller_specific = total_calls - cross_caller_supported

    background_supported = 0
    background_checked = 0
    if background_cols:
        background_checked = total_calls
        for sample_id, var_id in calls:
            row = binary_rows.get(var_id, {})
            if any(row.get(bc) == "1" for bc in background_cols):
                background_supported += 1

    depths, alts, vafs = [], [], []
    for sample_id, var_id in calls:
        col = f"{sample_id}__{caller}"
        depths.append(to_float(depth_rows.get(var_id, {}).get(col)))
        alts.append(to_float(alt_rows.get(var_id, {}).get(col)))
        vafs.append(to_float(vaf_rows.get(var_id, {}).get(col)))

    return {
        "caller": caller,
        "n_ctcs": n_ctcs,
        "total_calls": total_calls,
        "union_variants": len(union_variants),
        "shared_all_n_ctcs": shared_all_n_ctcs,
        "phylo_informative_variants": phylo_informative,
        "phylo_informative_fraction": phylo_informative / len(union_variants) if union_variants else 0.0,
        "private_calls": private_calls,
        "private_fraction": private_calls / total_calls if total_calls else 0.0,
        "mean_pairwise_jaccard": statistics.mean(pairwise_jaccards) if pairwise_jaccards else "",
        "median_pairwise_jaccard": statistics.median(pairwise_jaccards) if pairwise_jaccards else "",
        "mean_overlap_coefficient": statistics.mean(pairwise_overlaps) if pairwise_overlaps else "",
        "cross_caller_supported_calls": cross_caller_supported,
        "cross_caller_supported_fraction": cross_caller_supported / total_calls if total_calls else 0.0,
        "caller_specific_calls": caller_specific,
        "caller_specific_fraction": caller_specific / total_calls if total_calls else 0.0,
        "any_background_overlap_fraction": (background_supported / background_checked) if background_checked else "",
        "median_total_depth": median(depths),
        "median_altread": median(alts),
        "median_vaf": median(vafs),
        "_pairwise_jaccards": pairwise_jaccards,
    }


def write_scorecard(rows, out_path: Path):
    columns = [
        "caller", "n_ctcs", "total_calls", "union_variants", "shared_all_n_ctcs",
        "phylo_informative_variants", "phylo_informative_fraction", "private_calls",
        "private_fraction", "mean_pairwise_jaccard", "median_pairwise_jaccard",
        "mean_overlap_coefficient", "cross_caller_supported_calls",
        "cross_caller_supported_fraction", "caller_specific_calls",
        "caller_specific_fraction", "any_background_overlap_fraction",
        "median_total_depth", "median_altread", "median_vaf",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(columns)
        for row in rows:
            w.writerow([row.get(c, "") for c in columns])


def write_statistical_tests(rows, out_path: Path):
    tests = []
    callers = [r["caller"] for r in rows]

    jaccard_groups = [r["_pairwise_jaccards"] for r in rows if len(r["_pairwise_jaccards"]) >= 2]
    if len(jaccard_groups) >= 2:
        stat, p = scipy_stats.kruskal(*jaccard_groups)
        tests.append(("Kruskal-Wallis",
                      "Do pairwise CTC Jaccard distributions differ between callers?",
                      stat, p, f"{len(jaccard_groups)} callers with >=2 CTC pairs; descriptive/supportive."))

    if len(rows) >= 2:
        contingency = [[r["total_calls"] - r["private_calls"], r["private_calls"]] for r in rows]
        if all(sum(row) > 0 for row in contingency):
            stat, p, dof, _ = scipy_stats.chi2_contingency(contingency)
            tests.append(("Chi-square",
                          "Does private/shared call composition differ between callers?",
                          stat, p, f"dof={dof}; rows={callers}; cols=[shared,private]"))

        contingency_bg = [[r["cross_caller_supported_calls"], r["caller_specific_calls"]] for r in rows]
        if all(sum(row) > 0 for row in contingency_bg):
            stat, p, dof, _ = scipy_stats.chi2_contingency(contingency_bg)
            tests.append(("Chi-square",
                          "Does cross-caller-support composition differ between callers?",
                          stat, p, f"dof={dof}; rows={callers}; cols=[cross_caller_supported,caller_specific]"))

    for r in rows:
        if r.get("any_background_overlap_fraction") == "":
            continue
        others = [o for o in rows if o["caller"] != r["caller"] and o.get("any_background_overlap_fraction") != ""]
        if not others:
            continue
        caller_bg = round(r["any_background_overlap_fraction"] * r["total_calls"])
        caller_not_bg = r["total_calls"] - caller_bg
        other_total = sum(o["total_calls"] for o in others)
        other_bg = round(sum(o["any_background_overlap_fraction"] * o["total_calls"] for o in others))
        other_not_bg = other_total - other_bg
        table = [[caller_bg, caller_not_bg], [other_bg, other_not_bg]]
        if all(sum(row) > 0 for row in table) and min(min(table[0]), min(table[1])) >= 0:
            odds, p = scipy_stats.fisher_exact(table)
            tests.append(("Fisher exact",
                          f"Are {r['caller']} variants enriched/depleted for background support vs all other callers?",
                          odds, p, "odds ratio >1 means more background-supported than other callers."))

    with out_path.open("w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["test", "question", "statistic", "p_value", "notes"])
        for test, question, stat, p, notes in tests:
            w.writerow([test, question, f"{stat:.6g}", f"{p:.6g}", notes])


def write_summary(rows, out_path: Path):
    if not rows:
        out_path.write_text("No callers found in the input matrices.\n", encoding="utf-8")
        return
    by_shared = sorted(rows, key=lambda r: r["shared_all_n_ctcs"], reverse=True)
    by_jaccard = sorted(
        [r for r in rows if r["mean_pairwise_jaccard"] != ""],
        key=lambda r: r["mean_pairwise_jaccard"], reverse=True,
    )
    best = by_shared[0]
    sparse = [r for r in rows if r["shared_all_n_ctcs"] == 0 and r["n_ctcs"] >= 2]

    lines = ["# Caller scorecard summary", "", "Main conclusions:", ""]
    lines.append(
        f"- **{best['caller']}** has the most variants shared across all "
        f"{best['n_ctcs']} CTCs ({best['shared_all_n_ctcs']}), the criterion "
        f"this scorecard weighs most heavily as a proxy for caller reliability."
    )
    if by_jaccard:
        top_j = by_jaccard[0]
        lines.append(
            f"- **{top_j['caller']}** also has the highest mean pairwise "
            f"CTC Jaccard overlap ({top_j['mean_pairwise_jaccard']:.4f})."
        )
    if sparse:
        names = ", ".join(r["caller"] for r in sparse)
        lines.append(
            f"- {names} share essentially no variants across all CTCs "
            f"(shared_all_n_ctcs = 0) and should be treated as too sparse "
            f"for all-CTC trunk inference on this data."
        )
    most_private = max(rows, key=lambda r: r["private_fraction"])
    lines.append(
        f"- **{most_private['caller']}** has the highest private-call "
        f"fraction ({most_private['private_fraction']:.1%}), i.e. the "
        f"largest share of its calls are seen in only one CTC."
    )
    lines.append("")
    lines.append(
        "This is a relative comparison (caller-vs-caller agreement), not an "
        "absolute accuracy measurement against an independent truth set."
    )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(
        description="Per-caller reliability scorecard (cross-CTC sharing, cross-caller "
                    "agreement) plus statistical tests, computed from the comparison matrices."
    )
    ap.add_argument("--matrix-dir", required=True, help="Directory with mutation_binary.tsv, alt_reads.tsv, total_depth.tsv, vaf.tsv")
    ap.add_argument("--output-dir", required=True, help="Directory to write caller_scorecard.tsv, statistical_tests.tsv, scorecard_summary.md")
    args = ap.parse_args()

    matrix_dir = Path(args.matrix_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    columns, binary_rows = read_matrix(matrix_dir / "mutation_binary.tsv")
    _, alt_rows = read_matrix(matrix_dir / "alt_reads.tsv")
    _, depth_rows = read_matrix(matrix_dir / "total_depth.tsv")
    _, vaf_rows = read_matrix(matrix_dir / "vaf.tsv")

    caller_columns_by_var_sample = {}
    for var_id, row in binary_rows.items():
        for col, val in row.items():
            if val != "1":
                continue
            sample_id, caller = split_column(col)
            if not is_ctc_sample(sample_id):
                continue
            caller_columns_by_var_sample.setdefault((var_id, sample_id), set()).add(caller)

    by_caller_ctc_sets, by_caller_calls, background_samples_by_caller, ctc_samples = build_caller_data(
        columns, binary_rows, alt_rows, depth_rows, vaf_rows
    )

    rows = []
    for caller in sorted(by_caller_ctc_sets):
        row = compute_scorecard_row(
            caller,
            by_caller_ctc_sets[caller],
            by_caller_calls[caller],
            background_samples_by_caller[caller],
            binary_rows, alt_rows, depth_rows, vaf_rows,
            caller_columns_by_var_sample,
        )
        rows.append(row)

    write_scorecard(rows, output_dir / "caller_scorecard.tsv")
    write_statistical_tests(rows, output_dir / "statistical_tests.tsv")
    write_summary(rows, output_dir / "scorecard_summary.md")

    print(f"Scorecard: {len(rows)} callers, {len(ctc_samples)} CTC samples")
    print(f"Wrote {output_dir / 'caller_scorecard.tsv'}")
    print(f"Wrote {output_dir / 'statistical_tests.tsv'}")
    print(f"Wrote {output_dir / 'scorecard_summary.md'}")


if __name__ == "__main__":
    main()
