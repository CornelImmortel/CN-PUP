#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


BRANCH_RE = re.compile(r":(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)")


def quantile(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    values = sorted(values)
    pos = (len(values) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return values[lo]
    return values[lo] * (hi - pos) + values[hi] * (pos - lo)


def main() -> None:
    ap = argparse.ArgumentParser(description="Clip extreme Newick/NEXUS branch lengths for plotting only.")
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--quantile", type=float, default=0.75)
    ap.add_argument("--multiplier", type=float, default=10.0)
    ap.add_argument("--min-cap", type=float, default=1.0)
    args = ap.parse_args()

    text = Path(args.input).read_text(encoding="utf-8", errors="replace")
    matches = list(BRANCH_RE.finditer(text))
    lengths = [float(m.group(1)) for m in matches if math.isfinite(float(m.group(1))) and float(m.group(1)) > 0]
    baseline = quantile(lengths, args.quantile)
    cap = max(args.min_cap, baseline * args.multiplier)

    clipped = 0
    max_original = max(lengths) if lengths else 0.0

    def replace(match: re.Match[str]) -> str:
        nonlocal clipped
        value = float(match.group(1))
        if math.isfinite(value) and value > cap:
            clipped += 1
            return f":{cap:.10g}"
        return match.group(0)

    Path(args.output).write_text(BRANCH_RE.sub(replace, text), encoding="utf-8")
    Path(args.summary).write_text(
        "\t".join(["branch_count", "positive_branch_count", "max_original_branch_length", "clip_quantile", "clip_multiplier", "clip_cap", "clipped_branch_count"]) + "\n"
        + "\t".join([
            str(len(matches)),
            str(len(lengths)),
            f"{max_original:.10g}",
            str(args.quantile),
            str(args.multiplier),
            f"{cap:.10g}",
            str(clipped),
        ])
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
