#!/usr/bin/env python3
"""Merge VEP tab outputs while keeping one normalized header."""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args()

    header_written = False
    with open(args.output, "w", newline="") as out:
        for item in args.inputs:
            with open(item, "r", newline="") as handle:
                for line in handle:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#"):
                        if not header_written:
                            out.write(line[1:])
                            header_written = True
                        continue
                    if line.strip():
                        out.write(line)

    if not header_written:
        fields = [
            "Uploaded_variation", "Location", "Allele", "Gene", "SYMBOL", "Feature",
            "Consequence", "IMPACT", "HGVSc", "HGVSp", "Existing_variation", "MAX_AF",
            "MAX_AF_POPS", "AF", "gnomADe_AF", "gnomADg_AF",
        ]
        Path(args.output).write_text("\t".join(fields) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
