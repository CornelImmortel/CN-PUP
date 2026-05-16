#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip, re
from pathlib import Path

ANN_FIELDS = ["allele","effect","impact","gene_symbol","gene_id","feature_type","transcript","biotype","rank","hgvs_c","hgvs_p"]

def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")

def split_info(s):
    out = {}
    for item in s.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1); out[k] = v
        else:
            out[item] = True
    return out

def parse_sample(fmt, sample):
    if not fmt or not sample:
        return {}
    keys, vals = fmt.split(":"), sample.split(":")
    return {k: vals[i] if i < len(vals) else "" for i, k in enumerate(keys)}

def num(x):
    try:
        if x in (None, "", "."):
            return ""
        return str(float(x)).rstrip("0").rstrip(".")
    except ValueError:
        return ""

def ann_best(info, alt):
    raw = str(info.get("ANN", ""))
    anns = []
    for entry in raw.split(","):
        if not entry:
            continue
        parts = entry.split("|")
        parts += [""] * max(0, 16 - len(parts))
        d = dict(zip(ANN_FIELDS, parts[:len(ANN_FIELDS)]))
        d["effect"] = d.get("effect") or "unannotated"
        d["impact"] = d.get("impact") or "UNANNOTATED"
        d["gene_symbol"] = d.get("gene_symbol") or "Unknown"
        anns.append(d)
    matching = [x for x in anns if x.get("allele") == alt] or anns
    order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3, "UNANNOTATED": 4}
    if not matching:
        return {k: "" for k in ANN_FIELDS} | {"effect": "unannotated", "impact": "UNANNOTATED", "gene_symbol": "Unknown"}
    return min(matching, key=lambda x: (order.get(x.get("impact", ""), 9), 0 if x.get("biotype") == "protein_coding" else 1, x.get("gene_symbol", "")))

def read_vcf(path: Path, sample_id: str, caller: str, source_vcf: str):
    rows = []
    sample_name = sample_id
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#CHROM"):
                cols = line.split("\t")
                if len(cols) > 9 and sample_id == "AUTO":
                    sample_name = cols[9]
                continue
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 8:
                continue
            chrom, pos, _id, ref, alts, qual, filt, info_s = p[:8]
            fmt = p[8] if len(p) > 8 else ""
            sm = parse_sample(fmt, p[9]) if len(p) > 9 else {}
            im = split_info(info_s)
            ad = sm.get("AD", "")
            ad_vals = []
            for x in ad.split(","):
                try: ad_vals.append(int(float(x)))
                except ValueError: ad_vals.append(None)
            refread = ad_vals[0] if len(ad_vals) > 0 and ad_vals[0] is not None else ""
            total_depth = sm.get("DP") or im.get("DP") or (sum(x for x in ad_vals if x is not None) if ad_vals else "")
            for alt_i, alt in enumerate(alts.split(","), start=1):
                if len(ref) != 1 or len(alt) != 1 or alt in {".", "*"}:
                    continue
                altread = ad_vals[alt_i] if len(ad_vals) > alt_i and ad_vals[alt_i] is not None else ""
                try:
                    vaf = float(altread) / float(total_depth) if str(altread) != "" and str(total_depth) != "" and float(total_depth) > 0 else ""
                except ValueError:
                    vaf = ""
                ann = ann_best(im, alt)
                rows.append({
                    "sample_id": sample_name, "caller": caller,
                    "var_id": f"{chrom}:{pos}_{ref}/{alt}", "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                    "qual": qual, "filter": filt, "gt": sm.get("GT", ""),
                    "refread": refread, "altread": altread, "total_depth": total_depth,
                    "vaf": vaf, "gq": sm.get("GQ", ""),
                    "gene_symbol": ann.get("gene_symbol", "Unknown"), "gene_id": ann.get("gene_id", ""),
                    "effect": ann.get("effect", "unannotated"), "impact": ann.get("impact", "UNANNOTATED"),
                    "transcript": ann.get("transcript", ""), "hgvs_c": ann.get("hgvs_c", ""), "hgvs_p": ann.get("hgvs_p", ""),
                    "source_vcf": source_vcf,
                })
    return rows

def write_gz_tsv(path: Path, rows):
    cols = ["sample_id","caller","var_id","chrom","pos","ref","alt","qual","filter","gt","refread","altread","total_depth","vaf","gq","gene_symbol","gene_id","effect","impact","transcript","hgvs_c","hgvs_p","source_vcf"]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

def infer_sample_caller(path: Path):
    caller = path.parent.name
    m = re.match(r"(.+)\.([^.]+)\.comparable\.vcf\.gz$", path.name)
    sample = m.group(1) if m else path.stem.split(".")[0]
    return sample, caller

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--input-root", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()
    in_root, out_dir = Path(args.input_root), Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    vcfs = sorted(in_root.glob("*/*.comparable.vcf.gz"))
    if not vcfs:
        raise SystemExit(f"No comparable VCFs found under {in_root}")
    for vcf in vcfs:
        sample, caller = infer_sample_caller(vcf)
        rows = read_vcf(vcf, sample, caller, str(vcf))
        out = out_dir / f"{sample}.{caller}.long.tsv.gz"
        write_gz_tsv(out, rows)
        print(f"Wrote {out} ({len(rows)} rows)")
if __name__ == "__main__":
    main()
