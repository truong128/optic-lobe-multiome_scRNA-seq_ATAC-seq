#!/usr/bin/env python3
"""Verify published LinkPeaks labels, and retest GSE156455 at the genotype level."""
import csv
import gzip
import os
import statistics
from collections import defaultdict

import openpyxl

root = os.environ.get("PROJECT_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = os.path.join(root, "results")
cand_path = os.path.join(out, "candidate_enhancers_revised.csv")

published = {}
wb = openpyxl.load_workbook(
    os.path.join(root, "data", "GSE305940_Links.xlsx"), read_only=True, data_only=True
)
for stage, sheet in (("P24", "P24_links"), ("P48", "P48_links")):
    ws = wb[sheet]
    header = [c.value for c in next(ws.iter_rows(min_row=1, max_row=1))]
    idx = {name: i for i, name in enumerate(header)}
    for row in ws.iter_rows(min_row=2, values_only=True):
        try:
            pvalue = float(row[idx["pvalue"]])
        except (TypeError, ValueError):
            continue
        if pvalue >= 0.05:
            continue
        gene = row[idx["gene"]]
        peak = row[idx["peak"]]
        published.setdefault((stage, peak), set()).add(gene)
wb.close()

alias = {}
with open(os.path.join(root, "results", "gene_alias_map.csv")) as handle:
    for row in csv.DictReader(handle):
        alias[row["alias"]] = row["symbol"]

def symbol(name):
    return alias.get(name, name)

rows = []
with open(cand_path) as handle:
    for row in csv.DictReader(handle):
        genes = published.get((row["stage"], row["peak"]), set())
        symbols = {symbol(g) for g in genes}
        claimed = "published_LinkPeaks" in row["link"]
        confirmed = row["gene"] in symbols
        rows.append({
            "stage": row["stage"],
            "peak": row["peak"],
            "gene": row["gene"],
            "claimed_published": claimed,
            "confirmed_in_workbook": confirmed,
            "workbook_genes": ";".join(sorted(symbols)),
        })

with open(os.path.join(out, "published_link_check.csv"), "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
claimed = [r for r in rows if r["claimed_published"]]
print(f"claimed published {len(claimed)}; confirmed {sum(r['confirmed_in_workbook'] for r in claimed)}")

# Genotype-level means in GSE156455.
meta = {}
with gzip.open(os.path.join(root, "data", "external", "GSE156455_metadata_main.tsv.gz"), "rt") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["type"] == "Dm8":
            meta[row["barcode"]] = row

cells = []
with open(os.path.join(out, "gse156455_dm8_cells.csv")) as handle:
    for row in csv.DictReader(handle):
        info = meta.get(row["barcode"])
        if info is None:
            continue
        row["genotype"] = info["genotype"]
        row["rep"] = info["rep"]
        cells.append(row)

genes = ["nAChRalpha6", "CG44153", "Tl", "Sulf1", "Drl-2", "kek6"]
summary = []
for time in ("24h", "48h"):
    sub = [r for r in cells if r["time"] == time and r["subtype"] in ("pDm8", "yDm8")]
    genotypes = sorted({r["genotype"] for r in sub})
    for gene in genes:
        deltas = []
        used = 0
        for geno in genotypes:
            p = [float(r[gene]) for r in sub if r["genotype"] == geno and r["subtype"] == "pDm8"]
            y = [float(r[gene]) for r in sub if r["genotype"] == geno and r["subtype"] == "yDm8"]
            if len(p) < 3 or len(y) < 3:
                continue
            used += 1
            # log1p mean, y minus p: positive means higher in yDm8
            import math
            mp = sum(math.log1p(v) for v in p) / len(p)
            my = sum(math.log1p(v) for v in y) / len(y)
            deltas.append(my - mp)
        n_pos = sum(d > 0 for d in deltas)
        summary.append({
            "time": time,
            "gene": gene,
            "n_genotypes": used,
            "n_higher_in_yDm8": n_pos,
            "median_log1p_diff_y_minus_p": (statistics.median(deltas) if deltas else ""),
        })

with open(os.path.join(out, "gse156455_by_genotype.csv"), "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()))
    writer.writeheader()
    writer.writerows(summary)
for row in summary:
    print(row)
print("DONE")
