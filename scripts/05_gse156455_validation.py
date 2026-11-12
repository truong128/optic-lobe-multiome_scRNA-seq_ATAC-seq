#!/usr/bin/env python3
"""Independent check of the pDm8/yDm8 genes in GSE156455 (Kurmangaliyev optic-lobe atlas).

Dm8 is one published type there. Cells are split with the same raw-count rule
used on GSE305940: beat-IIIc only = pDm8, DIP-gamma (CG14521) only = yDm8.
"""
import csv
import gzip
import os
from collections import defaultdict

root = os.environ.get("PROJECT_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ext = os.path.join(root, "data", "external")
out = os.path.join(root, "results")
os.makedirs(out, exist_ok=True)

want = {
    "beat-IIIc", "CG15138", "CG14521", "DIP-gamma",
    "nAChRalpha6", "CG44153", "Tl", "Toll", "Sulf1", "Drl-2", "kek6",
}

dm8 = {}
with gzip.open(os.path.join(ext, "GSE156455_metadata_main.tsv.gz"), "rt") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["type"] == "Dm8":
            dm8[row["barcode"]] = row["time"]

barcodes = []
with gzip.open(os.path.join(ext, "GSE156455_barcodes_main.tsv.gz"), "rt") as handle:
    for line in handle:
        barcodes.append(line.rstrip("\n"))
col_of = {}
for i, barcode in enumerate(barcodes, start=1):
    if barcode in dm8:
        col_of[i] = barcode
print(f"Dm8 cells in matrix index: {len(col_of)} of {len(dm8)} metadata rows", flush=True)

row_of = {}
with gzip.open(os.path.join(ext, "GSE156455_features_main.tsv.gz"), "rt") as handle:
    for i, line in enumerate(handle, start=1):
        parts = line.rstrip("\n").split("\t")
        symbol = parts[1] if len(parts) > 1 else ""
        if symbol in want:
            row_of[i] = symbol
print("genes found:", sorted(set(row_of.values())), flush=True)

counts = defaultdict(lambda: defaultdict(float))
matrix = os.path.join(ext, "GSE156455_matrix_main.mtx.gz")
with gzip.open(matrix, "rt") as handle:
    header = handle.readline()
    while header.startswith("%"):
        header = handle.readline()
    nrow, ncol, nnz = (int(x) for x in header.split())
    print(f"matrix {nrow} x {ncol}, {nnz} entries", flush=True)
    seen = 0
    for line in handle:
        i_s, j_s, v_s = line.split()
        i = int(i_s)
        j = int(j_s)
        if i in row_of and j in col_of:
            counts[col_of[j]][row_of[i]] = float(v_s)
            seen += 1
print(f"stored nonzero gene-cell values: {seen}", flush=True)

def gene_count(barcode, *symbols):
    bag = counts.get(barcode, {})
    return sum(bag.get(s, 0.0) for s in symbols)

rows = []
for barcode, time in dm8.items():
    beat = gene_count(barcode, "beat-IIIc", "CG15138")
    dip = gene_count(barcode, "CG14521", "DIP-gamma")
    if beat > 0 and dip == 0:
        label = "pDm8"
    elif dip > 0 and beat == 0:
        label = "yDm8"
    elif beat > 0 and dip > 0:
        label = "double_positive"
    else:
        label = "double_negative"
    rec = {"barcode": barcode, "time": time, "subtype": label, "beat": beat, "dip": dip}
    for symbol in ("nAChRalpha6", "CG44153", "Tl", "Toll", "Sulf1", "Drl-2", "kek6"):
        rec[symbol] = counts.get(barcode, {}).get(symbol, 0.0)
    rows.append(rec)

with open(os.path.join(out, "gse156455_dm8_cells.csv"), "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

genes = ["nAChRalpha6", "CG44153", "Tl", "Toll", "Sulf1", "Drl-2", "kek6"]

def summarize(subset, time_label):
    out_rows = []
    for label in ("pDm8", "yDm8", "double_positive", "double_negative"):
        n = sum(1 for r in subset if r["subtype"] == label)
        out_rows.append({"time": time_label, "gene": "(cells)", "subtype": label, "n": n,
                         "mean": "", "pct": ""})
    p = [r for r in subset if r["subtype"] == "pDm8"]
    y = [r for r in subset if r["subtype"] == "yDm8"]
    for gene in genes:
        for label, group in (("pDm8", p), ("yDm8", y)):
            vals = [r[gene] for r in group]
            n = len(vals)
            out_rows.append({
                "time": time_label,
                "gene": gene,
                "subtype": label,
                "n": n,
                "mean": (sum(vals) / n) if n else "",
                "pct": (sum(v > 0 for v in vals) / n) if n else "",
            })
    return out_rows

summary = []
summary.extend(summarize(rows, "all_times"))
for time in ("24h", "48h", "72h", "96h"):
    summary.extend(summarize([r for r in rows if r["time"] == time], time))

with open(os.path.join(out, "gse156455_gene_summary.csv"), "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["time", "gene", "subtype", "n", "mean", "pct"])
    writer.writeheader()
    writer.writerows(summary)

print("DONE", flush=True)
