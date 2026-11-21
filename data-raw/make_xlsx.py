"""Assemble the per-analyte Bio-Plex example workbook from CSVs.

Run from the package root after data-raw/make_examples.R:
    python3 data-raw/make_xlsx.py
"""
import csv
import glob
import os

from openpyxl import Workbook

src = "data-raw/per_analyte"
out = "inst/extdata/bioplex_per_analyte_example.xlsx"

wb = Workbook()
wb.remove(wb.active)
for path in sorted(glob.glob(os.path.join(src, "*.csv"))):
    name = os.path.splitext(os.path.basename(path))[0][:31]
    ws = wb.create_sheet(title=name)
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.reader(fh):
            ws.append(row)
wb.save(out)
print("wrote", out)
