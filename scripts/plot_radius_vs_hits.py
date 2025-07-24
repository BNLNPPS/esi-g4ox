#!/usr/bin/env python3
"""
Pandas-free plot of sphere diameter vs. photon hits, with √N error bars.

Works even if the CSV header has weird spacing/capitalisation, or if you
just have two columns with no header at all.
"""

import csv
import math
from pathlib import Path
import matplotlib.pyplot as plt

file_path = Path("sphere_hits.csv")     # adjust if needed
marker_style = "o"

# ── STEP 1: read the file ─────────────────────────────────────────────
diameters_mm, hits = [], []

with file_path.open(newline="") as f:
    peek = f.readline()
    f.seek(0)

    # Does the first non–newline char look like a digit? If yes → no header
    first_token = peek.strip().split(",")[0]
    has_header = not first_token.replace(".", "", 1).isdigit()

    if has_header:
        reader = csv.reader(f)
        raw_header = next(reader)
        header = [h.strip().lower() for h in raw_header]

        # try to locate diameter & hits columns by common names
        try:
            x_idx = header.index("rmin")
        except ValueError:
            try:
                x_idx = header.index("radius")
            except ValueError:
                x_idx = 0  # guess: first column
                print("⚠️  Unrecognised diameter column; assuming column 1")

        try:
            y_idx = header.index("photonhits")
        except ValueError:
            try:
                y_idx = header.index("hits")
            except ValueError:
                y_idx = 1  # guess: second column
                print("⚠️  Unrecognised hits column; assuming column 2")

        for row in reader:
            if not row:               # skip empty lines
                continue
            diameters_mm.append(float(row[x_idx]))
            hits.append(int(row[y_idx]))
    else:
        # No header: treat whole file as numeric rows
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 2 or not row[0].strip():
                continue
            diameters_mm.append(float(row[0]))
            hits.append(int(row[1]))

# ── STEP 2: compute errors and plot ────────────────────────────────────
y_err = [math.sqrt(n) for n in hits]

plt.errorbar(
    diameters_mm, hits,
    yerr=y_err,
    fmt=marker_style,
    linestyle=""
)
plt.xlabel("radius of spherical mirror [mm]")
plt.ylabel("hits")
plt.tight_layout()
plt.show()
