#!/usr/bin/env python3
# Save two plots from Opticks.txt, both with log-scaled X:
# Need to run batching_performance_run.py beforehand to produce the to be plotted
#  1) time_per_photon_opticks.png  — Y (µs) on log-scale
#  2) speedup_batching.png         — baseline = first Y value, Y on log-scale

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless backend for saving PNGs without a display
import matplotlib.pyplot as plt
import sys
from pathlib import Path

INPUT_FILE = "Opticks.txt"
OUT1 = "time_per_photon_opticks.png"
OUT2 = "speedup_batching.png"

def main():
    p = Path(INPUT_FILE)
    if not p.exists():
        print(f"[!] File not found: {INPUT_FILE}", file=sys.stderr)
        sys.exit(1)

    # Expect two whitespace-separated columns: X (events), Y (time per photon [s])
    try:
        data = np.loadtxt(INPUT_FILE, dtype=float)
    except Exception as e:
        print(f"[!] Could not read {INPUT_FILE}: {e}", file=sys.stderr)
        sys.exit(1)

    if data.ndim != 2 or data.shape[1] < 2:
        print(f"[!] Expected 2 columns in {INPUT_FILE} (got shape {data.shape})", file=sys.stderr)
        sys.exit(1)

    # Sort by X to make log-scale plots clean
    order = np.argsort(data[:, 0])
    data = data[order]

    x = data[:, 0]
    y = data[:, 1]  # seconds per photon

    # -------- Plot 1: X vs (Y * 1e6) (log-scale on X and Y) --------
    y_us = y * 1e6  # convert seconds → microseconds
    plt.figure(figsize=(7.5, 4.8))
    plt.plot(x, y_us, marker="o")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("number of Geant4 Events in a single GPU call")
    plt.ylabel("Simulation time per photon [µs]") 
    plt.title("Per-photon simulation time vs. batch size")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUT1, dpi=200)
    plt.close()

    # -------- Plot 2: Speedup relative to first Y value (log-scale on X and Y) --------
    baseline = y[0]
    if baseline <= 0:
        print("[!] First Y value must be positive to compute speedup.", file=sys.stderr)
        sys.exit(1)

    # Keep only positive Y entries (needed for log-scale and division)
    mask = y > 0
    x2 = x[mask]
    speedup = baseline / y[mask]

    plt.figure(figsize=(7.5, 4.8))
    plt.plot(x2, speedup, marker="o")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Number of Geant4 Events")
    plt.ylabel("Speedup due to batching")
    plt.title("Batching speedup (relative to simulating one event in a GPU call)")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUT2, dpi=200)
    plt.close()

    print(f"[ok] Saved: {OUT1} and {OUT2}")

if __name__ == "__main__":
    main()
