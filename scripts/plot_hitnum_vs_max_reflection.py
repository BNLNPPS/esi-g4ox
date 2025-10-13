#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from pathlib import Path
import sys

# ---- Inputs (produced by your sweep script) ----
RAW_TXT = "hits_vs_bounce.txt"                 # raw per-run: "bounce hits"
AGG_CSV = "hits_vs_bounce_aggregated.csv"      # optional aggregated CSV
AGG_NPZ = "hits_vs_bounce_aggregated.npz"      # optional aggregated NPZ

# ---- Outputs (what we plot & save for later) ----
OUT_PNG = "hits_vs_bounce_points.png"
OUT_CSV = "hits_vs_bounce_plotted.csv"
OUT_NPZ = "hits_vs_bounce_plotted.npz"

def load_from_csv(path):
    """Load aggregated CSV: columns = bounce, repeats, mean_hits, std_hits."""
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    if arr.ndim == 1:
        arr = arr[None, :]
    if arr.shape[1] < 4:
        raise ValueError(f"Unexpected columns in {path}: {arr.shape}")
    bounce = arr[:, 0].astype(int)
    repeats = arr[:, 1].astype(int)
    mean = arr[:, 2].astype(float)
    stdev = arr[:, 3].astype(float)
    order = np.argsort(bounce)
    return bounce[order], repeats[order], mean[order], stdev[order]

def load_from_npz(path):
    d = np.load(path)
    bounce = d["bounce"].astype(int)
    repeats = d["repeats"].astype(int)
    mean = d["mean"].astype(float)
    stdev = d["stdev"].astype(float)
    order = np.argsort(bounce)
    return bounce[order], repeats[order], mean[order], stdev[order]

def aggregate_from_raw(path):
    """Aggregate per-bounce mean & sample stdev from raw TXT."""
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    xs = data[:, 0].astype(int)
    ys = data[:, 1].astype(float)
    uniq = np.unique(xs)
    mean = np.zeros_like(uniq, dtype=float)
    stdev = np.zeros_like(uniq, dtype=float)
    repeats = np.zeros_like(uniq, dtype=int)
    for i, b in enumerate(uniq):
        v = ys[xs == b]
        repeats[i] = v.size
        mean[i] = np.mean(v) if v.size else np.nan
        stdev[i] = np.std(v, ddof=1) if v.size > 1 else 0.0
    return uniq, repeats, mean, stdev

def load_data():
    if Path(AGG_CSV).exists():
        return load_from_csv(AGG_CSV)
    if Path(AGG_NPZ).exists():
        return load_from_npz(AGG_NPZ)
    if Path(RAW_TXT).exists():
        return aggregate_from_raw(RAW_TXT)
    print("[!] No input files found.", file=sys.stderr)
    sys.exit(1)

def fixed_e6_formatter(y, _pos):
    """Format ticks as mantissa × 10^6 (always with exponent 6 in mathtext)."""
    if y == 0:
        return "0"
    m = y / 1e6
    return rf"${m:.3g}\times 10^{{6}}$"

def main():
    bounce, repeats, mean, stdev = load_data()

    # Save exactly what we will plot (for later reuse)
    np.savetxt(
        OUT_CSV,
        np.column_stack([bounce, repeats, mean, stdev]),
        delimiter=",",
        header="bounce,repeats,mean_hits,std_hits",
        comments="",
        fmt="%.10g",
    )
    np.savez(OUT_NPZ, bounce=bounce, repeats=repeats, mean=mean, stdev=stdev)

    # Plot: points with vertical error bars; NO line; NO title
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.errorbar(bounce, mean, yerr=stdev, fmt="o", linestyle="none", capsize=3)

    # Axis labels
    ax.set_xlabel("Maximum allowed photon reflection")
    ax.set_ylabel("Number of hits (Opticks)")

    # Y-axis ticks as mantissa × 10^6 with exponent using mathtext
    ax.yaxis.set_major_formatter(FuncFormatter(fixed_e6_formatter))

    # X-axis ticks: start from 5 and show every third
    if bounce.size:
        xmax = int(np.max(bounce))
        ticks = np.arange(5, xmax + 1, 3, dtype=int)
        ax.set_xticks(ticks)

    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=220)
    print(f"[ok] Saved plot: {OUT_PNG}")
    print(f"[ok] Saved plotted arrays: {OUT_CSV}, {OUT_NPZ}")

if __name__ == "__main__":
    main()
