#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless-friendly
import matplotlib.pyplot as plt
from pathlib import Path
import sys

X_MAX = 0.05  # plot only points with x <= 0.05

# Inputs produced by the previous run script
RAW_CSV = "leak_results_raw.csv"            # t_min,rep,leaked_photon_ppm
AGG_CSV = "leak_results_aggregated.csv"     # t_min,repeats,mean_ppm,std_ppm
AGG_NPZ = "leak_results_aggregated.npz"     # t_min,repeats,mean_ppm,std_ppm

# Outputs (plotted arrays + figure)
OUT_CSV = "leak_plot_plotted.csv"
OUT_NPZ = "leak_plot_plotted.npz"
OUT_PNG = "leak_plot_points_loglog.png"

def load_from_raw(path):
    """Load per-run ppm, group by t_min; return arrays for mean, stdev, repeats, MC error."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data[None, :]
    t_all = data[:, 0].astype(float)
    ppm_all = data[:, 2].astype(float)

    uniq = np.unique(t_all)
    mean_ppm = np.zeros_like(uniq, dtype=float)
    std_ppm  = np.zeros_like(uniq, dtype=float)
    reps     = np.zeros_like(uniq, dtype=int)
    mc_ppm   = np.zeros_like(uniq, dtype=float)

    for i, tv in enumerate(uniq):
        vals_ppm = ppm_all[t_all == tv]
        reps[i] = vals_ppm.size
        mean_ppm[i] = vals_ppm.mean() if vals_ppm.size else np.nan
        std_ppm[i]  = vals_ppm.std(ddof=1) if vals_ppm.size > 1 else 0.0
        # Monte Carlo error from total hits across repeats: sqrt(sum(hits))/60
        total_hits = np.sum(vals_ppm * 60.0)
        mc_ppm[i] = np.sqrt(max(total_hits, 0.0)) / 60.0

    order = np.argsort(uniq)
    return uniq[order], reps[order], mean_ppm[order], std_ppm[order], mc_ppm[order]

def load_from_agg_csv(path):
    arr = np.loadtxt(path, delimiter=",", skiprows=1)
    if arr.ndim == 1:
        arr = arr[None, :]
    t = arr[:, 0].astype(float)
    r = arr[:, 1].astype(int)
    m = arr[:, 2].astype(float)
    s = arr[:, 3].astype(float)
    order = np.argsort(t)
    t, r, m, s = t[order], r[order], m[order], s[order]
    # Estimate MC error from mean & repeats: total_hits ≈ repeats * mean_ppm * 60
    mc_ppm = np.sqrt(np.clip(r * m * 60.0, 0.0, None)) / 60.0
    return t, r, m, s, mc_ppm

def load_from_agg_npz(path):
    d = np.load(path)
    t = d["t_min"].astype(float)
    r = d["repeats"].astype(int)
    m = d["mean_ppm"].astype(float)
    s = d["std_ppm"].astype(float)
    order = np.argsort(t)
    t, r, m, s = t[order], r[order], m[order], s[order]
    mc_ppm = np.sqrt(np.clip(r * m * 60.0, 0.0, None)) / 60.0
    return t, r, m, s, mc_ppm

def load_data():
    if Path(RAW_CSV).exists():
        return load_from_raw(RAW_CSV)
    if Path(AGG_CSV).exists():
        return load_from_agg_csv(AGG_CSV)
    if Path(AGG_NPZ).exists():
        return load_from_agg_npz(AGG_NPZ)
    print("[!] No input files found (need raw or aggregated leak results).", file=sys.stderr)
    sys.exit(1)

def clipped_err(mean, err, frac=0.999):
    """Clip error so lower bound stays > 0 on log scale."""
    mean = np.asarray(mean, float)
    err  = np.asarray(err, float)
    return np.minimum(err, frac * np.maximum(mean, 1e-300))

def main():
    # Load arrays
    t_vals, reps, mean_ppm, std_ppm, mc_ppm = load_data()

    # First: filter to X <= X_MAX
    mask_x = t_vals <= X_MAX
    if not np.any(mask_x):
        print(f"[!] No points with x <= {X_MAX}.", file=sys.stderr)
        sys.exit(1)
    t_vals = t_vals[mask_x]
    reps    = reps[mask_x]
    mean_ppm = mean_ppm[mask_x]
    std_ppm  = std_ppm[mask_x]
    mc_ppm   = mc_ppm[mask_x]

    # Total Y error (ppm): combine run-to-run stdev and Poisson MC in quadrature
    yerr_ppm = np.sqrt(std_ppm**2 + mc_ppm**2)

    # For log–log plot keep only positive means
    mask_pos = mean_ppm > 0
    if not np.any(mask_pos):
        print("[!] No positive mean values; cannot plot on log scale.", file=sys.stderr)
        sys.exit(1)

    x = t_vals[mask_pos]
    r = reps[mask_pos]
    y = mean_ppm[mask_pos]
    std_sel = std_ppm[mask_pos]
    mc_sel  = mc_ppm[mask_pos]
    yerr = clipped_err(y, yerr_ppm[mask_pos])

    # Save exactly what we plot
    np.savetxt(
        OUT_CSV,
        np.column_stack([x, r, y, std_sel, mc_sel, yerr]),
        delimiter=",",
        header="offset_mm,repeats,mean_leaked_photons(std=ppm),run2run_stdev_ppm,mc_error_ppm,yerr_ppm_combined",
        comments="",
        fmt="%.10g",
    )
    np.savez(
        OUT_NPZ,
        offset_mm=x, repeats=r, mean=y,
        stdev_ppm=std_sel, mc_ppm=mc_sel, yerr_ppm=yerr
    )

    # Plot: markers with vertical error bars; NO line; NO title; log–log
    plt.figure(figsize=(7.5, 5.2))
    plt.errorbar(x, y, yerr=yerr, fmt="o", linestyle="none", capsize=3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Ray offset [mm]")
    plt.ylabel("Leaked photons")
    plt.grid(True, which="both", linestyle="--", alpha=0.4)
    plt.xlim(left=None, right=X_MAX)
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=220)
    print(f"[ok] Saved plot: {OUT_PNG}")
    print(f"[ok] Saved plotted arrays: {OUT_CSV}, {OUT_NPZ}")

if __name__ == "__main__":
    main()
