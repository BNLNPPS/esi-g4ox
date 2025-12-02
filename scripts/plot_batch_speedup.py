import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless backend for saving PNGs without a display
import matplotlib.pyplot as plt
import sys
from pathlib import Path

INPUT_FILE = "Opticks.txt"
OUT1 = "time_per_photon_opticks.png"
OUT2 = "speedup_batching.png"
CSV_OUT = "Opticks_aggregated.csv"
NPZ_OUT = "Opticks_aggregated.npz"

# Scale factor for fonts, numbers, line widths, etc.
STYLE_SCALE = 1.6  # keep this

def scale_style(factor=2.0):
    """
    Scale fonts, tick label sizes, and common line widths by 'factor'.

    Works even if some rcParams are strings (e.g. 'medium') by
    falling back to numeric defaults in that case.
    """
    def _scaled(key, fallback):
        val = matplotlib.rcParams.get(key, fallback)
        if isinstance(val, (int, float)):
            return val * factor
        else:
            # e.g. 'medium', 'large' → use numeric fallback instead
            return fallback * factor

    matplotlib.rcParams.update({
        # Fonts
        "font.size":        _scaled("font.size", 10.0),
        "axes.titlesize":   _scaled("axes.titlesize", 12.0),
        "axes.labelsize":   _scaled("axes.labelsize", 10.0),
        "xtick.labelsize":  _scaled("xtick.labelsize", 10.0),
        "ytick.labelsize":  _scaled("ytick.labelsize", 10.0),
        "legend.fontsize":  _scaled("legend.fontsize", 10.0),

        # Lines & markers
        "lines.linewidth":  _scaled("lines.linewidth", 1.5),
        "lines.markersize": _scaled("lines.markersize", 6.0),
        "axes.linewidth":   _scaled("axes.linewidth", 1.0),
        "grid.linewidth":   _scaled("grid.linewidth", 0.8),

        # Tick width/length
        "xtick.major.width": _scaled("xtick.major.width", 0.8),
        "ytick.major.width": _scaled("ytick.major.width", 0.8),
        "xtick.major.size":  _scaled("xtick.major.size", 3.5),
        "ytick.major.size":  _scaled("ytick.major.size", 3.5),
    })


def load_data(path):
    try:
        data = np.loadtxt(path, dtype=float)
    except Exception as e:
        print(f"[!] Could not read {path}: {e}", file=sys.stderr)
        sys.exit(1)
    if data.ndim != 2 or data.shape[1] < 2:
        print(f"[!] Expected 2 columns in {path} (got shape {data.shape})", file=sys.stderr)
        sys.exit(1)
    return data

def aggregate_by_x(x_all, y_all):
    """Group repeated Y values for identical X and compute mean, sample stdev."""
    # Ensure integer-ish X for grouping
    x_all = x_all.astype(int)
    uniq = np.unique(x_all)
    means = np.zeros_like(uniq, dtype=float)
    stds  = np.zeros_like(uniq, dtype=float)
    counts = np.zeros_like(uniq, dtype=int)
    for i, xv in enumerate(uniq):
        ys = y_all[x_all == xv]
        counts[i] = ys.size
        means[i]  = np.mean(ys) if ys.size else np.nan
        stds[i]   = np.std(ys, ddof=1) if ys.size > 1 else 0.0
    return uniq.astype(float), means, stds, counts

def clipped_sym_err(mean, std, frac=0.999):
    """Clip symmetric error so the lower bound stays > 0 for log-scale plotting."""
    # If mean==0 this will become 0; caller should mask out nonpositive values.
    return np.minimum(std, frac * np.maximum(mean, 1e-300))

def save_plotted_data_csv(x, y_mean, y_std, y_mean_us, y_std_us, speed_mean, speed_std, counts, path):
    header = (
        "events,repeats,mean_sec_per_photon,std_sec_per_photon,"
        "mean_us_per_photon,std_us_per_photon,speedup_mean,speedup_std"
    )
    arr = np.column_stack([x, counts, y_mean, y_std, y_mean_us, y_std_us, speed_mean, speed_std])
    np.savetxt(path, arr, delimiter=",", header=header, comments="", fmt="%.10g")

def main():
    # Apply 2x styling
    scale_style(STYLE_SCALE)

    p = Path(INPUT_FILE)
    if not p.exists():
        print(f"[!] File not found: {INPUT_FILE}", file=sys.stderr)
        sys.exit(1)

    data = load_data(INPUT_FILE)
    # Sort by X for clean plotting
    order = np.argsort(data[:, 0])
    data = data[order]

    x_all = data[:, 0]
    y_all = data[:, 1]  # seconds per photon (repeated)

    # ---- Aggregate repeats ----
    x, y_mean, y_std, counts = aggregate_by_x(x_all, y_all)

    # ---- Plot 1: time per photon (µs) with error bars, points only ----
    mask_pos = y_mean > 0
    if not np.any(mask_pos):
        print("[!] No positive mean Y values; cannot make log-scale plot.", file=sys.stderr)
        sys.exit(1)

    x1 = x[mask_pos]
    y_mean_us = y_mean[mask_pos] * 1e6
    y_std_us  = y_std[mask_pos]  * 1e6
    yerr_us   = clipped_sym_err(y_mean_us, y_std_us)

    # 2x capsize (was 3)
    capsize = 3 * STYLE_SCALE

    plt.figure(figsize=(7.5, 4.8))
    plt.errorbar(
        x1, y_mean_us, yerr=yerr_us,
        fmt="o",
        linestyle="none",
        capsize=capsize
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("number of Geant4 Events in a single GPU call")
    plt.ylabel("Simulation time per photon [µs]")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUT1, dpi=200)
    plt.close()

    # ---- Plot 2: speedup (baseline = first mean) with error bars, points only ----
    # Baseline group = first (smallest X) mean
    if y_mean[0] <= 0:
        print("[!] First mean Y must be positive to compute speedup.", file=sys.stderr)
        sys.exit(1)

    b  = y_mean[0]
    sb = y_std[0]

    m2_mask = y_mean > 0  # only positive means allowed for ratio/log-scale
    x2 = x[m2_mask]
    m  = y_mean[m2_mask]
    s  = y_std[m2_mask]

    speedup_mean = b / m
    # First-order error propagation for S = b/m:
    # Var(S) ≈ (S^2)*[(sb/b)^2 + (s/m)^2]
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_var = (sb / b)**2 + (s / m)**2
    rel_var = np.nan_to_num(rel_var, nan=np.inf, posinf=np.inf, neginf=np.inf)
    speedup_std = speedup_mean * np.sqrt(rel_var)
    speedup_err = clipped_sym_err(speedup_mean, speedup_std)

    plt.figure(figsize=(7.5, 4.8))
    plt.errorbar(
        x2, speedup_mean, yerr=speedup_err,
        fmt="o",
        linestyle="none",
        capsize=capsize
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Number of Geant4 Events")
    plt.ylabel("Speedup due to batching")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUT2, dpi=200)
    plt.close()

    # ---- Save arrays used for plotting ----
    # Fill full-length speedup arrays (NaN where invalid)
    speedup_mean_full = np.where(m2_mask, b / y_mean, np.nan)
    # recompute std on full length to align CSV columns
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_var_full = (sb / b)**2 + (y_std / y_mean)**2
    rel_var_full = np.nan_to_num(rel_var_full, nan=np.inf, posinf=np.inf, neginf=np.inf)
    speedup_std_full = np.where(m2_mask, (b / y_mean) * np.sqrt(rel_var_full), np.nan)

    save_plotted_data_csv(
        x, y_mean, y_std,
        y_mean * 1e6, y_std * 1e6,
        speedup_mean_full, speedup_std_full,
        counts, CSV_OUT
    )
    np.savez(
        NPZ_OUT,
        x=x,
        repeats=counts,
        y_mean=y_mean,
        y_std=y_std,
        y_mean_us=y_mean * 1e6,
        y_std_us=y_std * 1e6,
        speedup_mean=speedup_mean_full,
        speedup_std=speedup_std_full,
        # convenience: the exact arrays used for each plot
        x_plot_time=x1,
        y_mean_us_plot=y_mean_us,
        yerr_us_plot=yerr_us,
        x_plot_speed=x2,
        speedup_mean_plot=speedup_mean,
        speedup_err_plot=speedup_err,
    )

    print(f"[ok] Saved: {OUT1}, {OUT2}")
    print(f"[ok] Plotted data: {CSV_OUT}, {NPZ_OUT}")

if __name__ == "__main__":
    main()
