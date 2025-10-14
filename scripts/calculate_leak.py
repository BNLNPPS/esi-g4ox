# plot error is not correct here need to check plotting script for that!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!/usr/bin/env python3
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import re
import os
from pathlib import Path

# --- sweep setup ---
REPEATS = 10
# 50 logarithmically spaced t_min values between 1e-5 and 1
t_values = np.logspace(-5, 0, 50)

EXE = "./build/src/simg4ox"
GDML = "esi-g4ox/geom/sphere_leak.gdml"
MAC  = "esi-g4ox/run.mac"

# --- outputs ---
RAW_CSV = "leak_results_raw.csv"           # per-run rows
AGG_CSV = "leak_results_aggregated.csv"    # arrays that will be plotted
AGG_NPZ = "leak_results_aggregated.npz"
PLOT_PNG = "leak_plot_points.png"

# --- helpers ---
def run_once(t_min: float) -> float:
    """Run one simulation, return leaked photons in ppm (NumHits/60.0)."""
    env = os.environ.copy()
    env["OPTICKS_PROPAGATE_EPSILON"] = f"{t_min:.10g}"
    try:
        proc = subprocess.run(
            [EXE, "-g", GDML, "-m", MAC],
            env=env, capture_output=True, text=True, check=True
        )
        combined = (proc.stdout or "") + "\n" + (proc.stderr or "")
        m = re.search(r"Opticks:\s*NumHits:\s*(\d+)", combined)
        if not m:
            print(f"[!] NumHits not found for t_min={t_min:.3g}; recording 0")
            return 0.0
        num_hits = int(m.group(1))
        return num_hits / 60.0  # ppm (as in your original script)
    except subprocess.CalledProcessError as e:
        print(f"[!] Simulation error for t_min={t_min:.3g}: {e}; recording 0")
        return 0.0

def clipped_sym_err(mean, std, frac=0.999):
    """Prevent negative lower error on log-y plots."""
    mean = np.asarray(mean, float)
    std  = np.asarray(std, float)
    return np.minimum(std, frac * np.maximum(mean, 1e-300))

# --- main sweep (collect raw per-run data) ---
raw_rows = []  # (t_min, rep_idx, ppm)
for t in t_values:
    for rep in range(1, REPEATS + 1):
        ppm = run_once(t)
        raw_rows.append((t, rep, ppm))
        print(f"[ok] t_min={t:.5g} rep={rep}/{REPEATS}  ppm={ppm:.6g}")

# save raw results
with open(RAW_CSV, "w") as f:
    f.write("t_min,rep,leaked_photon_ppm\n")
    for t, rep, ppm in raw_rows:
        f.write(f"{t:.10g},{rep},{ppm:.10g}\n")
print(f"[ok] saved raw per-run CSV: {RAW_CSV}")

# --- aggregate by t_min (mean & sample stdev) ---
t_arr = np.array([r[0] for r in raw_rows], float)
ppm_arr = np.array([r[2] for r in raw_rows], float)
uniq = np.unique(t_arr)
mean = np.zeros_like(uniq)
stdev = np.zeros_like(uniq)
count = np.zeros_like(uniq, int)

for i, tv in enumerate(uniq):
    vals = ppm_arr[t_arr == tv]
    count[i] = vals.size
    mean[i] = np.mean(vals) if vals.size else np.nan
    stdev[i] = np.std(vals, ddof=1) if vals.size > 1 else 0.0  # sample stdev

# save arrays that WILL be plotted
np.savetxt(
    AGG_CSV,
    np.column_stack([uniq, count, mean, stdev]),
    delimiter=",",
    header="t_min,repeats,mean_ppm,std_ppm",
    comments="",
    fmt="%.10g",
)
np.savez(AGG_NPZ, t_min=uniq, repeats=count, mean_ppm=mean, std_ppm=stdev)
print(f"[ok] saved aggregated CSV: {AGG_CSV}")
print(f"[ok] saved aggregated NPZ: {AGG_NPZ}")

# --- plotting: points + vertical error bars, no line, no title, log Y ---
# mask out non-positive means for log-y plotting
mask = mean > 0
if not np.any(mask):
    raise SystemExit("[!] No positive mean values; nothing to plot on log Y.")

x = uniq[mask]
y = mean[mask]
yerr = clipped_sym_err(y, stdev[mask])

plt.figure(figsize=(7.2, 5.0))
plt.errorbar(x, y, yerr=yerr, fmt="o", linestyle="none", capsize=3)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("X: t_min")
plt.ylabel("Y: leaked photon number [ppm]")
plt.grid(True, which="both", linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig(PLOT_PNG, dpi=220)
plt.close()
print(f"[ok] saved plot: {PLOT_PNG}")

