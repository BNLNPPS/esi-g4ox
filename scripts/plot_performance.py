import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import defaultdict

# --- Config ---
N_RUNS = 10
G4_FILE_TEMPLATE = "timings{idx}.txt"
OPTICKS_FILE = "Opticks.txt"
USE_SEM = False               # True -> standard error of the mean; False -> sample std dev
ERR_STYLE = "linear_on_log"   # "linear_on_log" or "log_symmetric"

def parse_thread_value_lines(path):
    """Return list of (thread, value) from 't val' lines; ignores malformed lines."""
    data = []
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            try:
                t = int(parts[0]); v = float(parts[1])
                data.append((t, v))
            except ValueError:
                continue
    return data

# --- Load Geant4 data from 10 files ---
g4_samples = defaultdict(list)  # thread -> list of values across runs
missing = []
for i in range(N_RUNS):
    fn = G4_FILE_TEMPLATE.format(idx=i)
    if not os.path.exists(fn):
        missing.append(fn)
        continue
    for t, v in parse_thread_value_lines(fn):
        g4_samples[t].append(v)

if missing:
    print(f"[!] Missing G4 files: {', '.join(missing)}")

# --- Load Opticks data (take last 10 values per thread if more exist) ---
opticks_samples = defaultdict(list)
if not os.path.exists(OPTICKS_FILE):
    raise FileNotFoundError(f"{OPTICKS_FILE} not found")

for t, v in parse_thread_value_lines(OPTICKS_FILE):
    opticks_samples[t].append(v)

# keep only last N_RUNS per thread (most recent 10)
for t in list(opticks_samples.keys()):
    if len(opticks_samples[t]) >= N_RUNS:
        opticks_samples[t] = opticks_samples[t][-N_RUNS:]

# --- Align threads present in both sets ---
threads = sorted(set(g4_samples.keys()) & set(opticks_samples.keys()))
if not threads:
    raise RuntimeError("No overlapping thread IDs between G4 and Opticks data.")

def mean_and_err(vals):
    arr = np.asarray(vals, dtype=float)
    mean = float(np.mean(arr)) if arr.size > 0 else float("nan")
    if arr.size > 1:
        std = float(np.std(arr, ddof=1))
    else:
        std = 0.0
    err = std / np.sqrt(arr.size) if (USE_SEM and arr.size > 0) else std
    return mean, err, int(arr.size)

# --- Compute per-thread stats ---
g4_mean, g4_err, opt_mean, opt_err = [], [], [], []
for t in threads:
    m, e, _ = mean_and_err(g4_samples[t])
    g4_mean.append(m); g4_err.append(e)
    m2, e2, _ = mean_and_err(opticks_samples[t])
    opt_mean.append(m2); opt_err.append(e2)

g4_mean = np.array(g4_mean, dtype=float)
g4_err  = np.array(g4_err,  dtype=float)
opt_mean = np.array(opt_mean, dtype=float)
opt_err  = np.array(opt_err,  dtype=float)

# --- Ratio and error propagation ---
ratio = g4_mean / opt_mean
rel_err = np.sqrt((g4_err / g4_mean)**2 + (opt_err / opt_mean)**2)

if ERR_STYLE == "log_symmetric":
    # equal in log-space (multiplicative)
    yerr_down = ratio * (1.0 - np.exp(-rel_err))
    yerr_up   = ratio * (np.exp(rel_err) - 1.0)
else:
    # symmetric in linear space (appears longer downward on a log axis)
    yerr_lin  = ratio * rel_err
    # guard to avoid zero/negative lower bound on a log axis
    yerr_down = np.minimum(yerr_lin, ratio * 0.999999)
    yerr_up   = yerr_lin

yerr = np.vstack([yerr_down, yerr_up])

# --- Print a stats table and also save it ---
header = "thread  G4_mean  G4_err  Opt_mean  Opt_err  Ratio  rel_err  yerr_down  yerr_up"
print(header)
lines = [header]
for t, a, ea, b, eb, r, re, yd, yu in zip(threads, g4_mean, g4_err, opt_mean, opt_err, ratio, rel_err, yerr_down, yerr_up):
    line = f"{t:6d}  {a:7.3f}  {ea:6.3f}  {b:8.3f}  {eb:7.3f}  {r:5.3f}  {re:7.4f}  {yd:9.3f}  {yu:8.3f}"
    print(line)
    lines.append(line)
with open("ratio_stats.txt", "w") as outf:
    outf.write("\n".join(lines) + "\n")

# --- Plot: points with vertical error bars only, log y-axis, integer x ticks, no title ---
threads_arr = np.array(threads, dtype=int)
valid = np.isfinite(ratio) & np.isfinite(yerr_down) & np.isfinite(yerr_up) & (ratio > 0)
x = threads_arr[valid]
y = ratio[valid]
yerr_plot = yerr[:, valid]

fig, ax = plt.subplots(figsize=(8, 5))
ax.errorbar(x, y, yerr=yerr_plot, fmt='o', linestyle='none', capsize=3)
ax.set_yscale('log')
ax.set_xlabel('Number of G4 threads')
ax.set_ylabel('G4 time / Opticks time')
ax.grid(True, which='both', alpha=0.3)

# Integer x-axis tick labels
ax.set_xticks(x.tolist())
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

fig.tight_layout()
fig.savefig('g4_opticks_ratio_log.png', dpi=200)
print("Plot saved as g4_opticks_ratio_log.png; table also written to ratio_stats.txt")
