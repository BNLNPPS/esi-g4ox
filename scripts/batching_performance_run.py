#!/usr/bin/env python3
import subprocess
import re
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless backend for saving PNGs without a display
import matplotlib.pyplot as plt

# Raw per-run output (two columns: events, time_per_photon_sec)
OPTICKS_OUT = "Opticks.txt"

# Aggregated + plotted outputs
CSV_OUT = "Opticks_aggregated.csv"
NPZ_OUT = "Opticks_aggregated.npz"
PLOT_TIME = "time_per_photon_opticks.png"
PLOT_SPEED = "speedup_batching.png"

REPEATS = 10

RUN_MAC_TEMPLATE = """
/run/verbose 1
/process/optical/cerenkov/setStackPhotons true
/run/initialize
/run/beamOn {events}
"""

def parse_sim_time(output: str):
    """Parse: 'Simulation time: 12.345 seconds' from combined stdout+stderr."""
    m = re.search(r"Simulation time:\s*([\d.]+)\s*seconds", output)
    return float(m.group(1)) if m else None

def _find_all_num_collected(text: str):
    """Return all NumCollected integers found, in order of appearance."""
    vals = re.findall(r"Opticks:\s*NumCollected:\s*([0-9,]+)", text)
    return [int(v.replace(",", "")) for v in vals]

def parse_num_collected_preferred(output: str):
    """
    Prefer the last 'Opticks: NumCollected:' that appears BEFORE 'Opticks: NumHits:'.
    Fallback to the last 'NumCollected' anywhere if needed.
    Return (num_collected:int or None, all_matches:list[int]).
    """
    head, sep, _ = output.partition("Opticks: NumHits:")
    candidates = _find_all_num_collected(head if sep else output)
    if candidates:
        return candidates[-1], candidates
    all_matches = _find_all_num_collected(output)
    if all_matches:
        return all_matches[-1], all_matches
    return None, []

def logspace_events(start=1, stop=100000):
    """1–2–5 per decade pattern up to <= stop, inclusive endpoints if missing."""
    vals = []
    for exp in range(0, 6):  # 10^0 through 10^5
        base = 10 ** exp
        for m in (1, 2, 5):
            n = m * base
            if start <= n <= stop:
                vals.append(n)
    vals = sorted(set(vals))
    if vals and start < vals[0]:
        vals.insert(0, start)
    if vals and stop > vals[-1]:
        vals.append(stop)
    return vals or [start, stop]

def write_run_mac(n_events: int, path: str = "run.mac"):
    Path(path).write_text(RUN_MAC_TEMPLATE.format(events=n_events), encoding="utf-8")

def aggregate(recs):
    """
    recs: list of (events:int, per_collected:float) for all runs
    Returns arrays sorted by events: x, mean, stdev, counts
    """
    # group values by X
    by_x = {}
    for x, v in recs:
        by_x.setdefault(int(x), []).append(float(v))

    xs = np.array(sorted(by_x.keys()), dtype=float)
    means = np.zeros_like(xs, dtype=float)
    stds = np.zeros_like(xs, dtype=float)
    counts = np.zeros_like(xs, dtype=int)

    for i, x in enumerate(xs):
        arr = np.array(by_x[int(x)], dtype=float)
        counts[i] = arr.size
        means[i] = arr.mean() if arr.size else np.nan
        stds[i] = arr.std(ddof=1) if arr.size > 1 else 0.0  # sample stdev

    # warn if repeats != REPEATS
    bad = counts != REPEATS
    if np.any(bad):
        details = ", ".join(f"{int(xs[i])}→{counts[i]}" for i in np.where(bad)[0])
        print(f"[!] Warning: repeats per events differ from expected {REPEATS}: {details}", file=sys.stderr)
    return xs, means, stds, counts

def save_plotted_data_csv(x, y_mean, y_std, y_mean_us, y_std_us, speed_mean, speed_std, counts, path):
    header = (
        "events,mean_sec_per_photon,std_sec_per_photon,mean_us_per_photon,std_us_per_photon,"
        "speedup_mean,speedup_std,repeats"
    )
    arr = np.column_stack([x, y_mean, y_std, y_mean_us, y_std_us, speed_mean, speed_std, counts])
    np.savetxt(path, arr, delimiter=",", header=header, comments="", fmt="%.10g")

def main():
    gdml_path = "esi-g4ox/geom/pfrich_min_FINAL.gdml"
    exe = "./build/src/simg4ox"

    events_list = logspace_events(1, 100000)
    raw_records = []  # (events, per_collected_sec) for every repetition

    with open(OPTICKS_OUT, "w") as of:
        for n in events_list:
            write_run_mac(n)

            for r in range(REPEATS):
                result = subprocess.run(
                    [exe, "-g", gdml_path, "-m", "run.mac"],
                    capture_output=True, text=True
                )
                combined = (result.stdout or "") + "\n" + (result.stderr or "")

                sim_time = parse_sim_time(combined)
                num_collected, all_matches = parse_num_collected_preferred(combined)

                if sim_time is None:
                    print(f"[!] Could not parse Simulation time for beamOn={n} (rep {r+1}/{REPEATS})")
                    continue
                if not num_collected or num_collected <= 0:
                    print(f"[!] Could not parse valid NumCollected for beamOn={n} (rep {r+1}/{REPEATS}). "
                          f"Found matches={all_matches}")
                    continue

                per_collected = sim_time / num_collected
                raw_records.append((n, per_collected))

                # Keep backward-compatible raw output: two columns (events, per_collected_seconds)
                of.write(f"{n} {per_collected:.9e}\n")
                of.flush()

                print(f"[ok] beamOn={n:6d} rep={r+1:02d}/{REPEATS}  sim_time={sim_time:.4f}s  "
                      f"NumCollected={num_collected:9d}  per_collected={per_collected:.9e}s")

    # ---- Aggregate across repeats ----
    if not raw_records:
        print("[!] No valid records to aggregate.", file=sys.stderr)
        print("Done.")
        return

    x, y_mean, y_std, counts = aggregate(raw_records)

    # ---- Prepare data for Plot 1 (time per photon, µs) ----
    y_mean_us = y_mean * 1e6
    y_std_us = y_std * 1e6

    # ---- Plot 1: points with vertical error bars (no lines) ----
    plt.figure(figsize=(7.5, 4.8))
    plt.errorbar(x, y_mean_us, yerr=y_std_us, fmt="o", linestyle="none", capsize=3)
    plt.xlabel("number of Geant4 Events in a single GPU call")
    plt.ylabel("Simulation time per photon [µs]")
    plt.title("Per-photon simulation time vs. batch size (mean ± stdev)")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(PLOT_TIME, dpi=200)
    plt.close()

    # ---- Plot 2: speedup relative to first mean ----
    if y_mean[0] <= 0:
        print("[!] First mean Y must be positive to compute speedup.", file=sys.stderr)
        speedup_mean = np.full_like(y_mean, np.nan)
        speedup_std = np.full_like(y_mean, np.nan)
    else:
        b = y_mean[0]
        sb = y_std[0]
        speedup_mean = b / y_mean
        # First-order error propagation for ratio S=b/m:
        # Var(S) ≈ (S^2)*[(sb/b)^2 + (s/m)^2]
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_var = (sb / b)**2 + (y_std / y_mean)**2
        rel_var = np.nan_to_num(rel_var, nan=np.inf, posinf=np.inf, neginf=np.inf)
        speedup_std = speedup_mean * np.sqrt(rel_var)

    plt.figure(figsize=(7.5, 4.8))
    plt.errorbar(x, speedup_mean, yerr=speedup_std, fmt="o", linestyle="none", capsize=3)
    plt.xlabel("Number of Geant4 Events")
    plt.ylabel("Speedup due to batching")
    plt.title("Batching speedup (relative to simulating one event per GPU call)")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(PLOT_SPEED, dpi=200)
    plt.close()

    # ---- Save the arrays that were plotted (for later plotting) ----
    save_plotted_data_csv(
        x, y_mean, y_std,
        y_mean_us, y_std_us,
        speedup_mean, speedup_std,
        counts, CSV_OUT
    )
    np.savez(
        NPZ_OUT,
        x=x,
        repeats=counts,
        y_mean=y_mean,
        y_std=y_std,
        y_mean_us=y_mean_us,
        y_std_us=y_std_us,
        speedup_mean=speedup_mean,
        speedup_std=speedup_std,
    )

    print(f"[ok] Raw per-run data: {OPTICKS_OUT}")
    print(f"[ok] Aggregated CSV:   {CSV_OUT}")
    print(f"[ok] Aggregated NPZ:   {NPZ_OUT}")
    print(f"[ok] Saved plots:      {PLOT_TIME}, {PLOT_SPEED}")
    print("Done.")

if __name__ == "__main__":
    main()
