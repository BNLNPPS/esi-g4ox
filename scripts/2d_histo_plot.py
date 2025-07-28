#!/usr/bin/env python3
"""
Build XY and XZ 2D histograms (10×10, log‑scale) for Geant4 vs Opticks hits.

Geant4 hits are collected from consecutive files
      g4_photon_hits_thread0.txt, g4_photon_hits_thread1.txt, ...
  and filtered so that only the sensitive layer is hit
Opticks hits come from  opticks_hits_output.txt  (all kept).
Images saved:  hits_xy_hist.png  and  hits_xz_hist.png
Prints:
      - # of Geant4 hits kept
      - # of Geant4 hits discarded by the Z cut 
      - # of Opticks hits
      - Fractional Poisson error
"""

import os, re, math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ---- Configurable Z cut ------------------------------------------------
Z_MIN, Z_MAX = 187.125, 187.525   # mm

# ---- Helpers -----------------------------------------------------------
def parse_line(line):
    m = re.findall(r"\(([^)]+)\)", line)
    if not m:
        return None
    parts = m[0].split(",")
    if len(parts) < 3:
        return None
    return tuple(float(p) for p in parts[:3])

def read_positions(fname):
    x, y, z = [], [], []
    with open(fname) as f:
        for line in f:
            p = parse_line(line)
            if p:
                xi, yi, zi = p
                x.append(xi); y.append(yi); z.append(zi)
    return x, y, z

def gather_g4_hits(prefix="g4_photon_hits_thread"):
    kept_x, kept_y, kept_z = [], [], []

    reject_count = 0
    idx = 0
    while True:
        fname = f"{prefix}{idx}.txt"
        if not os.path.isfile(fname):
            break
        xi, yi, zi = read_positions(fname)

        # Apply Z window (inclusive)
        for X, Y, Z in zip(xi, yi, zi):
            if Z_MIN <= Z <= Z_MAX:
                kept_x.append(X); kept_y.append(Y); kept_z.append(Z)
            else:
                reject_count += 1
        idx += 1

    if idx == 0:
        raise FileNotFoundError("No g4_photon_hits_thread*.txt files found.")
    return kept_x, kept_y, kept_z, reject_count

# ---- Main --------------------------------------------------------------
def main():
    xy_out, xz_out = "hits_xy_hist.png", "hits_xz_hist.png"
    opticks_file   = "opticks_hits_output.txt"

    # --- Load & filter ---------------------------------------------------
    g4_x, g4_y, g4_z, g4_rejected = gather_g4_hits()
    ok_x, ok_y, ok_z              = read_positions(opticks_file)

    # --- Stats -----------------------------------------------------------
    N_g4 = len(g4_x)
    N_ok = len(ok_x)

    frac_err_g4 = math.sqrt(N_g4)/N_g4 if N_g4 else float("nan")
    frac_err_ok = math.sqrt(N_ok)/N_ok if N_ok else float("nan")

    # --- Print summary ---------------------------------------------------
    print("=== Hit statistics (after Z cut applied to Geant4) ===")
    print(f"Geant4 kept hits   : {N_g4:>10} ")
    print(f"Geant4 rejected Zs : {g4_rejected:>10}")
    print(f"Opticks hits       : {N_ok:>10} \n")

    # --- Plot helper -----------------------------------------------------
    def make_hist(title, x1, y1, x2, y2, xlabel, ylabel, outfile):
        fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True,
                                constrained_layout=True)

        h1 = axs[0].hist2d(x1, y1, bins=100, norm=LogNorm(), cmap="viridis")
        axs[0].set_title("Geant4 (filtered)")
        axs[0].set_xlabel(xlabel); axs[0].set_ylabel(ylabel)

        h2 = axs[1].hist2d(x2, y2, bins=100, norm=LogNorm(), cmap="viridis")
        axs[1].set_title("Opticks")
        axs[1].set_xlabel(xlabel)

        cbar = fig.colorbar(h1[3], ax=axs.ravel().tolist(), pad=0.02)
        cbar.set_label("Hits per bin (log scale)")

        fig.suptitle(title, fontsize=14)
        plt.savefig(outfile, dpi=300)
        plt.close()
        print(f"Saved {outfile}")

    # --- Produce plots ---------------------------------------------------
    make_hist("Photon hits XY plane (Z filtered Geant4)",
              g4_x, g4_y, ok_x, ok_y,
              "X position", "Y position", xy_out)

    make_hist("Photon hits XZ plane (Z filtered Geant4)",
              g4_x, g4_z, ok_x, ok_z,
              "X position", "Z position", xz_out)

if __name__ == "__main__":
    main()
