#!/usr/bin/env python3
"""
Build XY and XZ 2D histograms (10×10, log-scale) for Geant4 vs Opticks hits.
...
"""

import os, re, math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

BASE_FONTSIZE = plt.rcParams.get("font.size", 10.0)
plt.rcParams.update({
    "font.size"      : BASE_FONTSIZE * 2.0,   # base text
    "axes.titlesize" : BASE_FONTSIZE * 2.4,   # subplot titles
    "axes.labelsize" : BASE_FONTSIZE * 2.2,   # x/y labels
    "xtick.labelsize": BASE_FONTSIZE * 2.0,   # tick numbers
    "ytick.labelsize": BASE_FONTSIZE * 2.0,
    "legend.fontsize": BASE_FONTSIZE * 2.0,
})

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
        axs[0].set_title("Geant4")
        axs[0].set_xlabel(xlabel); axs[0].set_ylabel(ylabel)

        h2 = axs[1].hist2d(x2, y2, bins=100, norm=LogNorm(), cmap="viridis")
        axs[1].set_title("EIC-Opticks")
        axs[1].set_xlabel(xlabel)

        for ax in axs:
            ax.tick_params(axis="both", which="both",
                           labelsize=plt.rcParams["xtick.labelsize"])

        cbar = fig.colorbar(h1[3], ax=axs.ravel().tolist(), pad=0.02)
        cbar.set_label("Hits per bin",
                       fontsize=plt.rcParams["axes.labelsize"])         
        cbar.ax.tick_params(labelsize=plt.rcParams["xtick.labelsize"])   

        if title:
            fig.suptitle(
                title,
                fontsize=plt.rcParams["axes.titlesize"] * 1.2
            )

        plt.savefig(outfile, dpi=300)
        plt.close()
        print(f"Saved {outfile}")

    # --- Produce plots ---------------------------------------------------
    make_hist("",
              g4_x, g4_y, ok_x, ok_y,
              "X position [mm]", "Y position [mm]", xy_out)

    make_hist("",
              g4_x, g4_z, ok_x, ok_z,
              "X position [mm]", "Z position [mm]", xz_out)

if __name__ == "__main__":
    main()
