#!/usr/bin/env python3
import re
import matplotlib.pyplot as plt

# ---------- 1. Helpers -------------------------------------------------
def parse_line(line):
    """Extract  floats from a line that contains '(x, y, z)'."""
    m = re.findall(r"\(([^)]+)\)", line)
    if not m:
        return None
    parts = m[0].split(",")
    if len(parts) < 3:
        return None
    return tuple(float(p) for p in parts[:3])

def read_positions(fname, max_hits=100):
    """Return lists x, y, z with at most max_hits entries."""
    x, y, z = [], [], []
    with open(fname) as f:
        for line in f:
            if len(x) >= max_hits:
                break
            parsed = parse_line(line)
            if parsed:
                xi, yi, zi = parsed
                x.append(xi); y.append(yi); z.append(zi)
    return x, y, z

# ---------- 2. Main ----------------------------------------------------
def main():
    g4_file      = "g4_photon_hits_thread0.txt"
    opticks_file = "opticks_hits_output.txt"
    xy_out       = "hits_xy_100.png"
    xz_out       = "hits_xz_100.png"

    g4_x,  g4_y,  g4_z  = read_positions(g4_file,      100)
    ok_x,  ok_y,  ok_z  = read_positions(opticks_file, 100)

    # ---- XY plot ----
    plt.figure(figsize=(8, 6))
    plt.scatter(g4_x, g4_y, c="red",   marker="o", label="Geant4")
    plt.scatter(ok_x, ok_y, c="blue",  marker="x", label="Opticks")
    plt.xlabel("X position")
    plt.ylabel("Y position")
    plt.title("Geant4 vs Opticks XY (first 100 hits)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(xy_out)
    plt.close()
    print(f"XY plot saved to {xy_out}")

    # ---- XZ plot ----
    plt.figure(figsize=(8, 6))
    plt.scatter(g4_x, g4_z, c="red",   marker="o", label="Geant4")
    plt.scatter(ok_x, ok_z, c="blue",  marker="x", label="Opticks")
    plt.xlabel("X position")
    plt.ylabel("Z position")
    plt.title("Geant4 vs Opticks XZ (first 100 hits)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(xz_out)
    plt.close()
    print(f"XZ plot saved to {xz_out}")

if __name__ == "__main__":
    main()
