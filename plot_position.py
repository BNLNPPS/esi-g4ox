#!/usr/bin/env python3

import re
import matplotlib.pyplot as plt

def parse_g4_line(line):
    """
    Example line:
      Adding hit from Geant4: 5.44596 eV  (-398.649, -199.208, 165.5)  (...)
    """
    # Find all parenthetical groups: e.g. ["-398.649, -199.208, 165.5", ...]
    groups = re.findall(r"\((.*?)\)", line)
    if not groups:
        return None

    # The first group is the position
    position_str = groups[0]  # e.g. "-398.649, -199.208, 165.5"
    coords = position_str.split(",")  # ["-398.649", " -199.208", " 165.5"]
    if len(coords) < 2:
        return None

    # Convert x, y to float
    x = float(coords[0])
    y = float(coords[1])
    return x, y

def parse_opticks_line(line):
    """
    Example line:
      190.138  (-399.023, -217.986, 165.5)  (0.00282741, 0.216947, 0.976179) ...
    """
    groups = re.findall(r"\((.*?)\)", line)
    if not groups:
        return None

    # The first group is the position
    position_str = groups[0]  # e.g. "-399.023, -217.986, 165.5"
    coords = position_str.split(",")
    if len(coords) < 2:
        return None

    x = float(coords[0])
    y = float(coords[1])
    return x, y

def main():
    # Filenames (adjust as needed)
    g4_file = "g4_photon_hits.txt"
    opticks_file = "opticks_hits_output.txt"

    g4_x, g4_y = [], []
    opticks_x, opticks_y = [], []

    # --- Read Geant4 hits ---
    with open(g4_file, "r") as f:
        for line in f:
            parsed = parse_g4_line(line.strip())
            if parsed is not None:
                x, y = parsed
                g4_x.append(x)
                g4_y.append(y)

    # --- Read Opticks hits ---
    with open(opticks_file, "r") as f:
        for line in f:
            parsed = parse_opticks_line(line.strip())
            if parsed is not None:
                x, y = parsed
                opticks_x.append(x)
                opticks_y.append(y)

    # --- Plot the data ---
    plt.figure(figsize=(8, 6))
    plt.scatter(g4_x, g4_y, color='red', marker='o', label='Geant4 hits')
    plt.scatter(opticks_x, opticks_y, color='blue', marker='x', label='Opticks hits')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.title('Geant4 vs. Opticks Hit Positions')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
