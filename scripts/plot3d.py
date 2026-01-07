#!/usr/bin/env python3

import re
import matplotlib.pyplot as plt

def parse_g4_line(line):
    """
    Example line:
      Adding hit from Geant4: 5.44596 eV  (-398.649, -199.208, 165.5)  (...)
    """
    groups = re.findall(r"\((.*?)\)", line)
    if not groups:
        return None

    position_str = groups[0]
    coords = position_str.split(",")
    if len(coords) < 3:
        return None

    x = float(coords[0])
    y = float(coords[1])
    z = float(coords[2])
    return x, y, z

def parse_opticks_line(line):
    """
    Example line:
      190.138  (-399.023, -217.986, 165.5)  (0.00282741, 0.216947, 0.976179) ...
    """
    groups = re.findall(r"\((.*?)\)", line)
    if not groups:
        return None

    position_str = groups[0]
    coords = position_str.split(",")
    if len(coords) < 3:
        return None

    x = float(coords[0])
    y = float(coords[1])
    z = float(coords[2])
    return x, y, z

def main():
    g4_file = "g4.txt"
    opticks_file = "o.txt"

    g4_x, g4_y, g4_z = [], [], []
    opticks_x, opticks_y, opticks_z = [], [], []

    # --- Read Geant4 hits ---
    with open(g4_file, "r") as f:
        for line in f:
            parsed = parse_g4_line(line.strip())
            if parsed is not None:
                x, y, z = parsed
                g4_x.append(x)
                g4_y.append(y)
                g4_z.append(z)

    # --- Read Opticks hits ---
    with open(opticks_file, "r") as f:
        for line in f:
            parsed = parse_opticks_line(line.strip())
            if parsed is not None:
                x, y, z = parsed
                opticks_x.append(x)
                opticks_y.append(y)
                opticks_z.append(z)

    # --- Plot the data in 3D ---
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(g4_x, g4_y, g4_z, color='red', marker='o', label='Geant4 hits')
    ax.scatter(opticks_x, opticks_y, opticks_z, color='blue', marker='x', label='Opticks hits')

    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.set_title('Geant4 vs. Opticks Hit Positions (3D)')
    ax.legend()
    ax.grid(True)

    plt.show()

if __name__ == "__main__":
    main()
