import numpy as np
import matplotlib.pyplot as plt

# Read Geant4 real times
g4_threads = []
g4_times = []
with open('timings.txt', 'r') as f:
    for line in f:
        t, val = line.strip().split()
        g4_threads.append(int(t))
        g4_times.append(float(val))

# Calculate average, skipping the first entry since that includes the geometry upload
def compute_average(filename):
    values = []
    with open(filename, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    values.append(float(parts[1]))
                except ValueError:
                    continue
    if len(values) <= 1:
        print(f"Not enough values in {filename} to calculate average (excluding first entry).")
        return
    # Exclude the first entry
    avg = sum(values[1:]) / len(values[1:])
    print(f"Average (excluding first entry) for {filename}: {avg:.3f}")
    return avg

# Calculate average Opticks time
opticks_avg = compute_average("Opticks.txt")

# Calculate G4/Opticks ratio for each thread
ratios = [g4_time / opticks_avg for g4_time in g4_times]

# Plot and save
plt.figure(figsize=(8, 5))
plt.plot(g4_threads, ratios, marker='o')
plt.xlabel('Number of G4 threads')
plt.ylabel('G4 simulation time / Opticks simulation time')
plt.title('G4 vs Opticks Simulation Time Scaling')
plt.grid(True)
plt.tight_layout()
plt.savefig('g4_opticks_ratio.png', dpi=200)
print("Plot saved as g4_opticks_ratio.png")
