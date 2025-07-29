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

# Read Opticks times
opticks_times = []
with open('Opticks.txt', 'r') as f:
    for line in f:
        _, val = line.strip().split()
        opticks_times.append(float(val))

# Calculate average Opticks time
opticks_avg = np.mean(opticks_times)

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
