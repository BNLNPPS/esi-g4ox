import numpy as np
import subprocess
import matplotlib.pyplot as plt
import re
import os

# Generate 20 logarithmically spaced t_min values between 0.00001 and 1
values = np.logspace(-5, 0, 50)
results = []

for value in values:
    env = os.environ.copy()
    env["OPTICKS_PROPAGATE_EPSILON"] = str(value)

    try:
        # Run the simulation
        proc = subprocess.run(
            ["./build/src/simg4ox", "-g", "esi-g4ox/geom/sphere_leak.gdml", "-m", "esi-g4ox/run.mac"],
            env=env,
            capture_output=True,
            text=True,
            check=True
        )

        # Extract NumHits from stdout
        match = re.search(r'Opticks:\s+NumHits:\s+(\d+)', proc.stdout)
        if match:
            num_hits = int(match.group(1))
            results.append(num_hits / 60.0)  # Convert to ppm
        else:
            print(f"NumHits not found for t_min = {value}")
            results.append(0.0)
    except subprocess.CalledProcessError as e:
        print(f"Error running simulation for t_min = {value}")
        print(e.stderr)
        results.append(0.0)

# Save results to CSV
with open("leak_results.csv", "w") as f:
    f.write("t_min,leaked_photon_ppm\n")
    for v, r in zip(values, results):
        f.write(f"{v},{r}\n")

# Plotting
plt.figure()
plt.plot(values, results, marker='o')
plt.xscale('log')
plt.xlabel("X: t_min")
plt.ylabel("Y: leaked photon number [ppm]")
plt.title("Photon Leak vs t_min")
plt.grid(True)
plt.savefig("leak_plot.png")
plt.show()

