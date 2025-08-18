import subprocess
import re
from pathlib import Path

opticks_file = "Opticks.txt"

run_mac_template = """
/run/verbose 1
/process/optical/cerenkov/setStackPhotons true
/run/initialize
/run/beamOn {events}
"""

def parse_sim_time(output: str):
    """
    Parse lines like: 'Simulation time: 12.345 seconds'
    from combined stdout+stderr.
    """
    m = re.search(r"Simulation time:\s*([\d.]+)\s*seconds", output)
    return float(m.group(1)) if m else None

def parse_num_collected(output: str):
    """
    Parse lines like: 'Opticks: NumCollected:  27961029'
    from combined stdout+stderr.
    """
    m = re.search(r"Opticks:\s*NumCollected:\s*([0-9,]+)", output)
    if not m:
        return None
    return int(m.group(1).replace(",", ""))

def logspace_events(start=1, stop=50000):
    """
    Generate a log-spaced integer sequence using the common 1–2–5 per decade pattern:
    1, 2, 5, 10, 20, 50, 100, ..., up to <= stop.
    Ensures 'start' and 'stop' are included when within range.
    """
    vals = []
    for exp in range(0, 6):  # 10^0 through 10^5
        base = 10 ** exp
        for m in (1, 2, 5):
            n = m * base
            if start <= n <= stop:
                vals.append(n)
    # Ensure uniqueness and ascending order, and include endpoints if missing
    vals = sorted(set(vals))
    if start < vals[0]:
        vals.insert(0, start)
    if stop > vals[-1]:
        vals.append(stop)
    return vals

def write_run_mac(n_events: int, path: str = "run.mac"):
    with open(path, "w") as rm:
        rm.write(run_mac_template.format(events=n_events))

def main():
    gdml_path = "esi-g4ox/geom/pfrich_min_FINAL.gdml"
    exe = "./build/src/simg4ox"

    events_list = logspace_events(1, 100000)

    with open(opticks_file, "w") as of:
        for n in events_list:
            write_run_mac(n)

            # Run the simulation (no external 'time' wrapper)
            result = subprocess.run(
                [exe, "-g", gdml_path, "-m", "run.mac"],
                capture_output=True, text=True
            )

            combined = (result.stdout or "") + "\n" + (result.stderr or "")

            sim_time = parse_sim_time(combined)
            num_collected = parse_num_collected(combined)

            if sim_time is None:
                print(f"[!] Could not parse Simulation time for beamOn={n}")
                continue
            if not num_collected or num_collected <= 0:
                print(f"[!] Could not parse valid NumCollected for beamOn={n}")
                continue

            per_collected = sim_time / num_collected
            # Save: n_events  (sim_time / NumCollected)
            of.write(f"{n} {per_collected:.9e}\n")
            of.flush()
            print(f"[ok] beamOn={n:6d}  sim_time={sim_time:.4f}s  NumCollected={num_collected}  per_collected={per_collected:.9e}s")

    print("Done.")

if __name__ == "__main__":
    main()
