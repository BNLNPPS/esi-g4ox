#!/usr/bin/env python3
import subprocess
import re
from pathlib import Path

# Output file for (n_events, sim_time / NumCollected)
OPTICKS_OUT = "Opticks.txt"

RUN_MAC_TEMPLATE = """
/run/verbose 1
/process/optical/cerenkov/setStackPhotons true
/run/initialize
/run/beamOn {events}
"""

def parse_sim_time(output: str):
    """
    Parse: 'Simulation time: 12.345 seconds' from combined stdout+stderr.
    """
    m = re.search(r"Simulation time:\s*([\d.]+)\s*seconds", output)
    return float(m.group(1)) if m else None

def _find_all_num_collected(text: str):
    """
    Return all NumCollected integers found, in order of appearance.
    """
    vals = re.findall(r"Opticks:\s*NumCollected:\s*([0-9,]+)", text)
    return [int(v.replace(",", "")) for v in vals]

def parse_num_collected_preferred(output: str):
    """
    Prefer the last 'Opticks: NumCollected:' that appears BEFORE 'Opticks: NumHits:'.
    Fallback to the last 'NumCollected' anywhere if needed.
    Return (num_collected:int or None, all_matches:list[int]).
    """
    # Prefer matches before 'NumHits'
    head, sep, _ = output.partition("Opticks: NumHits:")
    candidates = _find_all_num_collected(head if sep else output)
    if candidates:
        return candidates[-1], candidates

    # Fallback: look in full output
    all_matches = _find_all_num_collected(output)
    if all_matches:
        return all_matches[-1], all_matches

    return None, []

def logspace_events(start=1, stop=100000):
    """
    1–2–5 per decade pattern up to <= stop, inclusive endpoints if missing.
    """
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

def main():
    gdml_path = "esi-g4ox/geom/pfrich_min_FINAL.gdml"
    exe = "./build/src/simg4ox"

    events_list = logspace_events(1, 100000)

    with open(OPTICKS_OUT, "w") as of:
        for n in events_list:
            write_run_mac(n)

            result = subprocess.run(
                [exe, "-g", gdml_path, "-m", "run.mac"],
                capture_output=True, text=True
            )
            combined = (result.stdout or "") + "\n" + (result.stderr or "")

            sim_time = parse_sim_time(combined)
            num_collected, all_matches = parse_num_collected_preferred(combined)

            if sim_time is None:
                print(f"[!] Could not parse Simulation time for beamOn={n}")
                continue
            if not num_collected or num_collected <= 0:
                print(f"[!] Could not parse valid NumCollected for beamOn={n}. "
                      f"Found matches={all_matches}")
                continue

            per_collected = sim_time / num_collected
            of.write(f"{n} {per_collected:.9e}\n")
            of.flush()

            print(f"[ok] beamOn={n:6d}  sim_time={sim_time:.4f}s  "
                  f"NumCollected={num_collected}  all={all_matches}  "
                  f"per_collected={per_collected:.9e}s")

    print("Done.")

if __name__ == "__main__":
    main()
