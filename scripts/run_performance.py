import subprocess
import re
import shutil

opticks_file = "Opticks.txt"

run_mac_template = """
/run/numberOfThreads {threads}
/run/verbose 1
/process/optical/cerenkov/setStackPhotons {flag}
/run/initialize
/run/beamOn 50000
"""

def parse_sim_time(output):
    m = re.search(r"Simulation time:\s*([\d.]+)\s*seconds", output)
    return float(m.group(1)) if m else None

def run_and_time(cmd):
    """
    Runs cmd and returns (stdout, stderr, real_seconds).
    Uses /usr/bin/time -f "%e" when available; falls back to parsing bash 'time'.
    """
    time_path = shutil.which("/usr/bin/time")
    if time_path:
        # time prints seconds to stderr as a single float with -f "%e"
        result = subprocess.run([time_path, "-f", "%e"] + cmd,
                                capture_output=True, text=True)
        stdout, stderr = result.stdout, result.stderr
        # last token on stderr should be the elapsed seconds
        try:
            real_sec = float(stderr.strip().splitlines()[-1])
        except Exception:
            real_sec = None
        return stdout, stderr, real_sec
    else:
        # Fallback: bash built-in 'time' outputs 'real 0mX.XXXs'
        result = subprocess.run(["bash", "-c", "time " + " ".join(cmd)],
                                capture_output=True, text=True)
        stdout, stderr = result.stdout, result.stderr
        m = re.search(r"real\s+(\d+)m([\d.]+)s", stderr)
        real_sec = int(m.group(1))*60 + float(m.group(2)) if m else None
        return stdout, stderr, real_sec

# Repeat the full sweep 10 times; write (true - false) to timings{run}.txt
with open(opticks_file, "a") as of:  # append across runs
    for run_idx in range(10):
        timings_file = f"timings{run_idx}.txt"
        with open(timings_file, "w") as tf:
            for threads in range(1, 21):
                times = {}

                for flag in ["true", "false"]:
                    # Prepare run.mac for this configuration
                    with open("run.mac", "w") as rm:
                        rm.write(run_mac_template.format(threads=threads, flag=flag))

                    cmd = ["./build/src/simg4ox",
                           "-g", "esi-g4ox/geom/pfrich_min_FINAL.gdml",
                           "-m", "run.mac"]

                    stdout, stderr, real_sec = run_and_time(cmd)

                    if real_sec is None:
                        print(f"[!] Could not parse wall time for run={run_idx} threads={threads} flag={flag}")
                    else:
                        times[flag] = real_sec

                    # Record Opticks simulation time (from program output) for 'true'
                    if flag == "true":
                        sim_time_true = parse_sim_time(stdout + stderr)
                        if sim_time_true is not None:
                            of.write(f"{threads} {sim_time_true}\n")
                            of.flush()

                # Write Geant4 result = (true - false)
                if "true" in times and "false" in times:
                    diff = times["true"] - times["false"]
                    tf.write(f"{threads} {diff}\n")
                    tf.flush()
                else:
                    print(f"[!] Missing times to compute diff for run={run_idx} threads={threads}")

print("Done.")
