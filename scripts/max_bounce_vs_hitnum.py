#!/usr/bin/env bash
# Sweep OPTICKS_MAX_BOUNCE and record Opticks NumHits parsed from stdout.

set -u
MACRO=esi-g4ox/run.mac
EXE=./build/src/simg4ox
GDML=esi-g4ox/geom/pfrich_min_FINAL.gdml


RESULTS=hits_vs_bounce.txt
PLOT=hits_vs_bounce.png
LOGDIR=logs

mkdir -p "$LOGDIR"

echo "#OPTICKS_MAX_BOUNCE Hits" > "$RESULTS"

# Function to extract the last NumHits from a log
extract_numhits () {
  local log="$1"
  # Grab the last occurrence of "Opticks: NumHits:  <int>" and print the int
  awk '/Opticks:[[:space:]]*NumHits:/ {val=$NF} END {if (val=="") exit 1; print val}' "$log"
}

for bounce in $(seq 5 31); do
  finished=0
  while [[ $finished -eq 0 ]]; do
    export OPTICKS_MAX_BOUNCE="$bounce"
    log="$LOGDIR/run_bounce_${bounce}.log"

    echo "[*] Running with OPTICKS_MAX_BOUNCE=$OPTICKS_MAX_BOUNCE"

    # capture both stdout+stderr to log; preserve the real exit code of the sim
    set -o pipefail
    timeout 1800 "$EXE" -g "$GDML" -m "$MACRO" 2>&1 | tee "$log"
    exitcode=${PIPESTATUS[0]}
    set +o pipefail

    if [[ $exitcode -eq 0 ]]; then
      # Parse NumHits from stdout log
      if hits=$(extract_numhits "$log"); then
        printf "%d %d\n" "$bounce" "$hits" | tee -a "$RESULTS"
        finished=1
      else
        echo "WARNING: Could not find 'Opticks: NumHits:' in $log — recording 0 and continuing."
        printf "%d %d\n" "$bounce" 0 | tee -a "$RESULTS"
        finished=1
      fi
    else
      echo "Run with OPTICKS_MAX_BOUNCE=$bounce failed or timed out (exit $exitcode), retrying..."
    fi
  done
done

# --- Plot relation (line + markers) --------------------------------------
python3 - <<'PY'
import matplotlib.pyplot as plt

xs, ys = [], []
with open("hits_vs_bounce.txt") as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        b, h = line.split()[:2]
        xs.append(int(b)); ys.append(int(h))

plt.figure(figsize=(7,5))
plt.plot(xs, ys, marker='o')
plt.xlabel("OPTICKS_MAX_BOUNCE")
plt.ylabel("Opticks NumHits")
plt.title("Opticks NumHits vs. OPTICKS_MAX_BOUNCE")
plt.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig("hits_vs_bounce.png", dpi=220)
print("Saved hits_vs_bounce.png")
PY

echo "Done. Data: $RESULTS  |  Plot: $PLOT  |  Logs: $LOGDIR/"
