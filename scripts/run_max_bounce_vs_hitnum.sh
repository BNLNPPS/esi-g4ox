#!/usr/bin/env bash
# Sweep OPTICKS_MAX_BOUNCE, repeat each setting 10x, and plot mean±stdev as points.

set -u
MACRO=esi-g4ox/run.mac
EXE=./build/src/simg4ox
GDML=esi-g4ox/geom/pfrich_min_FINAL.gdml

RESULTS=hits_vs_bounce.txt            # raw per-run (bounce, hits)
AGG_CSV=hits_vs_bounce_aggregated.csv # arrays used for plotting
AGG_NPZ=hits_vs_bounce_aggregated.npz # arrays used for plotting
PLOT=hits_vs_bounce.png
LOGDIR=logs
REPEATS=10

mkdir -p "$LOGDIR"

echo "#OPTICKS_MAX_BOUNCE Hits  (raw; each row is one repetition)" > "$RESULTS"

# Extract the last NumHits from a log (returns nonzero on failure)
extract_numhits () {
  local log="$1"
  awk '/Opticks:[[:space:]]*NumHits:/ {val=$NF} END {if (val=="") exit 1; print val}' "$log"
}

for bounce in $(seq 5 31); do
  for rep in $(seq 1 $REPEATS); do
    attempt=1
    while : ; do
      export OPTICKS_MAX_BOUNCE="$bounce"
      log="$LOGDIR/run_bounce_${bounce}_rep${rep}_attempt${attempt}.log"

      echo "[*] Running OPTICKS_MAX_BOUNCE=$OPTICKS_MAX_BOUNCE  rep=$rep/$REPEATS  attempt=$attempt"

      set -o pipefail
      timeout 1800 "$EXE" -g "$GDML" -m "$MACRO" 2>&1 | tee "$log"
      exitcode=${PIPESTATUS[0]}
      set +o pipefail

      if [[ $exitcode -eq 0 ]]; then
        if hits=$(extract_numhits "$log"); then
          printf "%d %d\n" "$bounce" "$hits" | tee -a "$RESULTS"
        else
          echo "WARNING: No 'Opticks: NumHits:' in $log — recording 0."
          printf "%d %d\n" "$bounce" 0 | tee -a "$RESULTS"
        fi
        break
      else
        echo "Run failed (exit $exitcode); retrying…"
        attempt=$((attempt+1))
      fi
    done
  done
done

# ---------- Plot (points + error bars). Also save plotted arrays to CSV/NPZ ----------
python3 - <<'PY'
import numpy as np, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAW="hits_vs_bounce.txt"
OUT_PNG="hits_vs_bounce.png"
OUT_CSV="hits_vs_bounce_aggregated.csv"
OUT_NPZ="hits_vs_bounce_aggregated.npz"

# Load raw (bounce, hits), ignoring comments/blank lines
xs, ys = [], []
with open(RAW) as f:
    for line in f:
        s=line.strip()
        if not s or s.startswith("#"): continue
        b, h = s.split()[:2]
        xs.append(int(b)); ys.append(float(h))
xs = np.asarray(xs); ys = np.asarray(ys)

if xs.size == 0:
    print("[!] No data to plot.", file=sys.stderr); sys.exit(1)

# Aggregate by bounce: mean and sample stdev (ddof=1)
uniq = np.unique(xs)
mean = np.zeros_like(uniq, dtype=float)
stdev = np.zeros_like(uniq, dtype=float)
count = np.zeros_like(uniq, dtype=int)
for i, b in enumerate(uniq):
    v = ys[xs==b]
    count[i] = v.size
    mean[i]  = np.mean(v) if v.size else np.nan
    stdev[i] = np.std(v, ddof=1) if v.size>1 else 0.0

# Save the arrays that are plotted
hdr = "OPTICKS_MAX_BOUNCE,repeats,mean_hits,std_hits"
np.savetxt(OUT_CSV, np.column_stack([uniq, count, mean, stdev]),
           delimiter=",", header=hdr, comments="", fmt="%.10g")
np.savez(OUT_NPZ, bounce=uniq, repeats=count, mean=mean, stdev=stdev)

# Plot: points with vertical error bars, no connecting line, no title
plt.figure(figsize=(7,5))
plt.errorbar(uniq, mean, yerr=stdev, fmt="o", linestyle="none", capsize=3)
plt.xlabel("OPTICKS_MAX_BOUNCE")
plt.ylabel("Opticks NumHits")
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=220)
print(f"Saved {OUT_PNG}")
print(f"Saved plotted data: {OUT_CSV}, {OUT_NPZ}")
PY

echo "Done. Raw: $RESULTS  |  Aggregated: $AGG_CSV $AGG_NPZ  |  Plot: $PLOT  |  Logs: $LOGDIR/"
