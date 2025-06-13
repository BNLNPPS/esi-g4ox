#!/bin/bash

MACRO=esi-g4ox/run.mac
EXE=./build/src/simg4ox
GDML=pfrich_min_FINAL.gdml
RESULTS=results.txt

# Backup macro file once
cp "$MACRO" "$MACRO.bak"

echo "#Threads Runtime(s)" > "$RESULTS"

for threads in $(seq 2 40); do
    # Edit the macro file in place
    sed -i "s|/run/numberOfThreads .*|/run/numberOfThreads $threads|" "$MACRO"

    finished=0
    while [ $finished -eq 0 ]; do
        start=$(date +%s)
        timeout 1800 $EXE -g "$GDML" -m "$MACRO"
        exitcode=$?
        end=$(date +%s)
        runtime=$((end - start))

        if [ $exitcode -eq 0 ]; then
            echo "$threads $runtime" | tee -a "$RESULTS"
            finished=1
        else
            echo "Run with $threads threads timed out or failed (exit code $exitcode), retrying..."
        fi
    done
done

# Restore the original macro file
mv "$MACRO.bak" "$MACRO"
