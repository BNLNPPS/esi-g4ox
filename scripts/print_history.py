#!/usr/bin/env python3
"""
dump_hit_tracks.py
------------------
Print the full step history for every photon that produced a hit 
(as listed in hits.npy), including its sequence of interactions.
"""

from __future__ import annotations
import argparse
import pathlib as pl
import sys
import numpy as np

# --------------------------------------------------------------------------- #
# Configurable nibble map (edit if your codes differ)                         #
# --------------------------------------------------------------------------- #
LUT = {
    1: "Cherenkov", 2: "Scint", 3: "Miss", 4: "Bulk absorb",
    5: "Bulk reemit", 6: "Bulk scatter", 7: "Surface detect"
}

# --------------------------------------------------------------------------- #
# Sequence decoding (identical to dump_escaped_tracks.py)                     #
# --------------------------------------------------------------------------- #
def decode_sequences(run: pl.Path, n_ph: int) -> np.ndarray:
    """Return ndarray[str] of length n_ph with full interaction strings."""
    nib_file = run / "seqnib.npy"
    if nib_file.is_file():
        arr = np.load(nib_file)
        if arr.ndim == 2:  # (N,16) nibble grid
            return np.array([
                "".join(LUT[n] for n in row if n)
                for row in arr
            ])
        packed = arr.reshape(arr.shape[0], -1)
    else:
        packed = np.load(run / "seq.npy").reshape(n_ph, -1)

    out = []
    for row in packed:
        chars = []
        for w in row:
            w = int(w)
            for _ in range(16):
                nib = w & 0xF
                if nib == 0:
                    w = 0
                    break
                chars.append(LUT.get(nib, "?"))
                w >>= 4
            if w == 0:
                break
        out.append("".join(chars))
    return np.array(out)

# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Dump step history for photons that hit the detector."
    )
    ap.add_argument("run_dir", help="folder with record.npy, hit.npy, and seq*.npy files")
    ap.add_argument("--nmax", type=int, default=0,
                    help="max photons to list (0=all)")
    args = ap.parse_args()

    run = pl.Path(args.run_dir)
    if not (run / "hit.npy").is_file():
        print("Error: hits.npy not found in", run, file=sys.stderr)
        sys.exit(1)

    # 1) Load hits and extract unique photon IDs
    hits = np.load(run / "hit.npy")
    # assume first column is photon index
    ph_ids = np.unique(hits[:, 0].astype(int))
    if ph_ids.size == 0:
        print("No hits found in hits.npy.")
        return

    # 2) Load record data
    rec = np.load(run / "record.npy")             # (N, MAX, 4, 4)
    N, MAX = rec.shape[:2]
    xyz = rec[:, :, 0, :3]                        # (N, MAX, 3)

    # 3) Decode all sequences
    seq_all = decode_sequences(run, N)

    # 4) Print each hit-photon’s track
    to_show = ph_ids if args.nmax == 0 else ph_ids[:args.nmax]
    print(f"{to_show.size} / {N} photons with hits\n")
    for pid in to_show:
        seq = seq_all[pid]
        print(f"Photon {pid}  sequence = {seq or '<empty>'}")
        for step in range(MAX):
            pos = xyz[pid, step]
            if not pos.any():
                break   # end of track
            # figure out process for this step
            proc = "?"
            nib_file = run / "seqnib.npy"
            if nib_file.is_file() and np.load(nib_file).ndim == 2:
                proc = LUT.get(np.load(nib_file)[pid, step], "?")
            else:
                word_idx, nib_idx = divmod(step, 16)
                packed_row = np.load(run / "seq.npy")[pid].ravel()
                if word_idx < packed_row.size:
                    nib = (int(packed_row[word_idx]) >> (4 * nib_idx)) & 0xF
                    proc = LUT.get(nib, "?")
            x, y, z = pos
            print(f"  step {step:2d}  {proc:12s}  (x={x:9.3f}, y={y:9.3f}, z={z:9.3f})")
        print()

    if ph_ids.size > to_show.size:
        more = ph_ids.size - to_show.size
        print(f"... {more} more hit-photons not shown (use --nmax to increase)")

if __name__ == "__main__":
    main()

