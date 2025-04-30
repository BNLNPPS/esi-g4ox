#!/usr/bin/env python3
"""
dump_escaped_tracks.py
----------------------
Print every recorded step for photons that have **any point** outside an
R-mm sphere centred at the origin.

• Works with either seqnib.npy (grid or packed) or seq.npy (packed rows).  
• Does *not* use domain.npy – the escape test is purely spherical.

Process-nibble → letter map is Opticks’ default:
    1:B 2:R 3:A 4:S 5:T 6:C 7:D   (edit LUT as needed).
"""

from __future__ import annotations
import argparse, pathlib as pl, sys, numpy as np

# --------------------------------------------------------------------------- #
# Configurable nibble map (edit if your codes differ)                         #
# --------------------------------------------------------------------------- #
LUT = { 1:"Cherenkov", 2:"Scint", 3:"Miss", 4:"Bulk absorb", 5:"Bulk reemit", 6:"Bulk scatter", 7:"Surface detect"}

# --------------------------------------------------------------------------- #
# Sequence decoding (full, handles all variants)                              #
# --------------------------------------------------------------------------- #
def decode_sequences(run: pl.Path, n_ph: int) -> np.ndarray:
    """Return ndarray[str] of length n_ph with full interaction strings."""
    nib_file = run / "seqnib.npy"
    if nib_file.is_file():
        arr = np.load(nib_file)
        if arr.ndim == 2:                    # (N,16) nibble grid
            return np.array(["".join(LUT[n] for n in row if n) for row in arr])
        packed = arr.reshape(arr.shape[0], -1)   # weird 1-D ints
    else:
        packed = np.load(run / "seq.npy").reshape(n_ph, -1)  # (N,k)

    out=[]
    for row in packed:
        chars=[]
        for w in row:
            w = int(w)
            for _ in range(16):              # 16 nibbles/word
                nib = w & 0xF
                if nib == 0: w = 0; break
                chars.append(LUT.get(nib,"?")); w >>= 4
            if w==0: break
        out.append("".join(chars))
    return np.array(out)

# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", help="folder with record.npy + seq*.npy files")
    ap.add_argument("--radius", type=float, default=11.0,
                    help="sphere radius mm (default 11)")
    ap.add_argument("--nmax", type=int, default=0,
                    help="max photons to list (0=all)")
    args = ap.parse_args()

    run = pl.Path(args.run_dir)
    rec = np.load(run / "record.npy")              # (N, MAX, 4, 4)
    N, MAX = rec.shape[:2]
    xyz = rec[:, :, 0, :3]

    # 1) escape mask --------------------------------------------------------
    r2_lim = args.radius**2
    esc    = ((xyz**2).sum(-1) > r2_lim).any(1)
    ids    = np.flatnonzero(esc)
    if ids.size == 0:
        print(f"No photon crosses {args.radius} mm sphere."); return
    print(f"{ids.size} / {N} photons escape R={args.radius} mm\n")

    # 2) decode sequences ---------------------------------------------------
    seq_all = decode_sequences(run, N)

    # 3) step-by-step dump --------------------------------------------------
    show = ids if args.nmax==0 else ids[:args.nmax]
    for pid in show:
        seq = seq_all[pid]
        print(f"Photon {pid}  sequence = {seq or '<empty>'}")
        for step in range(MAX):
            pos = xyz[pid, step]
            if not pos.any(): break               # zero row → end of track
            # retrieve nibble for this step (works for packed or grid)
            proc = "?"
            if (run / "seqnib.npy").is_file() and np.load(run / "seqnib.npy").ndim==2:
                proc = LUT.get(np.load(run / "seqnib.npy")[pid, step], "?")
            else:
                # step→word,nibble index
                word_idx, nib_idx = divmod(step, 16)
                packed_row = np.load(run / "seq.npy")[pid].ravel()
                if word_idx < packed_row.size:
                    nib = (int(packed_row[word_idx]) >> (4*nib_idx)) & 0xF
                    proc = LUT.get(nib, "?")
            x, y, z = pos
            print(f"  step {step:2d}  {proc}   (x={x:9.3f}, y={y:9.3f}, z={z:9.3f})")
        print()

    if ids.size > show.size:
        print(f"... {ids.size - show.size} more photons not shown "
              f"(use --nmax or 0 for all)")

if __name__ == "__main__":
    main()

