#!/usr/bin/env python3
"""Sum the per-ground-state trRIXS2 channels into the physical cross section.

Mode B writes one set of files per DEGENERATE GROUND STATE (trRIXS2_chan<n>_t<delay>.txt).
The observable is their sum. Reading chan0 alone drops the rest -- for the PRX d8 square-planar
set that is 2 of 3 channels, and the Sz = 0 channel differs from the Sz = +-1 ones.

Usage:  sum_trrixs2_channels.py <run_dir> [out_dir]

Writes trRIXS2_total_t<delay>.txt with the same three columns (incident, loss, intensity),
plus trRIXS2_total_delays.txt listing the delays. Reports which channels it found and how much
each contributes, so a missing channel is obvious rather than silent.
"""
import sys, os, re, glob
import numpy as np

def main():
    if len(sys.argv) < 2:
        print(__doc__); return 1
    d = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else d
    os.makedirs(out, exist_ok=True)
    files = glob.glob(os.path.join(d, "trRIXS2_chan*_t*.txt"))
    if not files:
        print(f"no trRIXS2_chan*_t*.txt in {d}"); return 1
    pat = re.compile(r"trRIXS2_chan(\d+)_t(-?\d+\.\d+)\.txt$")
    byd = {}
    for f in files:
        m = pat.search(os.path.basename(f))
        if m: byd.setdefault(m.group(2), {})[int(m.group(1))] = f
    chans = sorted({c for v in byd.values() for c in v})
    print(f"found channels {chans} across {len(byd)} delay(s)")
    incomplete = [t for t, v in byd.items() if sorted(v) != chans]
    if incomplete:
        print(f"WARNING: {len(incomplete)} delay(s) are missing some channels and were SKIPPED:")
        for t in sorted(incomplete)[:10]:
            print(f"    t={t} has only {sorted(byd[t])}")
    delays = []
    for t in sorted(byd, key=float):
        if t in incomplete: continue
        arrs = [np.loadtxt(byd[t][c], skiprows=1) for c in chans]
        for a in arrs[1:]:
            if not (np.allclose(a[:, 0], arrs[0][:, 0]) and np.allclose(a[:, 1], arrs[0][:, 1])):
                print(f"  ERROR: t={t} channels are on different (incident, loss) grids"); return 1
        per = [a[:, 2].sum() for a in arrs]
        tot = arrs[0].copy()
        for a in arrs[1:]: tot[:, 2] += a[:, 2]
        fo = os.path.join(out, f"trRIXS2_total_t{float(t):f}.txt")
        with open(fo, "w") as fh:
            fh.write(f"{'absorption energy':>20s}{'energy loss':>20s}{'intensity':>20s}\n")
            for r in tot:
                fh.write(f"{r[0]:20.10g}{r[1]:20.10g}{r[2]:20.10g}\n")
        delays.append(float(t))
        frac = " ".join(f"{p/sum(per)*100:.1f}%" for p in per)
        print(f"  t={float(t):+10.5f}  total={tot[:,2].sum():.6e}   per-channel share: {frac}")
    with open(os.path.join(out, "trRIXS2_total_delays.txt"), "w") as fh:
        for t in delays: fh.write(f"{t:.10g}\n")
    print(f"wrote {len(delays)} summed delay file(s) to {out}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
