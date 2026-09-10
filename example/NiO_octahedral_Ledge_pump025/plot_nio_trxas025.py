#!/usr/bin/env python3
"""NiO L3 trXAS at 0.25 eV: transient shift and recovery, two pulse envelopes."""
import glob, re, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUNS = [("$PUMPSIGMA = 1.2247$  (sub-cycle)", sys.argv[1] if len(sys.argv)>1 else "/tmp/nio_trxas025_sub", "#c0392b"),
        ("$PUMPSIGMA = 10.7$  (1 cycle)",     sys.argv[2] if len(sys.argv)>2 else "/tmp/nio_trxas025",     "#2471a3")]
fig, ax = plt.subplots(2, 1, figsize=(6.4, 5.4), sharex=True)
for lab, d, c in RUNS:
    fs = sorted(glob.glob(f"{d}/trXAS_Ledge_X_t*.txt"),
                key=lambda f: float(re.search(r"_t(-?\d+\.\d+)\.txt", f).group(1)))
    T, C, W = [], [], []
    for f in fs:
        t = float(re.search(r"_t(-?\d+\.\d+)\.txt", f).group(1))
        a = np.loadtxt(f, skiprows=1); m = (a[:,0] > -10) & (a[:,0] < 8)
        x, y = a[m,0], a[m,1]
        if y.sum() <= 0: continue
        T.append(t); C.append((x*y).sum()/y.sum()); W.append(y.sum())
    T, C, W = np.array(T), np.array(C), np.array(W)
    ax[0].plot(T, (C-C[0])*1000, "o-", ms=3, color=c, lw=1.3, label=lab)
    ax[1].plot(T, (W/W[0]-1)*100, "o-", ms=3, color=c, lw=1.3, label=lab)
ax[0].axhspan(-85, -65, color="grey", alpha=0.25, lw=0)
ax[0].text(52, -75, "expt ~100 meV", fontsize=7.5, color="dimgrey", va="center", ha="right")
ax[0].axhline(0, color="grey", lw=0.6)
ax[0].set_ylabel("$L_3$ centroid shift  (meV)")
ax[0].legend(fontsize=8, frameon=False, loc="lower right")
ax[0].set_title("NiO $L_3$ trXAS, 0.25 eV pump: transient shift, full recovery", fontsize=10)
ax[1].axhline(0, color="grey", lw=0.6)
ax[1].set_ylabel("$L_3$ intensity change  (%)"); ax[1].set_xlabel(r"$\Delta t$  (1/eV)")
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nio_trxas025.png")
fig.savefig(out, dpi=170); print("wrote", out)
