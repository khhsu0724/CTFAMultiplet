#!/usr/bin/env python3
"""PRX Fig. 3(d)-(f) analogue: 3d9 cluster TM L3-edge trXAS from this code."""
import glob, re, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DX = sys.argv[1] if len(sys.argv) > 1 else "/tmp/d9fine_X"
DY = sys.argv[2] if len(sys.argv) > 2 else "/tmp/d9fine_Y"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "d9_trxas_Ledge.png")

def load(d, pol):
    fs = sorted(glob.glob(f"{d}/trXAS_Ledge_{pol}_t*.txt"),
                key=lambda f: float(re.search(r"_t(-?\d+\.\d+)\.txt", f).group(1)))
    T, S = [], []
    for f in fs:
        t = float(re.search(r"_t(-?\d+\.\d+)\.txt", f).group(1))
        a = np.loadtxt(f, skiprows=1)
        T.append(t); S.append(a[:, 1]); E = a[:, 0]
    return np.array(T), E, np.array(S).T          # S[energy, delay]

T, E, SX = load(DX, "X"); _, _, SY = load(DY, "Y")
m = (E > -6) & (E < 8)
Em, SXm, SYm = E[m], SX[m], SY[m]
cen = lambda S: (Em[:, None]*S).sum(axis=0)/S.sum(axis=0)
wt  = lambda S: S.sum(axis=0)
OFF = 931.5 - cen(SXm)[0]                          # align pre-pump onto the published ~931.5 eV

fig = plt.figure(figsize=(6.2, 7.2))
gs  = fig.add_gridspec(3, 1, height_ratios=[2.1, 1, 1], hspace=0.42)

ax = fig.add_subplot(gs[0])
im = ax.pcolormesh(T, Em, SXm, cmap="afmhot_r", shading="nearest", rasterized=True)
ax.plot(T, cen(SXm), color="#1f6fb4", lw=1.3, label="L$_3$ centroid")
ax.set_ylabel(r"$\Omega_i$ (eV, code units)"); ax.set_xlim(T.min(), T.max())
ax.legend(fontsize=8, frameon=False, loc="upper right")
ax.set_title("3d$^9$ TM $L_3$-edge trXAS, X polarization  (cf. PRX Fig. 3d)", fontsize=10)
ax2 = ax.secondary_yaxis("right", functions=(lambda y: y+OFF, lambda y: y-OFF))
ax2.set_ylabel(r"$\Omega_i$ (eV, aligned)", fontsize=9)
fig.colorbar(im, ax=ax, pad=0.13, fraction=0.04).ax.tick_params(labelsize=7)

ax = fig.add_subplot(gs[1])
ax.plot(T, wt(SXm), color="k", lw=1.2)
ax.set_ylabel("integrated\n$L_3$ intensity", fontsize=9)
ax.set_xlim(T.min(), T.max()); ax.set_title("cf. PRX Fig. 3e", fontsize=9)

ax = fig.add_subplot(gs[2])
mm = T > 8
ax.plot(T[mm], wt(SXm)[mm]-wt(SXm)[mm].mean(), color="k", lw=1.2, label="X")
ax.plot(T[mm], wt(SYm)[mm]-wt(SYm)[mm].mean(), color="#c0392b", lw=1.2, label="Y")
ax.axhline(0, color="grey", lw=0.6)
ax.set_xlabel(r"$\Delta t$  (1/eV)"); ax.set_ylabel("intensity\n(mean removed)", fontsize=9)
ax.legend(fontsize=8, frameon=False, ncol=2)
r = np.corrcoef(wt(SXm)[mm], wt(SYm)[mm])[0, 1]
ax.set_title(f"cf. PRX Fig. 3f -- X/Y correlation {r:+.3f} (out of phase)", fontsize=9)

fig.tight_layout(); fig.savefig(OUT, dpi=170); print("wrote", OUT)
