#!/usr/bin/env python3
"""PRX Fig.7(c)/(d) style maps: incident photon energy vs energy transfer, TM L3 edge,
before (dt << 0) and after (dt >> 0) the pump -- new code vs the PRX Fortran reference."""
import numpy as np, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
B    = os.path.expanduser("~/Documents/GradSchool/Devereaux Group/Herr_Jost/tr_spectroscopy/Ni-L-RIXS")
REFD = {"XX": B+"/new_shape_function/XX/16_points_new_shape_function_XX",
        "XY": B+"/new_shape_function/XY"}
D    = np.linspace(-6.2831853071795865, 22.305305514018902, 15)
BOX  = dict(w0=-1.9, w1=0.5, a0=-8.702, a1=-5.4)     # Fig.7(d) dashed box
POL  = sys.argv[2] if len(sys.argv) > 2 else "XX"

def togrid(d):
    ab, ls = np.unique(d[:, 0]), np.unique(d[:, 1])
    G = np.zeros((len(ls), len(ab)))
    for a, b_, c in zip(d[:, 0], d[:, 1], d[:, 2]):
        G[np.searchsorted(ls, b_), np.searchsorted(ab, a)] = c
    return ab, ls, G

ref  = lambda p, i: togrid(np.loadtxt(os.path.join(REFD[p], "trSPEC_%s_%04d.dat" % (p, 100*i+10))))
mine = lambda p, i: togrid(np.loadtxt(os.path.join(HERE, "map_%s" % p,
                                      "trRIXS2_chan0_t%f.txt" % D[i]), skiprows=1))

CM = plt.cm.afmhot_r
fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
rows = [("PRX Fortran reference", lambda i: ref(POL, i), (0, 17)),
        ("New code (spin-flip + tpp fixes)", lambda i: mine(POL, i), (0, 14))]

for r, (name, get, (i_pre, i_post)) in enumerate(rows):
    ab0, ls0, G0 = get(i_pre)
    ab1, ls1, G1 = get(i_post)
    sel = lambda ab: (ab >= -8.75) & (ab <= -1.0)
    wsel = lambda ls: (ls >= -2.0) & (ls <= 6.0)
    s0, w0 = sel(ab0), wsel(ls0)
    A0, A1 = G0[np.ix_(w0, s0)], G1[np.ix_(w0, s0)]
    vmax = max(A0.max(), A1.max())            # shared scale within a row
    for c, (A, lab, dt) in enumerate([(A0, r"$\Delta t \ll 0$", "-6.28"),
                                      (A1, r"$\Delta t \gg 0$", "+22.31")]):
        ax = axes[r, c]
        im = ax.pcolormesh(ls0[w0], ab0[s0], A.T, cmap=CM, vmin=0, vmax=vmax,
                           shading="nearest", rasterized=True)
        ax.add_patch(Rectangle((BOX["w0"], BOX["a0"]), BOX["w1"]-BOX["w0"], BOX["a1"]-BOX["a0"],
                               fill=False, ec="k", ls=(0, (3, 2)), lw=1.2))
        ax.axvline(-0.9375, color="#1f6fb4", ls=":", lw=1.1)
        ax.set_xlim(-2.0, 6.0); ax.set_ylim(-8.75, -1.0)
        ax.text(0.03, 0.955, "%s   %s  (%s eV$^{-1}$)" % (POL, lab, dt), transform=ax.transAxes,
                fontsize=9.5, va="top", fontweight="bold")
        if c == 0:
            ax.set_ylabel(r"$\Omega_i$  (eV, code units)" + "\n" + name, fontsize=9.5)
        ax.set_xlabel(r"energy transfer  $\omega$  (eV)", fontsize=9.5)
        ax.tick_params(labelsize=8.5)
        cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
        cb.ax.tick_params(labelsize=7.5)
        cb.formatter.set_powerlimits((0, 0)); cb.update_ticks()
        cb.ax.yaxis.get_offset_text().set_fontsize(7)
fig.suptitle("3d$^8$ TM $L_3$-edge trRIXS maps before and after the pump (PRX Fig. 7c/7d), "
             "%s polarization, channel 0" % POL, fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.955])
out = os.path.join(HERE, "fig7_maps_%s.png" % POL)
fig.savefig(out, dpi=165); print("wrote", out)
