#!/usr/bin/env python3
"""PRX Fig. 7(c)/(d)-style incident-energy maps: 3d8 TM L3 edge, before and after the pump.

Fixes three things the earlier plotting got wrong or left out:

  * DELAY MAPPING. The reference's delays are an irregular 19-point list (`allowed_tm` in
    TD_Cluster_Xray.f90) and its file index is `count`, so trSPEC_*_1710.dat is t = +22.305 --
    NOT trSPEC_*_1410.dat (t = +17.239). Our uniform TRDELAY grid coincides with the reference
    at t = -6.283 (both index 0), +3.927 and +22.305 (ours 14 = ref 17) only.
  * CHANNEL SUM. Mode B writes one file set per degenerate ground state; the cross section is
    their sum. chan2 is identical to chan0 (Sz = +-1 symmetry, verified to 1e-10), so the sum is
    2*chan0 + chan1 when chan2 was not run.
  * INCIDENT-ENERGY AXIS. Each row gets its own eV axis, offset so that ITS pre-pump resonance
    sits at 855.0 eV (the published Fig. 7(c) position). The post-pump shift is then a real
    shift, comparable row to row. This is an alignment, not an absolute-energy calculation.

Usage: plot_fig7_maps.py [XX|XY]
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

POL  = sys.argv[1] if len(sys.argv) > 1 else "XX"
HERE = os.path.dirname(os.path.abspath(__file__))
B    = os.path.expanduser("~/Documents/GradSchool/Devereaux Group/Herr_Jost/tr_spectroscopy/Ni-L-RIXS")
REFD = {"XX": B + "/new_shape_function/XX/16_points_new_shape_function_XX",
        "XY": B + "/new_shape_function/XY"}[POL]

ALLOWED = [-320,-52,0,32,140,160,200,284,362,446,536,626,710,788,878,968,1052,1136,1214]
RT      = np.array([np.pi*tm/160 for tm in ALLOWED])      # reference delay times
D       = np.linspace(-6.2831853071795865, 22.305305514018902, 15)   # our TRDELAY grid
REF_PRE, REF_POST = 0, 17          # reference counts at t = -6.283 and +22.305
OUR_PRE, OUR_POST = 0, 14          # our indices at the same two times
BOX = dict(w0=-1.9, w1=0.5, a0=-8.702, a1=-5.4)
ANTISTOKES = -0.9375

def togrid(d):
    ab, ls = np.unique(d[:, 0]), np.unique(d[:, 1])
    G = np.zeros((len(ls), len(ab)))
    for a, b, c in zip(d[:, 0], d[:, 1], d[:, 2]):
        G[np.searchsorted(ls, b), np.searchsorted(ab, a)] = c
    return ab, ls, G

def ref(count):
    return togrid(np.loadtxt(os.path.join(REFD, "trSPEC_%s_%04d.dat" % (POL, 100*count+10))))

def ours(dirn, idx):
    """Channel-summed. chan2 == chan0 exactly, so fall back to 2*chan0 + chan1."""
    got = {}
    for c in (0, 1, 2):
        f = os.path.join(HERE, dirn, "trRIXS2_chan%d_t%f.txt" % (c, D[idx]))
        if os.path.exists(f): got[c] = np.loadtxt(f, skiprows=1)
    if not got: raise FileNotFoundError(dirn)
    ab, ls, G = togrid(got[0]); G = G.copy()
    if 1 in got: G += togrid(got[1])[2]
    G += togrid(got[2])[2] if 2 in got else togrid(got[0])[2]
    return ab, ls, G

ROWS = [("PRX Fortran reference", lambda w: ref(REF_PRE if w == "pre" else REF_POST)),
        ("New code (corrected)",  lambda w: ours("map_sf_%s" % POL, OUR_PRE if w == "pre" else OUR_POST))]
if os.path.isdir(os.path.join(HERE, "map_exhop_%s" % POL)):
    ROWS.append(("New code + legacy\ncore-hole sign",
                 lambda w: ours("map_exhop_%s" % POL, OUR_PRE if w == "pre" else OUR_POST)))

fig, axes = plt.subplots(len(ROWS), 2, figsize=(9.6, 3.35*len(ROWS)), squeeze=False)
CM = plt.cm.afmhot_r
for r, (name, get) in enumerate(ROWS):
    ab0, ls0, G0 = get("pre"); ab1, ls1, G1 = get("post")
    s = (ab0 >= -8.75) & (ab0 <= -1.0); w = (ls0 >= -2.0) & (ls0 <= 6.0)
    A0, A1 = G0[np.ix_(w, s)], G1[np.ix_(w, s)]
    vmax = max(A0.max(), A1.max())                        # shared within the row
    res0 = ab0[s][np.argmax(A0.sum(axis=0))]              # this row's pre-pump resonance
    res1 = ab1[s][np.argmax(A1.sum(axis=0))]
    off  = 855.0 - res0                                   # align it to the published 855 eV
    for c, (A, lab, t) in enumerate([(A0, r"$\Delta t \ll 0$", -6.283),
                                     (A1, r"$\Delta t \gg 0$", +22.305)]):
        ax = axes[r][c]
        im = ax.pcolormesh(ls0[w], ab0[s], A.T, cmap=CM, vmin=0, vmax=vmax,
                           shading="nearest", rasterized=True)
        ax.add_patch(Rectangle((BOX["w0"], BOX["a0"]), BOX["w1"]-BOX["w0"], BOX["a1"]-BOX["a0"],
                               fill=False, ec="k", ls=(0, (3, 2)), lw=1.1))
        ax.axvline(ANTISTOKES, color="#1f6fb4", ls=":", lw=1.1)
        res = res0 if c == 0 else res1
        ax.axhline(res, color="#2e8b57", lw=0.9, ls="-.")
        ax.text(5.85, res, " %+.2f " % res, fontsize=7, color="#2e8b57", ha="right", va="bottom")
        ax.set_xlim(-2.0, 6.0); ax.set_ylim(-8.75, -1.0)
        ax.text(0.025, 0.96, "%s   %s  ($t=%+.2f$)" % (POL, lab, t), transform=ax.transAxes,
                fontsize=9, va="top", fontweight="bold")
        ax.set_xlabel(r"energy transfer  $\omega$  (eV)", fontsize=8.5)
        ax.tick_params(labelsize=7.5)
        if c == 0:
            ax.set_ylabel(r"$\Omega_i$ (eV, code units)" + "\n" + name, fontsize=8.5)
        ax2 = ax.secondary_yaxis("right", functions=(lambda y, o=off: y+o, lambda y, o=off: y-o))
        ax2.set_yticks([852, 854, 856]); ax2.tick_params(labelsize=7.5)
        if c == 1: ax2.set_ylabel(r"$\Omega_i$ (eV, aligned)", fontsize=8)
        cb = fig.colorbar(im, ax=ax, pad=0.16, fraction=0.045)
        cb.ax.tick_params(labelsize=6.5)
        cb.formatter.set_powerlimits((0, 0)); cb.update_ticks()
        cb.ax.yaxis.get_offset_text().set_fontsize(6)
    axes[r][1].text(0.025, 0.055, "row shares one colour scale;  shift %+.2f eV" % (res1-res0),
                    transform=axes[r][1].transAxes, fontsize=7.5, color="#2e8b57")
fig.suptitle("3d$^8$ TM $L_3$-edge trRIXS, %s polarization, channel-summed  "
             "(PRX Fig. 7c/7d)\ndashed box = Fig. 7(d) integration window;  "
             "dotted line = $-0.94$ eV anti-Stokes;  dash-dot = resonance" % POL, fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 1 - 0.055*3/len(ROWS)])
out = os.path.join(HERE, "fig7_maps_%s_corrected.png" % POL)
fig.savefig(out, dpi=170); print("wrote", out)
