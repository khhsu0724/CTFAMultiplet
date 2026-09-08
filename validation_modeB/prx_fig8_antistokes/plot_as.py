#!/usr/bin/env python3
"""Anti-Stokes cut of the 3d8 TM-L3 trRIXS: new code (corrected Mode B) vs the PRX Fortran reference.

For each code and polarization the cut is taken at the incident energy where the anti-Stokes
signal is strongest AT LONG DELAY (max intensity for loss in [-1.0,-0.75] eV). Those incident
energies are NOT the same in the two codes (-5.76 vs -7.39/-7.27), which is itself a result.

Run from this directory:  python3 plot_as.py
"""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
B    = os.path.expanduser("~/Documents/GradSchool/Devereaux Group/Herr_Jost/tr_spectroscopy/Ni-L-RIXS")
REFD = {"XX": B+"/new_shape_function/XX/16_points_new_shape_function_XX",
        "XY": B+"/new_shape_function/XY"}
MYD  = {"XX": "st_XX_dt2", "XY": "st_XY_dt2"}
D    = np.linspace(-6.2831853071795865, 22.305305514018902, 15)   # our delay grid
ASW  = (-1.0, -0.75)

def togrid(d):
    ab, ls = np.unique(d[:, 0]), np.unique(d[:, 1])
    G = np.zeros((len(ls), len(ab)))
    for a, b_, c in zip(d[:, 0], d[:, 1], d[:, 2]):
        G[np.searchsorted(ls, b_), np.searchsorted(ab, a)] = c
    return ab, ls, G

ref  = lambda p, i: togrid(np.loadtxt(os.path.join(REFD[p], "trSPEC_%s_%04d.dat" % (p, 100*i+10))))
mine = lambda p, i: togrid(np.loadtxt(os.path.join(HERE, MYD[p], "trRIXS2_chan0_t%f.txt" % D[i]), skiprows=1))

fig, axes = plt.subplots(2, 3, figsize=(15, 7.2), sharex=True)
C = {"pre": "#8a8a8a", "post": "#c0392b", "ours": "#1f6fb4", "ref": "#c0392b"}

for r, pol in enumerate(("XX", "XY")):
    cuts = {}
    for c, (name, key, pre_f, post_f) in enumerate([
            ("New code (corrected Mode B)", "ours", lambda p=pol: mine(p, 0), lambda p=pol: mine(p, 14)),
            ("PRX Fortran reference",       "ref",  lambda p=pol: ref(p, 0),  lambda p=pol: ref(p, 17))]):
        ax = axes[r, c]
        ab, ls, G0 = pre_f(); _, _, G1 = post_f()
        m  = (ls >= ASW[0]) & (ls <= ASW[1])
        ix = int(np.argmax(G1[m, :].max(axis=0)))          # incident energy of strongest AS
        w  = (ls >= -2.0) & (ls <= 0.30)
        ax.plot(ls[w], G0[w, ix], color=C["pre"],  lw=1.6, marker="o", ms=3,
                label=r"before pump  $\Delta t=-6.28$ eV$^{-1}$")
        ax.plot(ls[w], G1[w, ix], color=C["post"], lw=1.8, marker="o", ms=3,
                label=r"after pump  $\Delta t=+22.3$ eV$^{-1}$")
        ax.axvline(-0.9375, color="k", ls=":", lw=1.1)
        pk = G1[m, ix].max()
        ax.set_ylim(0, 1.45*pk); ax.set_xlim(-2.0, 0.30)   # elastic line deliberately off-scale
        ax.text(0.03, 0.94, "%s   %s" % (pol, name), transform=ax.transAxes,
                fontsize=9.5, va="top", fontweight="bold")
        ax.text(0.03, 0.855, r"$\Omega_i=%+.2f$ eV" % ab[ix], transform=ax.transAxes,
                fontsize=8.5, va="top", color="#444")
        ax.text(-0.9375, 1.38*pk, r"$-0.94$ eV", ha="center", fontsize=8)
        ax.tick_params(labelsize=8.5)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(8)
        if c == 0: ax.set_ylabel("RIXS intensity (arb. units)", fontsize=10)
        if r == 0 and c == 0:
            ax.legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(0.02, 0.84), frameon=False)
        cuts[key] = (ls[w], G1[w, ix], ab[ix])
    ax = axes[r, 2]
    for key, lbl, st in [("ours", "New code", "-"), ("ref", "Fortran reference", "--")]:
        x, y, a = cuts[key]
        mm = (x >= -1.6) & (x <= -0.45)
        ax.plot(x, y/y[mm].max(), color=C[key], ls=st, lw=1.9, marker="o", ms=3,
                label=r"%s ($\Omega_i=%+.2f$)" % (lbl, a))
    ax.axvline(-0.9375, color="k", ls=":", lw=1.1)
    ax.set_ylim(0, 1.6); ax.set_xlim(-2.0, 0.30)
    ax.text(0.03, 0.94, "%s   overlay, each normalized to its own AS peak" % pol,
            transform=ax.transAxes, fontsize=9.5, va="top", fontweight="bold")
    ax.legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(0.02, 0.85), frameon=False)
    ax.tick_params(labelsize=8.5)

for c in range(3):
    axes[1, c].set_xlabel(r"energy transfer  $\omega$  (eV)", fontsize=10)
fig.suptitle("3d$^8$ TM $L_3$-edge trRIXS, anti-Stokes side — cut at each code's own incident energy "
             "of strongest anti-Stokes signal (channel 0, $2S_z=-2$)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.955])
out = os.path.join(HERE, "antistokes_cut.png")
fig.savefig(out, dpi=165)
print("wrote", out)
