#!/usr/bin/env python3
"""PRX Fig. 2-style charge oscillations for the NiO cluster.

Usage: plot_charge_oscillations.py <run_dir> [out.png]

Panels follow PRX 15, 011012 Fig. 2:
    (a1) A(t)                      (a2) F[A(t)]
    (b1) <n>, e_g orbitals         (b2) F[<n>]
    (c1) <n>, t2g + ligands        (c2) total hole count (conservation check)

Occupations are HOLE occupations per orbital, summed over spin, so each runs 0..2 and they
sum to the input `Holes`. NiO has Holes = 2, the e_g orbitals carry ~0.95 each and everything
else is small -- hence the split y-ranges in (b1)/(c1) rather than one 0..1 axis: the paper can
use 0..1 because its occupations swing by ~0.7, whereas NiO at 0.25 eV moves by ~0.02.
Absolute hole units throughout; no scaling.
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL = {"dx2": "k", "dz2": "#4b2e83", "pxx": "#c0159b", "pyy": "#f07f3c",
       "pzz": "#7f7f7f", "dxy": "#2ca02c", "dxz": "#1f77b4", "dyz": "#17becf"}
LAB = {"dx2": r"$d_{x^2-y^2}$", "dz2": r"$d_{z^2}$", "pxx": r"$p_{xx}$", "pyy": r"$p_{yy}$",
       "pzz": r"$p_{zz}$", "dxy": r"$d_{xy}$", "dxz": r"$d_{xz}$", "dyz": r"$d_{yz}$"}
EG, REST = ["dx2", "dz2"], ["pxx", "pyy", "pzz", "dxy", "dxz", "dyz"]
MARKS = [(0.1046, "magnon"), (0.2097, "2-magnon"), (1.0, r"$dd$"), (4.8, "CT")]

def spec(t, y, fmax):
    dt = t[1]-t[0]
    f = np.fft.rfftfreq(len(t), d=dt) * 2*np.pi
    A = np.abs(np.fft.rfft((y - y.mean())*np.hanning(len(y)))) * dt
    m = (f <= fmax) & (f > 0)
    return f[m], A[m]

def main(dirn, out, fmax=6.0):
    d = np.genfromtxt(os.path.join(dirn, "occupation_vs_time.txt"), names=True)
    t = d["t"]
    orbs = [n for n in d.dtype.names if n not in ("t", "E_inst", "E_stat")]
    P = np.loadtxt(os.path.join(dirn, "pump_profile.txt"), skiprows=1)
    tp, Ap = P[:, 0], P[:, 1]
    w = tp[np.abs(Ap) > 0.02*np.abs(Ap).max()]
    tz = (w.min()-5, w.max()+60) if len(w) else (t.min(), t.max())

    fig, ax = plt.subplots(3, 2, figsize=(9.6, 8.2))
    ax[0][0].plot(tp, Ap, "k", lw=1.1); ax[0][0].set_xlim(*tz)
    ax[0][0].set_ylabel("$A(t)$"); ax[0][0].text(.02, .86, "(a1)", transform=ax[0][0].transAxes, fontsize=9)
    ax[0][0].set_title(f"NiO, 0.25 eV pump  ({os.path.basename(dirn)})", fontsize=10, loc="left")
    f, A = spec(tp, Ap, fmax)
    ax[0][1].semilogy(f, A, "k", lw=1.1, label="pump"); ax[0][1].legend(fontsize=8, frameon=False)
    ax[0][1].set_ylabel(r"$\mathcal{F}[A(t)]$")
    ax[0][1].text(.02, .86, "(a2)", transform=ax[0][1].transAxes, fontsize=9)

    for grp, a, tag in ((EG, ax[1][0], "(b1)"), (REST, ax[2][0], "(c1)")):
        for o in [x for x in grp if x in orbs and d[x].max() > 1e-4]:
            a.plot(t, d[o], color=COL.get(o), lw=1.0, label=LAB.get(o, o))
        a.set_xlim(*tz); a.set_ylabel(r"hole occupation  $\langle n \rangle$")
        a.legend(fontsize=7.5, frameon=False, ncol=3)
        a.text(.02, .86, tag, transform=a.transAxes, fontsize=9)
    ax[1][0].set_title("$e_g$ orbitals", fontsize=9, loc="right")
    ax[2][0].set_title("$t_{2g}$ + ligands", fontsize=9, loc="right")
    ax[2][0].set_xlabel(r"$\Delta t$  (1/eV)")

    for o in [x for x in EG+REST if x in orbs and d[x].max() > 1e-4]:
        f, A = spec(t[t > tz[0]], d[o][t > tz[0]], fmax)
        ax[1][1].semilogy(f, A, color=COL.get(o), lw=1.0, label=LAB.get(o, o))
    ax[1][1].set_ylabel(r"$\mathcal{F}[\langle n\rangle]$"); ax[1][1].legend(fontsize=7, frameon=False, ncol=4)
    ax[1][1].text(.02, .86, "(b2)", transform=ax[1][1].transAxes, fontsize=9)
    for e, lab in MARKS:
        for a in (ax[0][1], ax[1][1]): a.axvline(e, color="grey", ls=":", lw=0.8)
        ax[1][1].annotate(lab, (e, 1), xycoords=("data", "axes fraction"), xytext=(1.5, -2),
                          textcoords="offset points", fontsize=6.5, color="dimgrey",
                          rotation=90, va="top")
    for a in (ax[0][1], ax[1][1]): a.set_xlim(0, fmax)
    ax[1][1].set_xlabel("Frequency (eV)")

    tot = sum(d[o] for o in orbs)
    ax[2][1].plot(t, tot, "k", lw=1.0)
    ax[2][1].set_ylim(tot.mean()-5e-6, tot.mean()+5e-6)
    ax[2][1].set_ylabel("total holes"); ax[2][1].set_xlabel(r"$\Delta t$  (1/eV)")
    ax[2][1].text(.02, .86, "(c2)", transform=ax[2][1].transAxes, fontsize=9)
    ax[2][1].set_title(f"conserved to {np.abs(tot-tot[0]).max():.1e}", fontsize=9, loc="right")
    fig.tight_layout(); fig.savefig(out, dpi=170); print("wrote", out)

if __name__ == "__main__":
    dirs = [a for a in sys.argv[1:] if os.path.isdir(a)]
    out = ([a for a in sys.argv[1:] if a.endswith(".png")] or
           [os.path.join(os.path.dirname(os.path.abspath(__file__)), "nio_charge_oscillations.png")])[0]
    main(dirs[0], out)
