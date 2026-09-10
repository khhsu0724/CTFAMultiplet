#!/usr/bin/env python3
"""PRX Fig. 2-style charge oscillations for the NiO cluster: pump field + orbital
occupations, and their Fourier transforms.

Usage: plot_charge_oscillations.py <run_dir> [<run_dir2> ...] [out.png]

Panels follow PRX 15, 011012 Fig. 2 -- (a1) A(t), (a2) F[A(t)], (b1) occupations,
(b2) F[occupations] -- with two deliberate departures, because NiO at a 0.25 eV pump is a
much weaker perturbation than the PRX 3d9/3d8 cases at pi eV:

  * (b1) plots the CHANGE d<n> = <n>(t) - <n>(-inf), scaled by 1e3. The paper can use a raw
    0..1 axis because its occupations swing by ~0.7; here they move by ~0.02, which is
    invisible on that scale. Equilibrium values are given in the legend.
  * (b2) uses a log ordinate over 0-6 eV. The point of interest is which frequencies are
    ABSENT (the magnons) as much as which are present, and that needs decades.
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL = {"dx2": "k", "dz2": "#4b2e83", "pxx": "#c0159b", "pyy": "#f07f3c",
       "dxy": "#2ca02c", "dxz": "#1f77b4", "dyz": "#17becf", "pxy": "#8c564b"}
LAB = {"dx2": r"$d_{x^2-y^2}$", "dz2": r"$d_{z^2}$", "pxx": r"$p_{xx}$", "pyy": r"$p_{yy}$",
       "dxy": r"$d_{xy}$", "dxz": r"$d_{xz}$", "dyz": r"$d_{yz}$", "pxy": r"$p_{xy}$"}
# NiO reference energies: 1x and 2x MF (the magnon manifold), and the lowest dd / CT features
MARKS = [(0.1046, "magnon"), (0.2097, "2-magnon"), (1.0, r"$dd$"), (4.8, "CT")]

def spec(t, y, fmax):
    dt = t[1]-t[0]
    yy = (y - y.mean()) * np.hanning(len(y))
    f = np.fft.rfftfreq(len(t), d=dt) * 2*np.pi
    A = np.abs(np.fft.rfft(yy)) * dt
    m = (f <= fmax) & (f > 0)
    return f[m], A[m]

def panel(dirn, ax, title, nshow=5, fmax=6.0, tzoom=None):
    d = np.genfromtxt(os.path.join(dirn, "occupation_vs_time.txt"), names=True)
    t = d["t"]
    orbs = [n for n in d.dtype.names if n not in ("t", "E_inst", "E_stat")]
    amp = {o: d[o].max()-d[o].min() for o in orbs}
    show = [o for o in sorted(orbs, key=lambda o: -amp[o]) if amp[o] > 1e-5][:nshow]
    P = np.loadtxt(os.path.join(dirn, "pump_profile.txt"), skiprows=1)
    tp, Ap = P[:, 0], P[:, 1]
    if tzoom is None:
        w = tp[np.abs(Ap) > 0.02*np.abs(Ap).max()]
        tzoom = (w.min()-5, w.max()+60) if len(w) else (t.min(), t.max())

    ax[0][0].plot(tp, Ap, "k", lw=1.1); ax[0][0].set_xlim(*tzoom)
    ax[0][0].set_ylabel("$A(t)$"); ax[0][0].text(0.02, 0.88, "(a1)", transform=ax[0][0].transAxes, fontsize=9)
    ax[0][0].set_title(title, fontsize=10, loc="left")
    f, A = spec(tp, Ap, fmax)
    ax[0][1].semilogy(f, A, "k", lw=1.1, label="pump")
    ax[0][1].set_ylabel(r"$\mathcal{F}[A(t)]$"); ax[0][1].legend(fontsize=8, frameon=False)
    ax[0][1].text(0.02, 0.88, "(a2)", transform=ax[0][1].transAxes, fontsize=9)

    for o in show:
        n0 = d[o][0]
        ax[1][0].plot(t, (d[o]-n0)*1e3, color=COL.get(o), lw=1.0,
                      label=f"{LAB.get(o,o)}  ({n0:.3f})")
        f, A = spec(t[t > tzoom[0]], d[o][t > tzoom[0]], fmax)
        ax[1][1].semilogy(f, A, color=COL.get(o), lw=1.0, label=LAB.get(o, o))
    ax[1][0].set_xlim(*tzoom); ax[1][0].axhline(0, color="grey", lw=0.5)
    ax[1][0].set_ylabel(r"$\delta\langle n\rangle \times 10^{3}$"); ax[1][0].set_xlabel(r"$\Delta t$  (1/eV)")
    ax[1][0].text(0.02, 0.88, "(b1)", transform=ax[1][0].transAxes, fontsize=9)
    ax[1][0].legend(fontsize=7, frameon=False, ncol=2, title="orbital (equil. $\\langle n\\rangle$)",
                    title_fontsize=7)
    ax[1][1].set_ylabel(r"$\mathcal{F}[\langle n\rangle]$"); ax[1][1].set_xlabel("Frequency (eV)")
    ax[1][1].text(0.02, 0.88, "(b2)", transform=ax[1][1].transAxes, fontsize=9)
    ax[1][1].legend(fontsize=7, frameon=False, ncol=3)
    for a in (ax[0][1], ax[1][1]):
        for e, lab in MARKS:
            a.axvline(e, color="grey", ls=":", lw=0.8)
        a.set_xlim(0, fmax)
    for e, lab in MARKS:
        ax[1][1].annotate(lab, (e, 1), xycoords=("data", "axes fraction"),
                          xytext=(1.5, -2), textcoords="offset points",
                          fontsize=6.5, color="dimgrey", rotation=90, va="top")

if __name__ == "__main__":
    dirs = [a for a in sys.argv[1:] if os.path.isdir(a)]
    out = ([a for a in sys.argv[1:] if a.endswith(".png")] or
           [os.path.join(os.path.dirname(os.path.abspath(__file__)), "nio_charge_oscillations.png")])[0]
    fig, axes = plt.subplots(2*len(dirs), 2, figsize=(9.6, 3.0*2*len(dirs)))
    axes = np.atleast_2d(axes)
    for k, dn in enumerate(dirs):
        panel(dn, axes[2*k:2*k+2], os.path.basename(dn))
    fig.tight_layout(); fig.savefig(out, dpi=170); print("wrote", out)
