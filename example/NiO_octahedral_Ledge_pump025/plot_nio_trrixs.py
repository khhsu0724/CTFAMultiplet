#!/usr/bin/env python3
"""NiO trRIXS maps (Mode B) at selected delays, with on-resonance loss linecuts.

Usage: plot_nio_trrixs.py out.png "<label>::<file>" ...   (NROW env sets the map rows)

Rows 1..NROW are the maps. The final row adds one linecut panel per map row: the loss
profile integrated over incident energies within +-1.5 eV of that map's own resonance,
on a log ordinate, with the anti-Stokes side shaded.

The dotted grey curve is the MIRROR of the pre-pump Stokes side, I(-w), plotted on the
anti-Stokes side. At T = 0 there is no real energy-gain scattering, so whatever appears at
negative loss is the instrumental tail of the elastic line -- which is symmetric. Any
genuine anti-Stokes weight has to rise ABOVE that mirror.
"""
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAG = (0.1046, 0.2097)

def grid(path):
    d = np.loadtxt(path, skiprows=1)
    ab, ls = np.unique(d[:, 0]), np.unique(d[:, 1])
    G = np.zeros((len(ls), len(ab)))
    for a, b, c in d:
        G[np.searchsorted(ls, b), np.searchsorted(ab, a)] = c
    return ab, ls, G

def oncut(ab, ls, G, half=1.5):
    res = ab[np.argmax(G.sum(axis=0))]
    m = np.abs(ab - res) <= half
    return res, G[:, m].sum(axis=1)

def main(out, entries, nrow):
    ncol = len(entries)//nrow
    fig, axes = plt.subplots(nrow+1, ncol, figsize=(4.3*ncol, 3.3*nrow + 3.4), squeeze=False)
    cuts = [[] for _ in range(nrow)]
    for k, (lab, path) in enumerate(entries):
        r, c = divmod(k, ncol)
        ab, ls, G = grid(path)
        ax = axes[r][c]
        im = ax.pcolormesh(ls, ab, G.T, cmap="afmhot_r", shading="nearest", rasterized=True)
        ax.axvspan(ls.min(), 0, color="#2471a3", alpha=0.10, lw=0)
        for e in MAG:
            ax.axvline(e, color="#2471a3", ls=":", lw=0.8); ax.axvline(-e, color="#c0392b", ls=":", lw=0.8)
        res, cut = oncut(ab, ls, G)
        ax.axhline(res, color="#2e8b57", lw=0.9, ls="-.")
        ax.set_title(f"{lab}   ($\\Omega_i^{{res}}$ = {res:.2f})", fontsize=9)
        ax.set_xlabel(r"energy transfer $\omega$ (eV)", fontsize=8.5)
        if c == 0: ax.set_ylabel(r"$\Omega_i$ (eV, code units)", fontsize=8.5)
        ax.tick_params(labelsize=7.5)
        cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046); cb.ax.tick_params(labelsize=6.5)
        cb.formatter.set_powerlimits((0, 0)); cb.update_ticks()
        cb.ax.yaxis.get_offset_text().set_fontsize(6)
        cuts[r].append((lab, ls, cut))

    COLS = ["#111111", "#c0392b", "#2471a3", "#2e8b57"]
    for r in range(nrow):
        ax = axes[nrow][r]
        for i, (lab, ls, cut) in enumerate(cuts[r]):
            ax.semilogy(ls, cut, color=COLS[i % len(COLS)], lw=1.3, label=lab)
        ax.axvspan(-2.5, 0, color="#2471a3", alpha=0.10, lw=0)
        for e in MAG: ax.axvline(-e, color="#c0392b", ls=":", lw=0.9)
        ax.set_xlim(-2.0, 2.0)
        top = max(c.max() for _, _, c in cuts[r])
        ax.set_ylim(top*3e-6, top*2)
        ax.set_xlabel(r"energy transfer $\omega$ (eV)"); ax.set_ylabel("on-resonance intensity")
        ax.set_title("loss cut, $\Omega_i^{res}\pm1.5$ eV (anti-Stokes shaded)", fontsize=8.5)
        ax.legend(fontsize=7, frameon=False, loc="lower center")

    # ratio to each amplitude's own pre-pump map: at T = 0 the pre-pump map has no real
    # energy-gain scattering, so it IS the instrumental baseline and this divides it out.
    ax = axes[nrow][nrow] if nrow < ncol else axes[nrow][ncol-1]
    for r in range(nrow):
        lab0, ls0, ref = cuts[r][0]
        amp = lab0.split(",")[0]
        for i, (lab, ls, cut) in enumerate(list(enumerate(cuts[r]))[1:] and cuts[r][1:]):
            st = "-" if r == 0 else "--"
            ax.semilogy(ls, cut/np.maximum(ref, ref.max()*1e-12), ls=st,
                        color=COLS[(i+1) % len(COLS)], lw=1.3,
                        label=f"{amp}, {lab.split(', ')[1]}")
    ax.axhline(1.0, color="grey", lw=0.8)
    ax.axvspan(-2.5, 0, color="#2471a3", alpha=0.10, lw=0)
    for e in MAG:
        ax.axvline(-e, color="#c0392b", ls=":", lw=0.9)
        ax.annotate(f"$-${e:.3f}", (-e, 1), xycoords=("data", "axes fraction"),
                    xytext=(2, -2), textcoords="offset points", fontsize=6.5,
                    color="#c0392b", rotation=90, va="top")
    ax.set_xlim(-2.0, 1.0); ax.set_ylim(0.3, 100)
    ax.set_xlabel(r"energy transfer $\omega$ (eV)")
    ax.set_ylabel(r"$I(\Delta t)\,/\,I_{\rm pre}$")
    ax.set_title("pump-induced enhancement (solid $A_0$=5, dashed $A_0$=25.3)", fontsize=8.5)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right", ncol=1)

    for c in range(ncol):
        if c > nrow: axes[nrow][c].axis("off")
    fig.suptitle("NiO $L_3$ trRIXS (Mode B), 0.25 eV pump", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(out, dpi=170); print("wrote", out)

if __name__ == "__main__":
    out = sys.argv[1]
    entries = [(a.split("::", 1)[0], a.split("::", 1)[1]) for a in sys.argv[2:] if "::" in a]
    main(out, entries, int(os.environ.get("NROW", "2")))
