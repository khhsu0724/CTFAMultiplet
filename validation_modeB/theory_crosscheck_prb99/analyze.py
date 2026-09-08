#!/usr/bin/env python3
"""Reproduces every number quoted in THEORY_CROSSCHECK_PRB99.md section 3.

Run from this directory: python3 analyze.py
"""
import numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))

def grid(fn):
    d = np.loadtxt(os.path.join(HERE, fn), skiprows=1)
    ab, ls = np.unique(d[:, 0]), np.unique(d[:, 1])
    g = np.zeros((len(ls), len(ab)))
    for r in d:
        g[np.argmin(abs(ls - r[1])), np.argmin(abs(ab - r[0]))] = r[2]
    return ab, ls, g

def maxrel(f1, f2):
    a, b = grid(f1)[2], grid(f2)[2]
    return np.abs(a - b).max() / np.abs(a).max()

DELAYS = ["t-3.141593", "t0.000000", "t3.141593"]

print("3.1  Basis invariance: max relative change under TRRIXS2_FLIP_VAL=3")
print("     delay        corrected sum_l|A_l|^2      as-was |sum_l A_l|^2")
for t in DELAYS:
    f = "new_noflip/trRIXS2_chan0_%s.txt" % t, "new_flip3/trRIXS2_chan0_%s.txt" % t
    c = "coh_noflip/trRIXS2_chan0_%s.txt" % t, "coh_flip3/trRIXS2_chan0_%s.txt" % t
    print("     %-12s %8.2e                   %8.2e" % (t, maxrel(*f), maxrel(*c)))

print("\n3.2  Run-to-run reproducibility of the corrected result (two identical runs)")
for t in DELAYS:
    print("     %-12s det5a vs det5b: %8.2e" % (
        t, maxrel("det5a/trRIXS2_chan0_%s.txt" % t, "det5b/trRIXS2_chan0_%s.txt" % t)))

print("\n3.3  Pump-off Mode B vs equilibrium Kramers-Heisenberg (same single ground state)")
abE, lsE, EQf = grid("eq_fine/RIXS_Ledge_X_Y.txt")
abM, lsM, NEW = grid("poff_new/trRIXS2_chan0_t0.000000.txt")
_,   _,   COH = grid("poff_coh/trRIXS2_chan0_t0.000000.txt")
EQ = EQf[:, [np.argmin(abs(abE - a)) for a in abM]]
print("     incident profile (summed over loss, normalized):")
print("       ab    " + " ".join("%5.1f" % v for v in abM))
for nm, G in [("EQ    ", EQ), ("ModeB ", NEW)]:
    p = G.sum(axis=0); p /= p.max()
    print("       " + nm + " " + " ".join("%5.2f" % v for v in p))

def blur(G, s):
    x = lsM[:, None] - lsM[None, :]
    K = np.exp(-x**2 / (2 * s * s)); K /= K.sum(axis=0, keepdims=True)
    return K @ G

print("     2D correlation vs probe-broadened equilibrium:")
for s in (0.6, 0.8, 1.0, 1.2):
    B = blur(EQ, s); a = (B / B.sum()).ravel()
    print("       blur %.1f eV:  corrected r=%+.4f   as-was r=%+.4f" % (
        s, np.corrcoef(a, (NEW / NEW.sum()).ravel())[0, 1],
           np.corrcoef(a, (COH / COH.sum()).ravel())[0, 1]))
