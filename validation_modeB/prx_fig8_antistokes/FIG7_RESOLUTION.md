# PRX Fig. 7: both remaining discrepancies are one sign in the legacy core-hole sector

## Summary

After fixing two real bugs in this code (the square-planar inter-ligand hopping sign, and all
three RIXS paths dropping the spin-flip channel -- see `../sqpl_cluster_quanty/` and
`../rixs_spinflip_channel/`), two differences against the published Fig. 7 remained:

1. the reference's incident-energy resonance **red-shifts by ~2 eV** after the pump
   (-4.51 -> -6.55 in code units; Omega_i ~ 855 -> ~853 eV in the published figure), and
2. its total intensity **drops to ~35%**.

Both come from the same place: the legacy Fortran gives the **ligand-ligand hopping the opposite
sign in the core-hole sector to the valence sector**.

## The measurement

XX polarization, channel 0, delay 0 (`-6.28`, pre-pump) against delay 14 (`+22.31`, post-pump),
over the Fig. 7 window `Omega_i` in [-8.75, -1.0], `omega` in [-2, 6]. "box" is the Fig. 7(d)
dashed box (`omega` in [-1.9, 0.5], `Omega_i` in [-8.702, -5.4]).

| | total ratio | incident peak | box ratio |
|---|---|---|---|
| reference (Fortran production) | 0.346 | -4.51 -> **-6.55** | 1.137 |
| ours, both fixes | 1.011 | -4.55 -> -4.75 | 1.516 |
| **ours + legacy core-hole sign** (`EX_HOP_SIGN=-1`) | 0.824 | -4.55 -> **-6.53** | 3.406 |
| ours, coherent final states (diagnostic) | 3.318 | -4.55 -> -4.75 | 3.594 |

The red shift is reproduced to **0.02 eV** by that sign alone, and by nothing else tried. In
particular the reference's *coherent* final-state summation -- the other known convention
difference, and the subject of `../../THEORY_CROSSCHECK_PRB99.md` -- produces **no shift at all**
and moves the intensity the wrong way (x3.3 instead of x0.35), so it explains neither.

The same sign also reproduces the reference's equilibrium XAS **exactly**, including its +1.9 eV
satellite: every window weight to four decimals, r = 1.0000
(`../sqpl_cluster_quanty/xas_vs_reference/`). One sign, both observables.

## Why it is a defect and not a convention

The reference's two sectors are, in its own notation,

    valence   :  xcef - 2*tmatx - 2*tmaty - 4*tmato      (line ~764)
    core-hole :  xcef + 2*tmatx + 2*tmaty + 4*tmato      (line ~1555)

byte-identical code apart from the three signs. Flipping the two TM-ligand terms alone is the
gauge transformation `G = diag(+1 x5 d, -1 x6 ligand)`, which leaves the ligand-ligand term
`tmato` invariant -- verified by explicit conjugation:

    G (xcef - 2TX - 2TY - 4TO) G  ==  xcef + 2TX + 2TY - 4TO      (identical spectrum)

The reference uses `+4TO`, so its core-hole Hamiltonian is **not** gauge-equivalent to its own
valence Hamiltonian. The one-body spectra differ by up to **0.70 eV**:

    valence / its own gauge :  -0.6381 -0.0502  0.9547  0.9796  0.9796  9.9453  9.9503 ...
    core-hole as written    :  -0.7511 -0.0419  0.9630  0.9796  0.9796  9.9000 10.6511 ...

An O-O hopping integral cannot depend on whether a core hole sits on the metal atom. The
Slater-Koster sign of `Vpppi` is a physical input to the model and must be the same in both
sectors; nothing in the geometry or in the core-hole potential changes it. (Our own inter-ligand
fix, by contrast, was settled *by* Slater-Koster: it concerned the relative sign of the two
surviving elements *within* `tmato`, which the Bloch sum fixes unambiguously.)

## What this means for the paper

The +1.9 eV XAS satellite and the ~2 eV pump-induced red shift of the L3 resonance in Fig. 7 both
appear to be artifacts of this sign rather than physics. The features the paper's text actually
interprets are unaffected by it: the elastic line, the dd excitation at ~1 eV (the eg-t2g site
energy difference), and the photoinduced anti-Stokes feature at ~-0.9 eV, which our code produces
either way (box ratio 1.5 with the physical sign, 3.4 with the legacy one).

This is recorded as an argument from the structure of the model, not as a measurement: no
experimental Ni L3 XAS is compared here. Anyone relying on the satellite or the red shift should
decide the point deliberately rather than inherit it.

## Reproducing

    # our physical convention
    cd map_sf_XX     && OMP_NUM_THREADS=4 main-test
    # the legacy core-hole convention
    cd map_exhop_XX  && OMP_NUM_THREADS=4 EX_HOP_SIGN=-1 main-test
    # figures
    python3 sf/plot_fig7.py    "$PWD/sf"    XX      # ours vs reference
    python3 exhop/plot_fig7.py "$PWD/exhop" XX      # ours + legacy sign vs reference

All runs use `FG = 6.0 4.635 6.18 2.63` (the corrected F0,G1,F2,G3 ordering) and `GSNEV = 400`,
which the spin-flip fix needs so that every ground block has a complete valence eigenbasis.
