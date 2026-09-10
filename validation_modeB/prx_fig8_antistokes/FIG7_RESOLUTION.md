# PRX Fig. 7: both remaining discrepancies are one sign in the legacy core-hole sector

## Summary

After fixing two real bugs in this code (the square-planar inter-ligand hopping sign, and all
three RIXS paths dropping the spin-flip channel -- see `../sqpl_cluster_quanty/` and
`../rixs_spinflip_channel/`), two differences against the published Fig. 7 remained:

1. the reference's incident-energy resonance **red-shifts by ~2 eV** after the pump
   (-4.51 -> -6.55 in code units; Omega_i ~ 855 -> ~853 eV in the published figure), and
2. its total intensity **drops to ~35%**.

**(1) is explained**: the legacy Fortran gives the **ligand-ligand hopping the opposite sign in
the core-hole sector to the valence sector**, and that reproduces the shift to 0.02 eV as well as
the equilibrium XAS satellite exactly. **(2) is still open** -- see the caveat below; it is not
truncation, not that sign, and not the coherent final-state sum.

## The measurement

XX polarization, delay 0 (`-6.28`, pre-pump) against delay 14 (`+22.31`, post-pump), over the
Fig. 7 window `Omega_i` in [-8.75, -1.0], `omega` in [-2, 6].

| | total ratio | incident peak |
|---|---|---|
| reference (Fortran production) | 0.346 | -4.51 -> **-6.55** |
| ours, both fixes | 1.011 | -4.55 -> -4.75 |
| **ours + legacy core-hole sign** (`EX_HOP_SIGN=-1`) | 0.824 | -4.55 -> **-6.53** |
| ours, coherent final states (diagnostic) | 3.318 | -4.55 -> -4.75 |

The red shift is reproduced to **0.02 eV** by that sign alone, and by nothing else tried. In
particular the reference's *coherent* final-state summation -- the other known convention
difference, and the subject of `../../THEORY_CROSSCHECK_PRB99.md` -- produces **no shift at all**
and moves the intensity the wrong way, so it explains neither.

The same sign also reproduces the reference's equilibrium XAS **exactly**, including its +1.9 eV
satellite: every window weight to four decimals, r = 1.0000
(`../sqpl_cluster_quanty/xas_vs_reference/`). One sign, two observables.

### Caveat on the intensity column: these delays are probe-window truncated

Our simulated window is `[delay_min, delay_max] +- 0.75*max(probe_sigma, tau_ch)` = [-10.10,
26.13], with `probe_sigma = 5.093`. **13 of the 15 delays sit less than 3 probe-sigma from an
edge**, including *both* Fig. 7 delays: delay 0 and delay 14 are each only **0.75 sigma** from one.
The reference's own window is [-25.13, +50.27] (`nb = 4`, `nt = 8`, `timesteps = 320`), so its
delays are fully converged. Only delays **6, 7 and 8** (+5.97, +8.01, +10.05) clear 3 sigma in our
runs, and measured against them our endpoint values are ~12% low.

So the *incident-peak position* above is robust -- it is a resonance energy, not a normalization --
but the total-intensity column is not, and a converged rerun is impractical: raising the window to
3x pushes `n_is` from 923 to 1506, and the K-row assembly scales as `n_is^2 * nval`, giving ~32 h
per channel (measured, not estimated: 7.8% of channel 0 in 2.5 h).

### What the converged delays do say

Total intensity in the Fig. 7 window, each curve normalized to its own delay-6 value (the first
converged delay), summed over the degenerate ground states as `2*chan0 + chan1` (`chan2` is
identical to `chan0` to r = 1.000000 by the Sz = +-1 symmetry):

| i | delay | conv. | reference | ours | ours + legacy sign | ours coherent |
|---|---|---|---|---|---|---|
| 0 | -6.28 | | 2.203 | 0.878 | 1.000 | 0.503 |
| 3 | -0.16 | | 1.933 | 1.099 | 1.314 | 0.942 |
| 4 | +1.88 | | 1.273 | 1.068 | 1.224 | 0.983 |
| 6 | +5.97 | **yes** | 1.000 | 1.000 | 1.000 | 1.000 |
| 7 | +8.01 | **yes** | 0.823 | 1.000 | 0.960 | 1.123 |
| 8 | +10.05 | **yes** | 0.800 | 0.986 | 0.944 | 1.157 |
| 14 | +22.31 | | 0.762 | 0.862 | 0.820 | 1.160 |

The reference decays monotonically by ~2.8x across the pump and saturates by delay 7. **Ours is
flat**, with a slight *rise* during the pump (1.07-1.10 at delays 3-4). That difference is not a
truncation artifact: in the converged region (6, 7, 8) the reference is already on its low plateau
while our curve is flat, and the legacy core-hole sign does not change this (0.944 vs the
reference's 0.800). The coherent final-state convention goes the wrong way entirely.

**This remains unexplained**, and it is the one quantity where we still disagree with the
reference. Note though that our behaviour is what the paper's own text describes: *"During the
pump, there is a slight increase in intensity... However, these changes are not pronounced,
suggesting that the immediate effect of the pump on the energy dynamics of the 3d8 cluster is
relatively modest."* A 2.8x collapse is not a modest change. Fig. 8's caption likewise reports a
flat pre-pump response and an intensity that "increases only marginally" at `t ~ t0`.

Also note the absolute normalizations differ by 12-38x (conventions in the (t1,t2) measure, the
probe normalization, the spin-flip channel, and the ground-state sum), so only shapes and ratios
are comparable.

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
