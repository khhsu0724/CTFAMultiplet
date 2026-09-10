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

### Delay dependence, with the reference's own (irregular) delay grid

**The reference's delays are not a uniform grid.** Its `allowed_tm` list picks 19 irregular
times, `tmeas = pi*2*tm/timesteps` with `timesteps = 320`:

    tm  = -320, -52, 0, 32, 140, 160, 200, 284, 362, 446, 536, 626, 710, 788, 878, 968, 1052, 1136, 1214
    t   = -6.283, -1.021, 0, +0.628, +2.749, +3.142, +3.927, +5.576, +7.108, +8.757, +10.524,
          +12.291, +13.941, +15.472, +17.239, +19.007, +20.656, +22.305, +23.837

and the file index is `count`, so `trSPEC_XX_1710.dat` is `count = 17`, i.e. **t = +22.305**. Our
`TRDELAY = -6.2832 22.3053 15` is a *uniform* 15-point grid, so the two coincide only at
t = -6.283 (both index 0), t = +3.927 (ours 5 = ref 6) and t = +22.305 (ours 14 = **ref 17**).
Comparing index-to-index is wrong everywhere else; `plot_fig7.py` gets this right (it uses
`(0, 17)` for the reference against `(0, 14)` for ours), but an earlier version of the table below
did not, and reported the reference at t = +17.24 as though it were t = +22.31.

Total intensity in the Fig. 7 window, reference interpolated onto our times, each curve
normalized to its own t = -6.283 value, channel 0:

| t | reference | ours | + legacy sign | + coherent | + BOTH |
|---|---|---|---|---|---|
| -6.283 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| -0.157 | 0.967 | 1.283 | 1.300 | 2.155 | 1.496 |
| +1.885 | 0.700 | 1.253 | 1.210 | 2.420 | 1.641 |
| +3.927 | 0.454 | 1.195 | 1.084 | 2.455 | 1.798 |
| +8.011 | 0.362 | 1.173 | 0.963 | 3.173 | 2.395 |
| +14.137 | 0.369 | 1.168 | 0.953 | 3.793 | 2.775 |
| +22.305 | 0.385 | 1.011 | 0.824 | 3.318 | 2.298 |

| pre -> post (t = -6.283 -> +22.305) | ratio | incident peak |
|---|---|---|
| reference | **0.385** | -4.51 -> **-6.55** |
| ours, both fixes | 1.011 | -4.55 -> -4.75 |
| ours + legacy core-hole sign | 0.824 | -4.55 -> **-6.53** |
| ours + coherent final states | 3.318 | -4.55 -> -4.75 |
| ours + **BOTH** legacy conventions | 2.298 | -4.55 -> **-6.53** |

Two clean conclusions:

- **The red shift is the core-hole sign, and only that.** Every variant carrying it lands at
  -6.53 against the reference's -6.55; every variant without it stays at -4.75. The coherent
  final-state sum contributes nothing to the incident axis.
- **The intensity decay is explained by none of them, including the two together.** The reference
  falls to 0.385. The coherent sum makes our intensity *grow* (3.318); combined with the legacy
  sign it still grows (2.298). The reference's delay dependence is qualitatively opposite to what
  its own two known conventions produce in this code. Also ruled out: probe-window truncation (the
  reference is already on its plateau at the delays where our window is converged) and the
  first-order propagator (a ~2% occupation error at the production dt, see
  `../pump_dynamics_crosscheck/`).

**So this stays open.** Note our flat behaviour is what the paper's own text describes ("During the
pump, there is a slight increase in intensity... these changes are not pronounced"), and a 2.7x
collapse is not; but our pumped tr-RIXS has no external validation, so that is a hint, not a
verdict.

### The -0.9 eV anti-Stokes feature survives, and is cleaner with the physical sign

Post-pump (t = +22.305) energy-gain loss profile in the Fig. 7(d) box, each normalized to its own
maximum over -1.9 < omega < -0.1:

| | -1.31 | -1.12 | **-0.94** | -0.75 | -0.19 (elastic tail) |
|---|---|---|---|---|---|
| reference | 0.02 | 0.09 | **0.15** | 0.11 | 0.78 |
| **ours, both fixes** | 0.32 | 0.74 | **1.00** | 0.85 | 0.89 |
| ours + legacy sign | 0.08 | 0.19 | **0.30** | 0.33 | 0.90 |

In the reference it is a weak shoulder on the elastic tail; in our corrected calculation it is a
distinct peak at **-0.938 eV**, matching the paper's quoted -0.9 eV and its attribution to J_H and
the 1.0 eV crystal-field splitting. The legacy core-hole sign suppresses it. So the paper's central
Fig. 7/8 claim is intact and better supported by the corrected code than by the original.

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
