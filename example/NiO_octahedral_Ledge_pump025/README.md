# NiO octahedral cluster, 0.25 eV pump, TM L3 edge

Same Hamiltonian as `../NiO_octahedral_Ledge/` (`MF = 0.11`, valence SOC on, `tpd = 0.35`,
`tpp = 0`, `MLCT = 4.5`), same pump amplitude and envelope shape as the PRX runs, with
`PUMPFREQ = 0.25` for comparison with a 0.25 eV pump experiment.

## Why 0.25 eV is an interesting choice here

The NiO ground sector has **low-energy excitations at 0.1048 and 0.2100 eV**, which are
1x and 2x `MF = 0.11` (ratio 0.955 after hybridization/SOC renormalization) -- the Sz ladder of
the S=1 Ni(2+) in the molecular field, i.e. the magnon manifold. The lowest charge/dd excitation
is at **0.934 eV**. So a 0.25 eV drive sits just above the two-magnon and far below any charge
excitation.

## Envelope: 0.25 eV needs a wider pulse than the PRX runs

`PUMPSIGMA` is a width in eV^-1, so the number of field cycles inside the pulse depends on
`PUMPFREQ`:

| `PUMPFREQ` | period | cycles in FWHM (`PUMPSIGMA = 1.2247`) | omega / bandwidth |
|---|---|---|---|
| pi (PRX) | 2.00 | 1.44 | 3.85 |
| 4.75 | 1.32 | 2.18 | 5.82 |
| **0.25** | **25.13** | **0.115** | **0.31** |

At 0.25 eV with the PRX envelope the field completes about a **tenth of a cycle** and the pulse
bandwidth (~0.82 eV) is 3x the carrier: `sin(w t) ~ w t` across the pulse, so it is a single
asymmetric sub-cycle kick rather than a 0.25 eV drive. For roughly one cycle in the FWHM use
`PUMPSIGMA = 10.7`; two cycles needs 21.3, three needs 32.0 (and the pulse then spans ~190 eV^-1,
which makes Mode B trRIXS very expensive since its cost goes as n_is^2).

Both are provided; neither is wrong, they are different experiments.

## Result: a transient red shift that fully recovers

L3 centroid over -10 < Omega_i < 8, `PUMPAMPX = 5`, `TAUCH = 2.5` (matching the equilibrium
`epsab = 0.4`):

| envelope | max transient shift | residual shift | peak intensity change |
|---|---|---|---|
| `PUMPSIGMA = 1.2247` (sub-cycle) | **-42.0 meV** at `dt = 0` | -0.2 meV | +1.0 % |
| `PUMPSIGMA = 10.7` (1 cycle) | **-96.5 meV** at `dt = -16` | +0.2 meV | +2.5 % |

Both **recover completely** -- residual under 1 meV. That is the signature of driving below the
charge gap: the pump modulates the Hamiltonian strongly while it is on (the instantaneous
ground-state energy swings by up to 0.53 eV for the 1-cycle pulse) but leaves almost nothing
behind (permanently deposited energy 5e-4 eV, versus 1.19 eV for the PRX 3.14 eV pump on the 3d9
cluster).

This is qualitatively unlike the PRX pump: there the deposited energy is large and the shift is
**permanent**; here it is adiabatic and **transient**. Both magnitudes land in the range of the
~100 meV measured on a square-planar 3d9 nickelate, so this is the regime to compare experiment in.

## Caution: the amplitude is nonperturbative and not monotonic

`PUMPAMPX` is a dimensionless Peierls phase, so `A > pi/2` already inverts the sign of the
modulated hopping. Deposited energy at 0.25 eV with the 1-cycle pulse:

| `PUMPAMPX` | deposited | transient max |
|---|---|---|
| 2.5 | 0.00122 eV | 0.32 eV |
| 5.0 | 0.00047 eV | 0.53 eV |
| 10.0 | 0.02424 eV | 0.70 eV |
| 20.0 | 3.49852 eV | 4.83 eV |

Non-monotonic between 2.5 and 5, then explosive by 20. There is no clean power law to fit, so
**always re-measure the deposited energy after changing the amplitude** (`TDOCCSCAN = True`, the
`E_stat` column) rather than scaling.

## To resolve the magnons in RIXS

The two low-energy excitations are 0.105 eV apart, so `epsloss` must be well under that -- the
default 0.1 will merge them. Use `epsloss ~ 0.02` and a loss grid fine enough to sample it
(`EM`/`NEDOS` giving <= 0.01 eV steps).

---

# Charge oscillations (PRX Fig. 2 style)

`plot_charge_oscillations.py` produces the Fig. 2 analogue for **one** pump configuration:
(a1) `A(t)`, (a2) `F[A(t)]`, (b1) e_g occupations, (c1) t2g + ligand occupations,
(b2) `F[<n>]`, (c2) total hole count as a conservation check. Input
`INPUT_charge_subcycle` (the envelope carried over unchanged from the earlier runs).
`INPUT_charge_1cycle` is kept for the envelope comparison but is a separate calculation, not
part of this figure.

**Occupations are HOLE occupations per orbital, summed over spin**, so each runs 0..2 and they
sum to the input `Holes`. NiO has `Holes = 2`:

| orbital | equilibrium `<n>` |
|---|---|
| d_x2-y2 | 0.95051 |
| d_z2 | 0.95981 |
| p_xx, p_yy | 0.02804 each |
| p_zz | 0.02149 |
| d_xy | 0.00614 |
| d_xz, d_yz | 0.00297 each |
| **sum** | **2.00000** (conserved to 1e-6 through the pump) |

The e_g pair carries ~0.95 holes each and everything else is small, which is why (b1) and (c1)
use split y-ranges rather than the paper's single 0..1 axis: the paper can use 0..1 because its
occupations swing by ~0.7, whereas NiO at 0.25 eV moves by ~0.02. Everything is in absolute hole
units -- an earlier version of this figure plotted `d<n> x 1e3`, whose +-20 axis reads like an
occupancy of 20 and was simply a bad choice.

Sampling: `TDDT = 0.00785` and a window of -20..400 eV^-1. `scan_td_occupation` hardcodes
`print_every = 20`, so the output step is `20*TDDT = 0.157` eV^-1 (Nyquist 20 eV), and the window
gives 0.015 eV resolution -- needed because the magnon period is 60 eV^-1. ~53000 propagation
steps.

## Result: the charge rings at the CT and dd scales, NOT at the magnon frequencies

Fraction of post-pulse (`t > 60`) spectral weight in each band:

| orbital | magnon 0.105 | 2-magnon 0.210 | dd ~1.0 | CT ~4.8 |
|---|---|---|---|---|
| **sub-cycle envelope** | | | | |
| d_x2-y2 | 0.13 % | 0.18 % | 34.4 % | **53.4 %** |
| d_z2 | 0.17 % | 0.22 % | 18.7 % | **65.1 %** |
| p_xx | 0.01 % | 0.24 % | 1.5 % | **94.4 %** |
| p_yy | 0.01 % | 0.35 % | 1.7 % | **94.5 %** |
| d_xy | 0.34 % | 0.05 % | **92.8 %** | 0.4 % |
| **1-cycle envelope** | | | | |
| d_x2-y2 | 1.42 % | 0.12 % | **69.1 %** | 23.9 % |
| d_z2 | 1.79 % | 0.14 % | **62.8 %** | 29.0 % |
| p_xx | 0.13 % | 0.27 % | 6.4 % | **89.1 %** |
| p_yy | 0.07 % | 0.37 % | 8.9 % | **86.6 %** |
| d_xy | 3.07 % | 0.08 % | **94.5 %** | 0.1 % |

- The **ligands ring at the charge-transfer scale** (~4.8 eV, 87-94 % of their weight), the
  **d_xy at the dd scale** (~1.0 eV, 93-95 %), and d_x2-y2 / d_z2 split between the two.
- **The magnon frequencies carry essentially nothing** (0.01-3.1 %). That is the expected
  consequence of the model's structure: the Peierls pump is a *charge* coupling, and the magnon
  splitting is set by the static `MF`, which `calc_ham` puts in `ham_static` -- the pump cannot
  modulate it. The longer 1-cycle pulse does couple slightly better to the low-energy sector
  (1.4-3.1 % vs 0.1-0.3 %), but it stays small.

Amplitudes: `p_xx` swings by up to `-25e-3` about an equilibrium value of `0.028`, i.e. the ligand
hole is nearly fully depleted at moments, while `d_x2-y2` moves `+20e-3` on `0.951` (2 %). Large
*relative* modulation of a small occupancy, small absolute charge transfer.

## A note on how the two envelopes excite

Panel (a2) shows the difference clearly. The **sub-cycle** pulse has a broad spectrum peaking near
0.7 eV with real weight past 1 eV, so it can excite the dd manifold *directly*. The **1-cycle**
pulse is narrow at 0.25 eV and has fallen eight decades by 1 eV -- yet its response still shows dd
and CT peaks. That excitation is therefore **nonlinear**, which is consistent with `PUMPAMPX = 5`
being a large Peierls phase rather than a perturbation.

---

# Amplitude: why A0 = 5 gave no oscillations, and what A0 = 25.3 does

## The pump is a hybridization modulation, and it enters as cos(A)

Peierls substitution puts `kx = pi/2 + A_x(t)`, so `sx = 2 sin(kx) = 2 cos(A_x(t))` and **every**
TM-ligand and ligand-ligand matrix element is multiplied by `cos(A(t))`. On-site energies (`CF`,
`MLCT`), Coulomb (`SC1`/`SC2`/`FG`), spin-orbit and `MF` are untouched. Two consequences:
`|cos| <= 1` means the pump can only *reduce* the hybridization magnitude or flip its sign, never
raise it; and the response is not linear in `A`.

| | peak &#124;A&#124; | cos(A) range | sign changes of t_eff in the pulse |
|---|---|---|---|
| A0 = 5 | 0.914 | +0.610 .. +1.000 | **0** |
| A0 = 25.3 | 4.627 | -1.000 .. +1.000 | **4** |

At `A0 = 5` the hopping is squeezed smoothly by at most 39 % with no sign reversal -- a monotonic
adiabatic pinch, which is why nothing oscillated. `A0 = 25.3` (chosen to match the PRX peak field
at this frequency, since a fixed `A0` gives a 5x weaker field at 0.25 eV than at pi eV) sweeps the
full range and reverses the hopping four times.

## What changes

| | A0 = 5 | A0 = 25.3 |
|---|---|---|
| deposited energy | 0.00012 eV | **3.5182 eV** |
| transient E excursion | 0.0278 eV | 4.5845 eV |
| d_x2-y2 peak-to-peak | 0.0187 | **0.6068** |
| p_xx peak-to-peak | 0.0249 | **0.8342** |
| d_x2-y2 post-pulse p-p | 0.0022 | 0.3526 |
| p_xx post-pulse p-p | 0.0018 | 0.6430 |

`p_xx` now runs from its equilibrium 0.028 up to 0.84 -- the ligand hole goes from 3 % to 42 % of
the cluster's two holes -- and the oscillations persist long after the pulse with a clear beat.
Hole number stays conserved to 2e-6.

## CAUTION: at this amplitude the spectral band labels stop being diagnostic

At `A0 = 25.3` the post-pulse weight near 0.20 eV rises to 28-41 %, which is tempting to read as
the two-magnon line at 0.2097 eV. **That reading is not justified here.** With 3.52 eV deposited,
the populated manifold is dense: among the 372 valence levels, **355 level pairs** have spacings
in [0.19, 0.23] eV, the spacing density is ~3000 pairs per eV, and the frequency resolution of a
340 eV^-1 window is 0.0185 eV -- about 56 pairs per resolution element. The measured peak sits at
0.203 eV, within one resolution element of the 2-magnon, but so are hundreds of unrelated pairs.

The band assignments *are* meaningful at `A0 = 5`, where only 1.2e-4 eV is deposited and the
dynamics stays in a small neighbourhood of the ground state. There the conclusion stands: the
charge rings at CT (~4.8 eV) and dd (~1.0 eV), and the magnon frequencies carry essentially
nothing.

## The useful amplitude window

`E_stat` is conserved after the pulse, so a short window suffices to scan this:

| `PUMPAMPX` | peak &#124;A&#124; | min cos(A) | E_dep | regime |
|---|---|---|---|---|
| 5 | 0.914 | +0.610 | 0.00012 eV | adiabatic, no real excitation |
| 8 | 1.463 | +0.108 | 0.00317 eV | adiabatic |
| **12** | **2.194** | **-0.584** | **0.1234 eV** | **within the magnon manifold** |
| 16 | 2.926 | -0.977 | 0.8508 eV | below the dd gap (0.9381) |
| 20 | 3.657 | -1.000 | 2.0717 eV | above the gap, dense manifold |
| 25.3 | 4.627 | -1.000 | 3.5182 eV | above the gap, dense manifold |

For magnon-specific physics `PUMPAMPX ~ 12` is the interesting setting: the hopping reverses sign
(so there is real driving, not an adiabatic squeeze) while the deposited energy 0.123 eV stays
inside the Sz ladder, below the 0.9381 eV dd gap. Then a feature at 0.105 or 0.210 eV can actually
be attributed to the magnon manifold.

# Hybridization: this parameter set vs the Fortran d8 set

| | NiO octahedral | PRX/Fortran d8 square-planar |
|---|---|---|
| geometry | 6 ligands, 9 p orbitals | 2 effective ligands, 6 p |
| `tpd` | 0.35 | 1.00 |
| tpd_z | `tpd`/sqrt(3) = 0.2021 | `tpd`*`tpdzr` = 0.25 |
| tpd_pi | `tpd`*`sigpi` = 0.0875 | `tpd`*`sigpi` = 0.225 |
| `tpp` | 0.0 | 0.0 |
| tpp_pi | not settable (`tpp_sig_pi` = 0) | `TPPPI` = 0.25 |
| Delta (`MLCT`) | 4.5 | 10.900004 |
| `EFFDEL` | True (effective Delta) | False (bare Delta) |
| equilibrium 2*`tpd` | 0.700 eV | 2.000 eV |
| **2*`tpd`/Delta** | **0.1556** | **0.1835** |
| ground-state ligand-hole character | **3.88 %** | **6.11 %** |

The bare hopping is 2.9x smaller in NiO, but Delta is 2.4x smaller too, so the dimensionless
mixing parameter comes out within 17 % and the ligand-hole character differs by only 1.6x -- not
the factor 3 the bare `tpd` suggests. Two asymmetries: the d8 set has a much harder pi channel
(tpd_pi 0.225 vs 0.0875, plus a live `TPPPI` that the octahedral geometry cannot reach), and the
two use different Delta conventions (`EFFDEL = True` vs `False`), so the `Delta` row is not
strictly like-for-like and the ratio inherits that caveat.
