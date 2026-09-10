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
