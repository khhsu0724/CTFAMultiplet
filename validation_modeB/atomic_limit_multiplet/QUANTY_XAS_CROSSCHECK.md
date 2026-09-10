# Atomic-limit L-edge XAS: our code vs Quanty 0.81

Purpose: before attributing the equilibrium-XAS satellite discrepancy (Fortran puts
12.4% of the L3 weight 4-9 eV above the main line; we put 1.3%) to the charge-transfer
part of the Hamiltonian, verify that the *single-ion* part is exact. This checks four
things that had never been independently validated:

1. the d^8 initial-state multiplet (d-d Slater F2, F4),
2. the 2p^5 3d^9 final-state multiplet (core-valence F0, F2, G1, G3),
3. the 2p spin-orbit coupling,
4. the XAS dipole matrix elements.

## Setup

Atomic limit of the PRX d8 parameter set: no coordination (no ligands), no crystal
field, no hybridization.

| quantity        | value    | source                                              |
|-----------------|----------|-----------------------------------------------------|
| F0_dd           | 6.0      | `SC2 = 6.0 0 0.12857 0 0.025`                       |
| F2_dd           | 6.2999   | `0.12857 * 49`  (main.cpp:619)                      |
| F4_dd           | 11.0250  | `0.025 * 441`   (main.cpp:620)                      |
| F0, G1, F2, G3 (core-valence) | 6.0, 4.635, 6.18, 2.63 | `FG` (upstream order F0,G1,F2,G3) |
| zeta_2p         | 11.5     | `SO = 11.5 0`                                       |

Note F4/F2 = 1.75 is unusual for a free ion (~0.62), but the legacy Fortran computes
`rk(2,...)=f2*49` and `rk(4,...)=f4*441` from the same input numbers, so both codes use
identical integrals.  This crosscheck is about whether the two codes agree, not about
whether the parameter set is physical.

Ours: `ours_d8xas/INPUT`, run with `solver = 1` (direct pole summation, no Lanczos) and
`NEDOS = 50000` so the output is an exact pole list rather than a broadened spectrum.

Quanty: `quanty_crosscheck/q_dip.lua`, which diagonalizes both configurations and
evaluates `sum_i sum_{pol} |<Phi_f| T_pol |Psi_i>|^2` exactly over the 21-fold 3F ground
term and the three polarizations.  (`q_d8xas.lua` is the broadened-spectrum version;
it agrees too, but the pole comparison is cleaner -- see "Pitfalls" below.)

## Result

Both codes give the same 11 dipole-allowed poles out of the 12 final multiplet levels.
The nondegenerate level at 1.844 eV is dipole-forbidden in both (Quanty weight 1.6e-26).

|  E_f (eV) rel. lowest | deg | Quanty weight | our weight | ours/Quanty |
|-----------|-----|----------|----------|----------|
|  0.00000  |  5  | 0.035616 | 0.035616 | 1.000000 |
|  0.06035  |  9  | 0.142857 | 0.142857 | 1.000000 |
|  0.66791  |  7  | 0.217753 | 0.217753 | 1.000000 |
|  1.31394  |  5  | 0.118871 | 0.118871 | 1.000000 |
|  1.60712  |  3  | 0.047214 | 0.047214 | 1.000000 |
|  2.57588  |  7  | 0.100181 | 0.100181 | 1.000000 |
|  2.91010  |  3  | 0.027796 | 0.027796 | 1.000000 |
| 17.30625  |  5  | 0.068528 | 0.068528 | 1.000000 |
| 18.49912  |  5  | 0.078573 | 0.078573 | 1.000000 |
| 18.95831  |  7  | 0.104288 | 0.104288 | 1.000000 |
| 19.45059  |  3  | 0.058324 | 0.058324 | 1.000000 |

- max |E_f difference|      = **4e-05 eV**
- max |weight ratio - 1|    = **2e-06**
- L3/L2 branching ratio     = **2.228806 both**

The full 45-level d^8 multiplet and the full 60-level 2p^5 3d^9 multiplet also match
level-by-level with identical degeneracies:

  d^8:        0 (x21), 0.053541 (x9), 1.767847 (x5), 1.792833 (x9), 6.203527 (x1)
  2p^5 3d^9:  0 (x5), 0.060351 (x9), 0.667908 (x7), 1.313933 (x5), 1.607120 (x3),
              1.844127 (x1), 2.575881 (x7), 2.910101 (x3), 17.306252 (x5),
              18.499118 (x5), 18.958306 (x7), 19.450588 (x3)

## Conclusion

The single-ion Hamiltonian and the XAS dipole operator are **exact**.  The satellite
discrepancy against the Fortran cannot come from the multiplet integrals, the core
spin-orbit coupling, the exchange sign, the `FG` channel assignment, or the dipole
matrix elements.  It must live in the charge-transfer / ligand sector: the definition of
the charge-transfer energy (`MLCT`, `EFFDEL`), the p-d hopping normalization, or whether
the core-hole potential acts on the ligand orbitals.

This also independently re-confirms the upstream core-valence exchange sign
(`if (s.first == s.second) meG *= -1;` in `multiplet.cpp`), which an earlier session
wrongly "fixed" and then retracted.

## Pitfalls (for anyone repeating this)

- **Do not compare broadened spectra.** Our `epsab` is a Lorentzian HWHM; Quanty's
  `Gamma` is a FWHM, so `epsab = 0.2` corresponds to `Gamma = 0.4`. Even matched, the
  best agreement was 1.5% RMS / 9% peak, which is entirely lineshape, not physics.
- **Do not shrink `epsab` to sharpen the poles.** `ContFracExpan` runs Lanczos without
  reorthogonalization, so at `epsab = 0.004` the true poles break up into ~40 ghost
  peaks split by a few meV. Use `solver = 1` (direct pole summation) instead.
- Quanty spectra are `-Im(G)/pi`; the `.spec` file interleaves Re and Im per spectrum,
  and `CreateSpectra` with n operators and m states writes n*m spectra.
