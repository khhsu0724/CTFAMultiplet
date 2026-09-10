# Atomic-limit validation of the core-valence interaction, settled against Quanty

Switching off hybridization, crystal field and core spin-orbit coupling leaves the core-hole
sector as the pure **2p^5 3d^9 multiplet**. That limit is directly comparable to Quanty
(Haverkort et al.), which is the arbiter used here.

## Result: the code's core-valence interaction is CORRECT as it stands

`quanty_crosscheck/q_pd.lua` builds the same 60-state multiplet in Quanty 0.81 with only the
2p-3d Coulomb interaction and the same unreduced Slater integrals (F0 = 6.0, F2 = 6.1760,
G1 = 4.6264, G3 = 2.6320 -- the NiO values after `HFscale = 0.8`):

| energy above lowest | degeneracy, Quanty | degeneracy, this code |
|---|---|---|
| 0.0000 | 5 | 5 |
| 0.4050 | 21 | 21 |
| 1.3994 | 15 | 15 |
| 2.1849 | 9 | 9 |
| 4.1553 | 3 | 3 |
| 4.1706 | 7 | 7 |

Identical degeneracy sequence, and the energies agree to **4.7e-05 eV**. No flags, no changes:
`calc_CV`'s exchange sign and the `FG = F0,G1,F2,G3` ordering are both right.

## Retracted: the "exchange-sign bug"

An earlier revision of this branch claimed `calc_CV` had an inverted core-valence exchange sign,
and briefly shipped a flag to flip it. **That was wrong and is fully retracted** -- the flag has
been removed and the code restored to the upstream expression unconditionally.

The error was in this directory's own analytic tool. `validation_modeB/pdterms.py` implemented the
Slater-Condon exchange coefficient as `g_k = -(-1)^S ...` when it should be `+(-1)^S ...`. That one
sign interchanges the singlet and triplet of every L, which produced a self-consistent but wrong
"exact theory" that the real code appeared to contradict. `pdterms.py` is now corrected and
reproduces Quanty's sequence and energies.

Two things should have raised suspicion earlier and did not:

- The supposed bug made agreement with **two independent references worse** (the stored NiO
  equilibrium run and the PRX Fortran), which is the wrong direction for a genuine fix.
- The Hund's-rule argument used to support it was misapplied. For two holes in *inequivalent*
  shells with a large F2, the lowest term genuinely need not be high-spin -- here it is 1D
  (degeneracy 5), as both Quanty and this code find.

Consequently the claim that the PRX Fortran and the stored NiO reference contain a
"Hund's-rule-violating exchange sign" is **also retracted**. Their core-valence exchange is correct.

## Files

| directory | what it is |
|---|---|
| `quanty_crosscheck/q_pd.lua`, `.out` | the Quanty input and output used to settle this |
| `at2_leg/` | this code in the same atomic limit -- matches Quanty |
| `at2_fix/` | the same with the retracted sign flip -- singlet/triplet interchanged, kept as the record |
| `atomic`, `at_noG`, `at_onlyG`, `at_noF2`, `x_*` | earlier runs from the retracted analysis; `at_noG` is still a useful check that the direct channel alone gives exactly degenerate singlet/triplet pairs (20, 28, 12) |

## What still stands from the atomic-limit work

The degeneracy structure itself is a genuine validation: this code produces exactly the
1P/1D/1F/3P/3D/3F terms with degeneracies 3, 5, 7, 9, 15, 21 (60 states) plus the ligand-hole
manifold, and with the exchange integrals switched off the singlet/triplet pairs are exactly
degenerate (20, 28, 12). Both are correct and were never in question.
