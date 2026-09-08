# Validation data for THEORY_CROSSCHECK_PRB99.md

All runs use `validation_modeB/modeb_pumpon_v2_fixed/INPUT` as the base (each run's exact INPUT is
kept alongside its output), the toy d8 case, channel 0 only (`2Sz=-2`, the 55-dim ground-state
block).

Diagnostic environment variables (all no-ops unless set; see `src/tr_photon.hpp` / `src/photon.cpp`):

| run | env | purpose |
|---|---|---|
| `new_noflip`  | (none) | corrected `sum_l \|A_l\|^2`, baseline |
| `new_flip3`   | `TRRIXS2_FLIP_VAL=3` | same, one valence eigenvector sign-flipped |
| `coh_noflip`  | `TRRIXS2_COHERENT_FINAL=1` | as-was `\|sum_l A_l\|^2`, baseline |
| `coh_flip3`   | `TRRIXS2_COHERENT_FINAL=1 TRRIXS2_FLIP_VAL=3` | same, sign-flipped |
| `det5a`,`det5b` | (none) | two identical runs, run-to-run reproducibility |
| `poff_new`,`poff_coh` | `PUMPAMPX=0` | pump-off, both conventions, for the equilibrium check |
| `eq_fine`     | `RIXS_GS_INDEX=0` | equilibrium Kramers-Heisenberg, 61-point incident grid |
| `eq_gs0`      | `RIXS_GS_INDEX=0` | same on the 16-point grid -- **returned zeros at 15 of 16 incident energies**; the equilibrium BiCGS path has a grid-dependent failure worth investigating separately. Kept as the record of that. |
| `poffhr_new`,`poffhr_coh` | `TRRIXS2_WINDOW_MULT=2.865`, `sigma_pr=3` | under-windowed (4 sigma) attempt; kept as the record of the truncation artifact described in the writeup |

`sign-flip invariance`: diff `new_noflip` vs `new_flip3` (agrees to <=4e-3) against `coh_noflip`
vs `coh_flip3` (differs by 57-172%).
