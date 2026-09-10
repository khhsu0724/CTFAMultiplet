# Atomic-limit validation of the core-valence interaction (calc_CV)

A reference-independent test. Switching off hybridization, crystal field and core spin-orbit
coupling leaves the core-hole sector as the pure **2p^5 3d^9 multiplet**, whose term energies are
known exactly from Slater-Condon theory. `pdterms.py` (in `validation_modeB/`) computes them from
Wigner 3j/6j symbols via Cowan eq. 12.5-12.6; its symbols are self-checked against known values.

Run any of these INPUTs and compare `exdos.txt` against `pdterms.pd_terms(F0,F2,G1,G3)`.

| directory | exchange sign | what it shows |
|---|---|---|
| `atomic`, `at_onlyG`, `at_noF2` | pre-fix | every term energy exactly right but **singlet and triplet of each L interchanged**. NOTE: these were run with the since-reverted k-map change and an INPUT whose F2/G1 were swapped; `at2_leg`/`at2_fix` below repeat the test cleanly |
| `at_noG` | n/a (G=0) | direct channel alone: singlet/triplet exactly degenerate (deg 20, 28, 12) -- the direct part was always correct |
| `x_atomic`, `x_at_onlyG`, `x_at_noF2` | fixed | matches exact theory to **5e-5** |
| `fortran_atomic_limit` | -- | the paper-1 Fortran run in the same limit |

## Result

Degeneracies of the 2p^5 3d^9 terms are 1P=3, 1D=5, 1F=7, 3P=9, 3D=15, 3F=21 (60 states; the
remaining 72 are the ligand-hole configurations, degenerate at MLCT).

Our fixed code reproduces the exact analytic energies **and** assignments. The pre-fix code
reproduced the energies with singlet/triplet swapped -- the signature of a sign error on the
core-valence exchange, since flipping `g_k = -(-1)^S ...` interchanges S=0 and S=1 at every L.

**The paper-1 Fortran reference has the same sign error.** Integrating its own `cDOS.dat` from an
atomic-limit run (`fortran_atomic_limit/`, γ=0.01, areas normalised on the 72-state ligand peak):

| level | degeneracy in the reference | correct sign predicts | opposite sign predicts |
|---|---|---|---|
| +4.062 | 5.4 | 3D (15) | **1D (5)** |
| +4.467 | 21.1 | 1F (7) | **3F (21)** |
| +5.466 | 15.1 | 1D (5) | **3D (15)** |
| +6.251 | 9.1 | 1P (3) | **3P (9)** |

Four independently resolved levels, all matching the wrong sign. (The 8.221/8.239 pair is split by
0.018 eV, below the resolution of the γ=0.01 Lorentzians, so it only confirms 3+7=10 either way.)

So the reference's 2p-3d exchange favours **low** spin in the core-hole multiplet -- a Hund's-rule
violation -- and that propagates into its XAS, trXAS, RIXS and trRIXS, including PRX Fig. 7/8.
The upstream sign is the **default**; `EXCH_SIGN_FIX=1` opts into the corrected one.

**Caveat worth raising with the paper's authors:** the values `cf2=6.18, cf3=4.635, cf4=2.63` may
have been chosen to match experimental Ni L-edge spectra *under* that sign. With the corrected sign
the same numbers give a different spectrum, so reproducing the published figures and doing correct
atomic physics are, for now, two different things.

## Re-done after reverting the k-map change (authoritative)

The k-channel change referenced above was **wrong and has been reverted** -- upstream documents
`FG` as `F0,G1,F2,G3` and `calc_U`'s indexing already matched that. The exchange-sign conclusion is
unaffected, and was re-established from scratch using the NiO INPUT exactly as written
(F2 = 6.176, G1 = 4.626, G3 = 2.632 after `HFscale = 0.8`):

| exchange sign | max |diff| vs exact Slater-Condon theory |
|---|---|
| upstream (**default**, no flag) | **3.77 eV**, singlet and triplet of each L interchanged |
| fixed (`EXCH_SIGN_FIX=1`) | **2.67e-05** |

Term by term with the fix: 3D +4.0651, 1F +4.4701, 1D +5.4645, 1P +6.2500, 3P +8.2204,
3F +8.2357 -- all exact. So the exchange sign was genuinely wrong, independently of the `FG`
ordering question, and both the stored NiO reference and the PRX Fortran carry it.
