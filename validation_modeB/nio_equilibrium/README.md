# NiO equilibrium XAS/RIXS regression against a stored reference

Reference: `~/Documents/GradSchool/Devereaux Group/Herr_Jost/NiO_ED_tr_magnons/equlibrium_NiO/`
-- output from an **earlier version of this code** (same file naming), octahedral 1-site NiO,
`HFscale = 0.8`, `EFFDEL = True`, `MF = 0.11`, `CROSS = True`.

## Result

| variant | XAS_X vs reference | RIXS_X_Y vs reference |
|---|---|---|
| reference INPUT + `FG_KMAP_LEGACY=1 EXCH_SIGN_LEGACY=1` | **0.00e+00**, r = 1.00000 | 1.4e-4, r = 1.00000 |
| reference INPUT, current code | 0.72, r = 0.749 | 0.99, r = 0.463 |
| F2/G1 unswapped, current code | 0.55, r = 0.861 | 0.88, r = 0.623 |
| F2/G1 unswapped + `EXCH_SIGN_LEGACY=1` | **0.00e+00**, r = 1.00000 | 1.4e-4, r = 1.00000 |

`eig.txt`, `dos.txt` and `tmat.txt` are **bit-identical** in every variant: the valence sector,
octahedral geometry, `EFFDEL`, `MF`, `HFscale` and the `SC2` reduced/unreduced conversion all
reproduce exactly. Both `calc_CV` fixes are core-hole-only, as intended. (The 1e-4 RIXS residual is
BiCGStab's iterative tolerance, not a discrepancy.)

## THE IMPORTANT PART: this INPUT has F2 and G1 in each other's slots

`FG` is documented as `F^0, F^2, G^1, G^3`. The stored INPUT has

    FG = 6.0  5.783  7.720  3.290

but the Hartree-Fock values for Ni2+ 2p3d are **F2 = 7.721, G1 = 5.787, G3 = 3.291** -- so F2 and
G1 are swapped. That swap **exactly cancels** the k-channel bug fixed in d44d4f0 (which fed FG[1]
into the k=1 exchange slot and FG[2] into the k=2 Coulomb slot): unswapping the INPUT and
restoring the old exchange sign reproduces the reference bit-for-bit.

So someone compensated for the bug by swapping the input values. **Consequence: with the fix in,
every INPUT written under that workaround is now wrong and must be corrected**, e.g. for NiO

    FG = 6.0  7.720  5.783  3.290     # F^0, F^2, G^1, G^3

Note the two datasets in this project used *opposite* conventions, which is why the k-channel fix
improved one and appeared to break the other:

| dataset | `FG` as written | ordering | effect of the fix |
|---|---|---|---|
| PRX / toy (`modeb_pumpon_v2_fixed`) | `6.0 6.18 4.635 2.63` | correct (F2, G1, G3) | improves it -- equilibrium XAS r 0.83 -> 0.93 |
| NiO (`equlibrium_NiO`) | `6.0 5.783 7.720 3.290` | **swapped** | breaks it until the INPUT is corrected |

After `HFscale = 0.8`, NiO's unswapped integrals are F2 = 6.176, G1 = 4.626, G3 = 2.632 --
identical to the PRX Fortran's `cf2/cf3/cf4`. Both datasets describe the same core-valence physics.

## Exchange sign, on these parameters

The remaining difference after unswapping is entirely the exchange-sign fix (1579bdf). Verified
here on the NiO parameter set: in the atomic limit (`nio_atomic/`, hybridization, CF, SOC and MF
switched off) the 2p^5 3d^9 multiplet matches exact Slater-Condon theory to **2.67e-05**:

| term | deg | analytic | ours |
|---|---|---|---|
| 3D | 15 | +4.0651 | +4.0651 |
| 1F | 7 | +4.4701 | +4.4701 |
| 1D | 5 | +5.4645 | +5.4645 |
| 1P | 3 | +6.2500 | +6.2500 |
| 3P | 9 | +8.2204 | +8.2204 |
| 3F | 21 | +8.2357 | +8.2357 |

So the stored NiO reference carries the same Hund's-rule-violating exchange sign as the PRX
Fortran, and `nio_unsw/` (unswapped `FG`, fixed sign) is the physically correct calculation.
