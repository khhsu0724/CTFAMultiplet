# NiO equilibrium XAS/RIXS regression against a stored reference

Reference: `~/Documents/GradSchool/Devereaux Group/Herr_Jost/NiO_ED_tr_magnons/equlibrium_NiO/`
-- output from an **earlier version of this code** (same file naming), octahedral 1-site NiO,
`HFscale = 0.8`, `EFFDEL = True`, `MF = 0.11`, `CROSS = True`.

## Result

| variant | XAS_X vs reference | RIXS_X_Y vs reference |
|---|---|---|
| reference INPUT + `FG_KMAP_LEGACY=1 EXCH_SIGN_FIX=1` (both flags since removed) | **0.00e+00**, r = 1.00000 | 1.4e-4, r = 1.00000 |
| reference INPUT, current code | 0.72, r = 0.749 | 0.99, r = 0.463 |
| F2/G1 unswapped, current code | 0.55, r = 0.861 | 0.88, r = 0.623 |
| F2/G1 unswapped + `EXCH_SIGN_FIX=1` (flag since removed) | **0.00e+00**, r = 1.00000 | 1.4e-4, r = 1.00000 |

`eig.txt`, `dos.txt` and `tmat.txt` are **bit-identical** in every variant: the valence sector,
octahedral geometry, `EFFDEL`, `MF`, `HFscale` and the `SC2` reduced/unreduced conversion all
reproduce exactly. Both `calc_CV` fixes are core-hole-only, as intended. (The 1e-4 RIXS residual is
BiCGStab's iterative tolerance, not a discrepancy.)

## CORRECTION: this INPUT is CORRECT as written

An earlier version of this file claimed the stored INPUT had F2 and G1 swapped. **That was wrong.**
Upstream CTFAMultiplet documents `FG` as **"F0,G1,F2,G3"**
(https://github.com/khhsu0724/CTFAMultiplet/wiki/Input-Parameters), so

    FG = 6.0  5.783  7.720  3.290     # F0, G1, F2, G3

means F2 = 7.720 and G1 = 5.783, which after `HFscale = 0.8` gives F2 = 6.176, G1 = 4.626,
G3 = 2.632 -- against Hartree-Fock for Ni2+ 2p3d x0.8 of 6.177, 4.630, 2.633. The INPUT is right,
and `calc_U`'s indexing by multipole order k was already right for that ordering: the Coulomb
channel reads FG[0] at k=0 and FG[2] at k=2, the exchange channel FG[1] at k=1 and FG[3] at k=3,
with the mismatched slots killed by Gaunt parity. The `Udp = F0 - G1/15 - 3G3/70` printout in
main.cpp encodes the same convention.

The "k-channel fix" that briefly lived in this branch has therefore been **reverted**: inputs
written for the upstream code work unchanged, as they must.

## Result with the k-map change reverted

| variant | XAS_X | RIXS_X_Y | eig/dos/exeig/exdos |
|---|---|---|---|
| INPUT as written, upstream exchange sign | **0.00e+00**, r = 1.00000 | 1.4e-4, r = 1.00000 | **all IDENTICAL** |
| INPUT as written, exchange-sign fix | 0.55, r = 0.861 | 0.88, r = 0.623 | valence identical |

So the code reproduces this NiO reference **completely** -- valence sector, core-hole sector, XAS
and RIXS -- and the only behavioural change from this branch is the exchange-sign fix.

## Translating Fortran parameters into this convention

The Fortran declares `cf2` = F^2, `cf3` = G^1, `cf4` = G^3, i.e. **F0, F2, G1, G3** order. Our
`FG` is **F0, G1, F2, G3**, so F2 and G1 must be exchanged when translating:

    Fortran: cf0=6.0  cf2=6.18 (F2)  cf3=4.635 (G1)  cf4=2.63 (G3)
    ours   : FG = 6.0  4.635  6.18  2.63

`validation_modeB/modeb_pumpon_v2_fixed/INPUT` (and every INPUT derived from it during the PRX
work) has `FG = 6.0 6.18 4.635 2.63`, which under our convention means G1 = 6.18, F2 = 4.635 --
i.e. **that translation was done wrong**. Correcting it improves the equilibrium XAS against the
Fortran from r = 0.832 to r = 0.932 with the main L3 peak moving from -5.10 to -4.600 (reference
-4.633). Use `FG = 6.0 4.635 6.18 2.63` for that dataset.

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

*(The paragraph that stood here, concluding that the stored reference "carries the same
Hund's-rule-violating exchange sign as the PRX Fortran", is retracted -- see the RETRACTED section
below. The analytic table above is from `pdterms.py`, which had the exchange sign wrong at the
time; it has since been corrected and now agrees with Quanty and with the code.)*

## Default behaviour (final)

**The upstream exchange sign is the default.** An input written for the previous code, run with no
environment variables set, reproduces its previous results:

| output | default run vs stored reference |
|---|---|
| `eig.txt`, `dos.txt`, `tmat.txt` | **bit-identical** |
| `exeig.txt`, `exdos.txt` | **bit-identical** |
| `XAS_Ledge_X/Y/Z.txt` | **bit-identical** |
| `RIXS_Ledge_*.txt` | max rel 1e-4 ... 1.3e-3 |

The RIXS residual is BiCGStab's own convergence tolerance, not a discrepancy: `CG_tol = 1e-8` is
tested against the **squared** residual norm, so it converges to `||r|| ~ 1e-4`. Iteration counts
match (7 per solve) and the maps correlate at r = 1.0000000 with median difference exactly 0.

Related caveat: RIXS is not bit-reproducible across thread counts for the same reason -- runs at
3/6/9 threads agree to ~1e-4 but not exactly. XAS and all eigen/DOS output *are* bit-reproducible.
Making RIXS bit-reproducible would mean tightening `CG_tol` and serialising the remaining
elementwise loops in `BiCGS`; the physics is unaffected at this level, since spectral features are
compared at the percent level.

## RETRACTED: the "exchange sign" discussion

An earlier version of this file said the stored reference "carries the same Hund's-rule-violating
exchange sign as the PRX Fortran" and pointed at an `EXCH_SIGN_FIX=1` flag. **Both claims are
withdrawn and the flag is gone.**

The upstream convention in `multiplet.cpp` (`if (s.first == s.second) meG *= -1;`) is **correct**.
It was verified against Quanty 0.81: the 2p^5 3d^9 multiplet on the NiO integrals reproduces
Quanty's degeneracy sequence 5, 21, 15, 9, 3, 7 and its energies to 4.7e-05 eV, and the full
atomic d^8 L-edge XAS agrees pole-for-pole to 4e-05 eV with every dipole weight matching to 2e-06
(`../atomic_limit_multiplet/QUANTY_XAS_CROSSCHECK.md`). The "Hund's-rule" argument behind the
retracted claim was misapplied: for two holes in *inequivalent* shells with a large F^2, the
lowest term need not be high-spin, and here it is genuinely 1D.

## Also verified here: the spin-flip RIXS fix does not touch this reference

`SO = 11.507 0.083 0.102` turns on **valence** spin-orbit coupling and `MF = 0.11` adds a mean
field, so this system's ground sector is a **single** 168-dimensional block rather than being split
into Sz blocks. Every RIXS final state was therefore already included, and the spin-flip fix
(`../rixs_spinflip_channel/`) leaves this reference reproduced identically -- `RIXS_Ledge_X_Y`
totals 172.896567 before and after, max|diff| 1.4e-4 against the reference in both cases (BiCGStab
tolerance). The bug only bites Sz-blocked ground sectors, i.e. valence SOC off and no field.
