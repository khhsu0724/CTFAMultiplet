# Validating the full trRIXS incident-energy map against PRX 15, 011012

Working assumption for this document: **the paper-1 Fortran reference is correct.** The goal is
either to make Mode B agree with it across the whole incident-energy map, or to prove our
implementation right and explain the reference's error. Status: **not yet resolved, but the
discrepancy is now localized to the core-hole (intermediate) sector.**

Companion documents: [`THEORY_CROSSCHECK_PRB99.md`](THEORY_CROSSCHECK_PRB99.md) (the final-state
sum derivation and fix), [`ANTISTOKES_FORTRAN_REFERENCE_CHECK.md`](ANTISTOKES_FORTRAN_REFERENCE_CHECK.md).
Data: `validation_modeB/prx_fig8_antistokes/`.

## 1. What is reproduced

The Hund's-exchange anti-Stokes feature of Fig. 8 comes out correctly:

| | reference | new code |
|---|---|---|
| anti-Stokes peak position | -0.9375 eV | **-0.9375 eV** (both XX and XY) |
| gain side before the pump | flat (contrast 1.00) | flat (contrast 1.00) |
| contrast after the pump, XX / XY | 1.87 / 5.48 | 5.75 / 5.66 |
| enhancement at the peak, XX / XY | 51.5x / 35.5x | 105x / 50x |
| XX:XY intensity ratio at the peak | 4.4x | 3.8x |

Normalized gain-side lineshapes are nearly superimposable (see `antistokes_cut.png`).

## 2. What is NOT reproduced: the pump-induced red shift of the absorption resonance

Peak and centroid of the map's incident-energy profile (summed over energy loss), vs delay:

| delay (eV^-1) | reference peak | reference centroid | new code peak | new code centroid |
|---|---|---|---|---|
| -6.28 | -4.51 | -4.52 | -4.75 | -4.65 |
| 0.00  | -4.63 | -4.60 | -4.95 | -5.00 |
| +3.9  | -4.51 | -5.08 | -4.95 | -5.04 |
| +8.8  | **-6.55** | -5.88 | -4.75 | -5.01 |
| +22.3 | **-6.55** | -6.08 | -4.75 | -5.00 |

The reference red-shifts its whole spectrum by **-2.0 eV** (setting in between delay +3.9 and +8.8)
and loses **2.8x** of its total spectral weight (0.049 -> 0.020 over its *full* incident range, so
this is not weight leaving a restricted window). We shift by -0.35 eV in centroid, 0.0 eV in peak,
and *gain* 1.5x in weight.

**This is not an artifact of the reference.** The red shift is visible in the published Fig. 7
itself: the dominant blob sits at Omega_i ~ 855 eV in panel (c) and ~852.7 eV in panel (d). The
paper, the reference code, and our code disagree in the same direction -- we are the outlier.

This single discrepancy explains the anti-Stokes incident-energy mismatch noted in
`ANTISTOKES_FORTRAN_REFERENCE_CHECK.md` (its finding 2): the reference's anti-Stokes sits 2.4 eV
below its equilibrium resonance because *everything* has red-shifted; ours sits 0.8 eV below
because nothing has.

## 3. Components verified -- do not re-check these

- **Valence Hamiltonian: verified to <0.05 eV.** Level pattern relative to each code's own ground
  state, with multiplicities: reference `0.0000 x3, 0.8338 x1, 1.0391 x3, 1.0617 x6, 1.5184 x1,
  1.8826 x3, 1.9694 x6, 2.5211 x1, 2.5903 x2`; ours `0.0000 x3, 0.8357 x1, 1.0653 x3, 1.0617 x6,
  1.5190 x1, 1.9122 x3, 1.9694 x6, 2.5702 x1, 2.5903 x2`. Multiplicities match exactly.
- **`MLCT = 10.900004` re-derived from the reference's own formula** `delxx = nd*aveu` with
  `aveu = raca - 14*racb/9 + 7*racc/9`, `nd=2` -> 10.900004. Exact.
- **Pump: exact.** Peak |A| ratio reference/ours = 2.0012 (the expected factor 2: their Peierls
  argument is `(k+A)/2`, ours `k/2+A`), same phase, same peak time to within a grid step.
- **Energy deposited by the pump: +7.5928 eV**, exactly conserved after the pulse (spread 0.0000
  eV over 24 eV^-1) -- the propagator is excellent.
- **Occupations** agree with the reference's own `Occ_Pulse.dat` to <=2.5% absolute with no drift.
- **Mode B's structure matches the reference's**, confirmed by reading its
  `tpivert`/`bige`/`tevolution` construction: both build a time-ordered product of pump-modulated
  single-step propagators in the core-hole sector, both propagate the initial state from before
  the pump, both use the same vertex factorization.
- **Time-step and window: converged.** dt=0.0392 and the reference's own dt=0.0196 give identical
  results (leakage 0.0007/0.0019/0.0098 vs 0.0007/0.0018/0.0097; peak, contrast, enhancement and
  anti-Stokes incident energy all unchanged).
- **The final-state-sum convention is NOT the cause**: running our code with the reference's
  coherent convention at converged settings puts the anti-Stokes incident maximum at -5.58,
  essentially where our corrected run puts it (-5.76), not at the reference's -7.39.

## 4. Where the discrepancy is: the core-hole sector

> **RESOLVED (see `validation_modeB/sqpl_cluster_quanty/README.md`).** The satellite
> discrepancy reported in sections 4 and 4b is **real** -- the reference is
> `Ni-L-RIXS/new_shape_function/XX/16_points_new_shape_function_XX/XAS_X.dat`, from the
> production run -- and it is now **fully explained**.
>
> Cause: the legacy Fortran gives the **ligand-ligand hopping the opposite sign in the
> core-hole sector to the valence sector** (its line ~1555 vs ~764). That is not a gauge
> choice; the p -> -p gauge on the ligands leaves that term invariant. Running our code with
> the same convention (`EX_HOP_SIGN=-1`) reproduces the reference XAS **exactly** -- all four
> window weights to four decimals and r = 1.0000, including w[+4..+9] = 0.1238 against our
> default's 0.0067. Since an O-O hopping cannot depend on whether a core hole sits on the
> metal, our default is the physical one and the reference's +1.9 eV "charge-transfer
> satellite" appears to be an artifact of that sign.
>
> Section 4's conclusion -- "an error in the core-hole Hamiltonian is the remaining
> candidate" -- was therefore correct. Its `EX_HOP_SIGN` test looked negative only because
> that diagnostic also flipped the ligand charge-transfer energy `del`, which lives on the
> diagonal of the same matrix (a 21.8 eV error); it has been fixed to flip only off-diagonals.
>
> Also confirmed since: the `FG` ordering correction is right; the "open structural
> difference" (a possible inter-site monopole in the Fortran) is **resolved, and the answer is
> no** -- its `checkint` gate enforces strictly on-site interactions, as our `calc_coulomb`
> does. And the square-planar cluster XAS is verified **exact** against Quanty 0.81 (120
> dipole-allowed poles to 4e-05 eV, every weight ratio 1.0000, sum rule to 1e-7), as is the
> atomic limit.
>
> Two real bugs in *our* code were found in the process: `SquarePlanar::get_tmat` gave the two
> inter-ligand hopping elements opposite signs where Slater-Koster requires them equal, and all
> three RIXS paths dropped the spin-flip channel (43-54% of the cross section). Both fixed;
> both backward compatible; NiO regression bit-identical.

The equilibrium XAS, with both spectra broadened identically (0.35 eV Gaussian), X polarization:

| | main L3 peak | satellite | L2 region | Pearson r vs ref |
|---|---|---|---|---|
| reference | -4.60 (1.00) | **+1.90 (0.29)** | 12.95, 13.95 | -- |
| ours | -5.10 (1.00) | **-3.00 (0.47)** | 12.00, 13.90 | 0.832 |

Our main L3 peak is 0.5 eV low, and -- more importantly -- **the reference has ~30% relative XAS
weight at +1.9 eV that we place at -3.0 eV instead**, a ~5 eV misplacement of a satellite. The
core-hole bandwidth agrees well (reference 27.4 eV, ours 26.9 eV from `cDOS.dat` vs `exdos.txt`)
but the spectral centroid above the lowest core-hole level differs by ~2 eV (9.08 vs 11.08).

Since the valence Hamiltonian, the pump, and the propagator are all verified, and the RIXS
assembly structure matches, **an error in the core-hole Hamiltonian (or in the dipole matrix
elements into it) is the remaining candidate**, and it is the right size and location to produce
a ~2 eV error in the pumped absorption resonance.

### Tested and rejected as the cause

The reference builds the core-hole sector's hopping with the **opposite overall sign** to the
valence sector -- its line ~1552 (`+2*tmatx`, `+2*tmaty`, `+4*tmato`) against line ~761
(`-2*tmatx`, `-2*tmaty`, `-4*tmato`), the two code blocks otherwise byte-identical. Flipping all
three together is a `p -> -p` gauge transformation only if the p-p term (`tmato`, the `tpp`
channel) vanishes; with `tpp_pi = 0.25` it does not, so these are genuinely different core-hole
Hamiltonians. Implemented as the `EX_HOP_SIGN=-1` diagnostic in `calc_HYB` (off by default; the
default path is bit-identical, verified). Result: it moves our XAS slightly toward the reference
(main peak -5.10 -> -4.95, satellite -3.00 -> -1.85) but does not reproduce it, and the overall
agreement is unchanged (r 0.832 -> 0.834). **Real convention difference, not the explanation.**

## 4b. Rung 2 of the ladder: equilibrium XAS -- one real bug found and fixed

### CORRECTED: this was an INPUT translation error, not a code bug

**The "k-channel bug" described below was a misreading and the code change has been reverted.**
Upstream CTFAMultiplet documents `FG` as `F0,G1,F2,G3`, and `calc_U`'s indexing by multipole
order already matched that convention. What was actually wrong is
`validation_modeB/modeb_pumpon_v2_fixed/INPUT`: translating the Fortran's `cf2`(F2)=6.18,
`cf3`(G1)=4.635 into our order requires `FG = 6.0 4.635 6.18 2.63`, but it was written
`FG = 6.0 6.18 4.635 2.63`. Correcting the INPUT gives numerically the same improvement the
code change did (r 0.832 -> 0.932, peak -5.10 -> -4.600), because swapping the inputs and
swapping the k-channels are the same operation. The section below is kept for the reasoning
trail; read "fix" as "corrected FG ordering in the INPUT".

### Original (superseded) description

`calc_CV` passed the `FG` array straight into `calc_U`, but **`calc_U` contracts its Slater array
by multipole order k** (`u += gaunt1[k]*gaunt2[k]*SC[k]`). `FG` holds `{F^0, F^2, G^1, G^3}` in
input order, and those four integrals belong at k = 0, 2, 1, 3. So `F^2` was applied at k=1 and
`G^1` at k=2: the Coulomb channel was fed `G^1 = 4.635` where it needed `F^2 = 6.18`, and the
exchange channel was fed `F^2` where it needed `G^1`.

It was silent because the Gaunt parity rules kill the mismatched slots -- for `l_v=2, l_c=1` the
Coulomb product is nonzero only for even k and the exchange only for odd k -- so the code simply
used the wrong two numbers with no warning.

**Fixed** in `src/multiplet.cpp` by giving each integral its own k channel
(`Fk = {F0, 0, F2}`, `Gk = {0, G1, 0, G3}`). Verified: the source fix reproduces a
swapped-input control run to `max|diff| = 0.000e+00`.

| | r vs reference | main L3 peak |
|---|---|---|
| reference | 1.0000 | -4.633 |
| before fix | 0.8318 | -5.10 |
| **after fix** | **0.9315** | **-4.600** |

**This affects every L-edge result this repository has produced**, not just Mode B. The valence
(GS) spectrum is bit-identical (`calc_CV` returns early unless `is_ex`), so the <0.05 eV valence
agreement of section 3 still stands. The K edge is unaffected: the code requires `FG[2]=FG[3]=0`
there, so `Gk` is all zeros.

### What remains: ~11% of the XAS weight is misplaced

**No dipole weight is missing.** The sum rule agrees to 0.1%: reference integrated XAS/pi/3 =
0.2429, ours 0.2426 (and the same normalisation recovers 131.1 ~ 132 core states from the
reference's unweighted DOS column, confirming the convention). It is purely a redistribution:

| window rel. to main L3 peak | reference | ours |
|---|---|---|
| -3 .. +1 | 0.428 | 0.495 |
| +1 .. +4 | 0.133 | 0.165 |
| **+4 .. +9** | **0.124** | **0.013** |
| +14 .. +22 (L2) | 0.251 | 0.301 |

The reference puts 12.4% of the total weight 4-9 eV above the main line -- its `+1.90 eV`
satellite, where a Ni L3 charge-transfer satellite belongs. We put 1.3% there.

The core-hole spectra themselves also differ: comparing cumulative distributions (dipole-free,
each aligned at its own lowest level) gives max |CDF difference| = 0.31. Both codes put exactly
2/3 of the 132 core states in the L3 manifold and 1/3 in L2, so **core SOC is right in both**, but
our L3 manifold has extra low-lying states.

### Ruled out as the cause of the missing satellite

All tested with the FG fix in place, tracking `w(+4..+9)` (target 0.124):

- `MLCT` 6 -> 16 eV: `w` stays 0.007-0.013. Nothing tracks the charge-transfer energy. (MLCT =
  10.9 does give the best overall agreement of the scan, r = 0.9315.)
- `tpd` 0.5/2.0/3.0, `sigpi` 0.5, `tpdzr` 1.0, `SO` 6.0/17.0: `w` never exceeds 0.05, and all
  degrade r.
- Uniform scaling of the core-valence multipoles (F2, G1, G3) by 0.5-4x: generates spurious peaks,
  r collapses to 0.03-0.87.
- `FG[0]` (core-valence monopole) 3 -> 9 eV: shifts the whole L3 region **rigidly** 1:1 and leaves
  the main-line-to-satellite splitting at ~2 eV. So the visible companion peak is a d-hole state,
  not the charge-transfer satellite.
- `EX_HOP_SIGN=-1` (the reference's opposite core-hole hopping sign, section 4): combined with the
  FG fix it makes agreement *worse* (r 0.9315 -> 0.8949) and leaves `w` at 0.0127.
- Missing dipole channels: excluded by the sum rule above.
- Our XAS summing over only one ground-state sublevel: checked, it loops over all of `gsi`, same
  as the reference's `do qqq / gsec(ws)==1 / do j=1,gdegen(ws)`.

### The open structural difference

Our `calc_coulomb` uses a **single site index for all four quantum numbers** -- it is strictly
on-site (metal d-d and ligand p-p within one ligand). The reference's interaction tensor `hint` is
indexed by `(l, m_l)` only, with no site label, so `rk(0,1,2,2,1) = cf0 = 6.0` nominally applies to
*any* l=1/l=2 hole pair -- including **ligand-2p <-> metal-3d**, an inter-site monopole we do not
have at all. That is exactly the term that sets where the charge-transfer configuration sits, and
it acts differently in the two sectors because the valence hole count differs (2 in GS, 1 in EX).

**Caveat before acting on this:** if the reference really applied a 6 eV inter-site monopole, its
d^1 L^1 configurations would sit ~6 eV higher than ours and the *valence* spectra could not agree
to <0.05 eV -- but they do. So either the reference's `mflag`/`mdiff` momentum bookkeeping excludes
inter-site pairs in practice, or the term is absorbed elsewhere. **Resolving this requires tracing
`lflag`/`mflag`/`mdiff` and the `hpairs`/`hothers` hole-pair enumeration** to see which site pairs
actually reach `hint`. That is the next concrete step, and it decides whether we are missing a term
or the reference is applying one it should not.

## 4c. Rung 3-4: did the fix improve the match? Yes for XAS; equilibrium RIXS has no usable reference

### XAS: clearly improved

r 0.8318 -> 0.9315, and the main L3 peak moved from -5.10 to **-4.600** against the reference's
-4.633 (0.03 eV). Sum rule agrees to 0.1%. Remaining defect: 12.4% of the weight belongs 4-9 eV
above the main line and we put 1.3% there (section 4b).

### Equilibrium RIXS: improved on the incident axis, and validated internally

**There is no usable equilibrium-RIXS reference in this dataset.** The only files with matching
parameters are `new_shape_function/YY/16_points_.../SPEC_*.dat`, and that output is unphysical: the
energy-loss profile decays monotonically across the whole range with **no elastic peak at loss=0
and no dd features at all** --

```
loss=-1.922 I=0.0216   loss=-0.984 I=0.0105   loss=-0.047 I=0.00622
loss=+1.047 I=0.00385  loss=+8.000 I=0.000696
```

-- i.e. a bare Lorentzian tail from poles placed far below the grid. Its
`omega - (deig(wsp,i) - a(1))` energy reference is evidently wrong, so `E_f - E_gs` comes out
large and negative. (`SPEC_XX_10.dat` in the same directory is identically zero in all 256x128
cells.) Consistent with the picture elsewhere: the paper's figures come from the `trSPEC`
(time-domain) path, and this equilibrium block was not the one used.

So equilibrium RIXS is validated **internally** instead, which is possible because the valence
spectrum is independently verified against the reference to <0.05 eV (section 3):

1. **Loss axis.** Our equilibrium RIXS loss poles must sit at the validated valence excitation
   energies `E_f - E_gs`. They do, over six peaks:

   | RIXS loss peak | nearest validated level | offset |
   |---|---|---|
   | 0.000 | 0.000 | +0.000 |
   | 1.062 | 1.062 | +0.001 |
   | 1.938 | 1.912 | +0.025 |
   | 2.562 | 2.570 | -0.008 |
   | 3.000 | 3.005 | -0.005 |
   | 3.625 | 3.620 | +0.005 |

   max |offset| = 0.025 eV. The loss axis was already correct and the `calc_CV` fix does not touch
   it (`calc_CV` is `is_ex`-only, so final-state energies are unaffected).

2. **Incident axis.** The RIXS incident profile peaks at **-4.633** -- exactly the reference's XAS
   resonance -- and tracks our own XAS with r = 0.8716.

3. **Effect of the fix.** The equilibrium RIXS absorption resonance moved from **-5.231 to -4.633**,
   i.e. onto the reference's resonance. The 2D map changed substantially (Pearson r between the
   pre-fix and post-fix maps is only 0.48), so the bug was materially distorting every equilibrium
   RIXS result, entirely through the incident axis.

`FG_KMAP_LEGACY=1` restores the pre-fix mapping for exactly this kind of before/after measurement.

## 5. Next steps, in order

1. **The core-hole attraction on the valence/ligand levels.** The reference passes the
   core-valence monopole into its interaction tensor via `gethint(f0,f2,f4,cf0,cf2,cf3,cf4,...)`
   with `cf0 = f0 = 6.0`. A charge-transfer satellite moves by several eV depending on whether,
   and how, that monopole shifts the ligand level relative to the metal in the *intermediate*
   state. Check how our `FG[0]`/`SC2EX`/`MLCT` handling reproduces that shift. This is the
   leading candidate: it is exactly a several-eV effect on a charge-transfer satellite.
2. Line-by-line comparison of the core-valence Coulomb/exchange assembly (`gethint` vs our
   `calc_coulomb`/Gaunt path) for the `cf2/cf3/cf4` terms.
3. Only after the equilibrium XAS matches: re-run the maps and re-check the pumped red shift.
   The equilibrium XAS is the right regression target because it is seconds to compute and needs
   no pump.
