# Square-planar cluster vs Quanty 0.81: our code is exact, and one real bug found

This closes out the equilibrium-XAS dispute against the legacy Fortran (PRX 15, 011012)
for the `PRX_d8_sqpl_Ledge` parameter set. Three results:

1. **Our full-cluster L-edge XAS is exact against Quanty** -- one-body, d-d, ligand p-p,
   core-valence, core spin-orbit and the dipole operator, all of it.
2. **A real bug in `SquarePlanar::get_tmat`**: the two surviving inter-ligand hopping
   elements were given *opposite* signs. Slater-Koster says they must be equal. Fixed.
   Backward compatible -- it only bites when `TPPPI != 0`.
3. **The "12.4% XAS satellite" discrepancy is retracted**: both Fortran reference runs
   set `flag_skip_XAS = 1`, so neither ever computed an XAS. There is no reference.

## Method

Quanty is set up in the **hole language**, mirroring our convention directly: 2 particles
in the 22 valence spin-orbitals (5 d + 6 ligand p), 1 particle in the 6 core 2p
spin-orbitals. In that language the core spin-orbit coupling enters as `-zeta*ldots`
(verified: it reproduces the electron-language 2p^5 3d^9 multiplet and our `exdos.txt`
level-for-level, see `../atomic_limit_multiplet/`).

Our 11x11 one-body matrix is taken verbatim from our own `tmat.txt` plus the `CF`
diagonal, rotated from the real-harmonic basis into Quanty's complex-harmonic basis with
the *same* `U_d`/`U_p` the code uses (`src/cluster.hpp`), and handed to Quanty as
`NewOperator("Number", NF, i_list, j_list, v_list)`. So the comparison isolates the
many-body assembly: identical one-body input, identical Slater integrals, independent code.

Our code is run with `solver = 1` (direct pole summation) and `NEDOS = 40000` so the
output is a pole list, not a broadened spectrum. Quanty evaluates
`sum_i |<Phi_f| T_x |Psi_i>|^2 / 3` over the 3-fold degenerate ground state, exactly.

## Result 1: exact agreement

**Valence levels** (2 valence holes, no core hole), ours vs Quanty, absolute energies:

| Quanty | ours | deg |   | Quanty | ours | deg |
|---|---|---|---|---|---|---|
| 3.674553 | 3.67455 | 3 | | 6.244690 | 6.24474 | 1 |
| 4.508540 | 4.51023 | 1 | | 6.264840 | 6.26488 | 2 |
| 4.736300 | 4.73620 | 6 | | 6.620370 | 6.62065 | 1 |
| 4.739960 | 4.73983 | 3 | | 6.653490 | 6.65348 | 6 |
| 5.192960 | 5.19359 | 1 | | 6.679370 | 6.67937 | 3 |
| 5.586720 | 5.58669 | 3 | | 7.294870 | 7.29488 | 1 |
| 5.643950 | 5.64394 | 6 | | 7.336250 | 7.33626 | 2 |

max difference 0.0017 eV, on the singly-degenerate levels only -- that is our block
Lanczos truncation, not a Hamiltonian difference.

**XAS**: 132 core-hole states, 120 dipole-allowed in both codes.

- max |pole energy difference| = **4e-05 eV**
- every pole weight ratio = **1.0000** (the apparent outliers are Quanty resolving
  degenerate pairs that our binned output merges; they sum exactly -- e.g. Quanty's
  1.210e-2 + 1.358e-2 at -4.67498 equals our single 2.568e-2 bin)
- dipole sum rule: ours 0.24429662, Quanty 0.24429665, **ratio 0.9999999**
- window-integrated weight agrees to **1e-7** in all six coarse windows

## Result 2: the inter-ligand hopping sign bug

Slater-Koster for the four O_x-O_y neighbour vectors (+-1/2, +-1/2), Bloch-summed:

    <O_x p_x| H |O_y p_y> = <O_x p_y| H |O_y p_x> = -2(Vsig - Vpi) sin(kx/2) sin(ky/2)

and the diagonal channels [xx<->yx], [xy<->yy], [xz<->yz] land on cos(kx/2)cos(ky/2),
which **vanishes** at the (pi,pi) momentum this geometry uses. So exactly two elements
survive, and they carry the **same** coefficient. Numerically, at k=(pi,pi) with
`vpps = 0`, `vppp = 0.25`: `(px,py) = (py,px) = +0.5`, all diagonal channels 0.

That also retires the old comment in `get_tmat` claiming the diagonal channels were
dropped "to keep the Hamiltonian real": they are dropped because they are zero here.

What the two codes had, as the ligand-x -> ligand-y block in the real basis:

|            | (px,py) | (py,px) | structure |
|---|---|---|---|
| Slater-Koster | +0.5 | +0.5 | symmetric |
| legacy Fortran | -1.0i | -1.0i | symmetric (correct, up to a gauge phase and a factor 2 the Fortran took on deliberately -- see below) |
| **ours, before** | **-1.0** | **+1.0** | **antisymmetric -- WRONG** |
| ours, after | -1.0 | -1.0 | symmetric |

The antisymmetry is not removable by a gauge phase: before the fix our one-body
eigenvalues were `0.96297, 9.9, 11.93703` where the Fortran has `0.954725, 9.945281,
11.900006`. After the fix the two one-body matrices agree to **6e-06**, which is exactly
our INPUT's `MLCT = 10.900004` against the Fortran's `nd*aveu = 10.900006`.

Equivalently, in the complex-harmonic basis: the correct block is *off*-diagonal in m
(m=-1 <-> m=+1), ours was *diagonal* in m. Antisymmetric in m <=> symmetric in the real
basis, because p_y carries a factor i.

The fix is one character (`tppsigma+tpp_pi` -> `tppsigma-tpp_pi` on the second element),
and it is **backward compatible by construction**: both forms reduce to
`tppsigma*sx*sy` when `tpp_pi = 0`, and `TPPPI` defaults to 0, so every input written for
the upstream code is bit-for-bit unaffected. Confirmed: the NiO regression
(`example/NiO_octahedral_Ledge`) is unchanged -- `XAS_Ledge_{X,Y,Z}` max|diff| = 0.000e+00,
`eig.txt` identical, `dos.txt` identical, RIXS 1.4e-4 (BiCGStab tolerance). The octahedral
geometry has its own, more complete `tpp` construction (it keeps the cos*cos channels and
uses the `tpp_sig_pi` ratio rather than `TPPPI`) and is untouched.

After the fix our valence spectrum matches the legacy Fortran on all 14 lowest levels to
<= 0.0017 eV, including the level at 4.7134 that we previously placed at 4.7398, and the
one at 5.5571 that we previously placed at 5.5867.

Verified again against Quanty after the fix, now with the Fortran's own one-body matrix as
the common input: 80/80 of our bins matched, sum rule ratio 0.99999988, window-integrated
weights agreeing to 1e-7.

**Effect on results**: negligible for XAS on this parameter set (the L3 window weights
move by <0.001 and the main-line position not at all), but it changes the valence
eigenvalues by up to 0.03 eV and it is the ligand-ligand channel, so it should be
re-checked for anything driven by ligand dynamics -- which includes the pumped tr-RIXS.

## Result 3: the satellite discrepancy is retracted

`PRX_MAP_VALIDATION.md` sections 4 and 4b reported that the legacy Fortran puts **12.4%**
of the L3 XAS weight 4-9 eV above the main line (a "+1.90 eV satellite") where we put
1.3%, and treated that as the outstanding defect in our code.

**Both Fortran reference runs set `flag_skip_XAS = 1`** (`fortran_toy_reference` and
`fortran_toy_reference_v2`, line 31), i.e. neither computed an XAS spectrum, and neither
directory contains one. The only core-hole output is `cDOS.dat`, which is the
**unweighted** core-hole density of states -- a plain sum of Lorentzians over the core-hole
eigenvalues with no dipole matrix elements -- on **16** grid points spanning -1.4 to 24.5 eV
(1.73 eV spacing). The quoted reference numbers cannot be reproduced from it, and no other
file in the dataset is a candidate. So that comparison had no reference behind it and the
discrepancy is withdrawn.

For the record, in the corrected Hamiltonian the window 4-10 eV above the main L3 line
carries **exactly zero** weight, in our code *and in Quanty*. The L3/L2 branching ratio is
2.2288 in both, matching the atomic-limit value verified in `../atomic_limit_multiplet/`.

## Also resolved: the "open structural difference"

`PRX_MAP_VALIDATION.md` section 4b flagged, as the leading remaining candidate, that the
Fortran's interaction tensor `hint` is indexed by `(l, m_l)` with no site label, so
`rk(0,1,2,2,1) = cf0 = 6.0` would nominally apply to ligand-2p <-> metal-3d pairs -- an
inter-site monopole we do not have.

It does not. The Fortran gates every interaction matrix element through
`call checkint(u1,u2,u3,u4,...)`, whose own header reads *"make sure that we do not include
interaction terms that are unphysical such as core TM-ligand interaction"* / *"make sure
that only interactions on the same atom are considered"*. It zeroes any quadruple mixing
TM with ligand, mixing the two ligands, or (for `ncore = 6`, the TM L edge) mixing the TM
core with a ligand. The Fortran is strictly on-site, exactly like our `calc_coulomb`. The
caveat recorded in that section -- that a real 6 eV inter-site monopole would have destroyed
the observed <0.05 eV valence agreement -- was the correct instinct.

## Incidental: a latent off-by-one in the Fortran's `checkint`

The `ncore == 6` (TM L edge) blocks test `u > norbstm .and. u < norbspc` where every other
block in the routine uses `u <= norbspc`. With `norbstm = 5`, `norbspc = 11` and ligand
`uflag` values 6..11, that leaves `uflag = 11` -- the last ligand orbital -- **not** excluded
from interacting with the TM core hole, in the index ordering where the core index comes
first. It picks up the p-p integrals `upp0 = 1.0`, `upp2 = 0.1`. Small (~1 eV, one of six
ligand orbitals) and it cannot produce a several-eV satellite, but it is a genuine defect
and it is asymmetric in the index ordering, so it is worth knowing about.
