# 3d8 / Ni / TM L edge -- the PRX 15, 011012 parameter set, FG ordering CORRECTED

Translated from the paper-1 Fortran reference
(`Ni-L-RIXS/new_shape_function`, `TD_Cluster_Xray.f90`), which is the dataset behind PRX Fig. 7/8.

## The correction

`validation_modeB/modeb_pumpon_v2_fixed/INPUT` -- used for all of the Fig. 7/8 work so far -- has

    FG = 6.0  6.18  4.635  2.63        # WRONG for this code

The Fortran declares `cf2` = F^2 = 6.18, `cf3` = G^1 = 4.635, `cf4` = G^3 = 2.63, i.e. it orders
them `F0, F2, G1, G3`. **This code orders them `F0, G1, F2, G3`**, so F2 and G1 must be exchanged:

    FG = 6.0  4.635  6.18  2.63        # correct

Read the wrong way round, G1 and F2 are interchanged (6.18 <-> 4.635), which silently changes the
whole core-hole multiplet. Verified against Quanty 0.81 in the atomic limit -- with the ordering
above, this code reproduces Quanty's 2p^5 3d^9 multiplet exactly:

| energy above lowest | Quanty deg | this code deg |
|---|---|---|
| 0.0000 | 5 | 5 |
| 0.4045 | 21 | 21 |
| 1.4031 | 15 | 15 |
| 2.1883 | 9 | 9 |
| 4.1589 | 3 | 3 |
| 4.1769 | 7 | 7 |

(`validation_modeB/atomic_limit_multiplet/quanty_crosscheck/q_prx_d8.lua`.)

Effect on the equilibrium XAS against the Fortran's own `XAS_X.dat`, identically broadened:

| `FG` line | Pearson r | main L3 peak (ref -4.633) |
|---|---|---|
| `6.0 6.18 4.635 2.63` (as used so far) | 0.832 | -5.10 |
| `6.0 4.635 6.18 2.63` (correct) | **0.932** | **-4.600** |

## Remaining known discrepancy

Even with the ordering corrected, the Fortran puts **12.4%** of the total XAS weight 4-9 eV above
the main L3 line -- its satellite at +1.9 eV -- and this code puts **1.3%** there. The dipole sum
rule agrees to 0.1%, so the weight is misplaced rather than missing. This is unresolved, and it is
worth settling *before* trusting trRIXS maps built on the same core-hole eigenstates. Ruled out so
far: MLCT, tpd, sigpi, tpdzr, SO, uniform multipole scaling, the `FG` ordering itself, the exchange
sign (Quanty-confirmed correct), and missing dipole channels. See
`PRX_MAP_VALIDATION.md` section 4b.

## Other translation notes for this Fortran

| Fortran | here | note |
|---|---|---|
| `pump_ampx = 10.0` | `PUMPAMPX = 5.0` | their Peierls argument is `(k+A)/2`, ours `k/2+A` |
| `pump_sigma = 1.5` | `PUMPSIGMA = 1.2247` | theirs is `exp(-t^2/(2*sigma))`, i.e. sigma^2 = 1.5 |
| `pump_freq = 0.5` | `PUMPFREQ = 3.14159` | theirs is `sin(2*pi*f*t)`, ours takes the angular frequency |
| `gauss = 160*pi/(320*8)` | `TRPROBESIGMA = 5.0930` | `sigma_probe = 1/gauss` |
| `tch = 2*pi/1.5` | `TAUCH = 4.18879` | |
| `delxx = nd*aveu` | `MLCT = 10.900004` | re-derived from their f0/f2/f4; `EFFDEL = False` here |
