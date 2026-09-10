# NiO, octahedral, TM L edge -- validated equilibrium parameter set

Eder's NiO equilibrium calculation, kept here because it is a **verified regression target**: with
no environment variables set, this code reproduces its stored output with `eig.txt`, `dos.txt`,
`tmat.txt`, `exeig.txt`, `exdos.txt` and all three `XAS_Ledge_*.txt` files **bit-identical**, and
the six cross-polarisation `RIXS_Ledge_*.txt` files to ~1e-4 (BiCGStab's own convergence
tolerance -- `CG_tol = 1e-8` is tested against the *squared* residual norm).

Stored reference output and the comparison scripts live in `validation_modeB/nio_equilibrium/`.
Source of the reference:
`~/Documents/GradSchool/Devereaux Group/Herr_Jost/NiO_ED_tr_magnons/equlibrium_NiO/`.
Plotting recipe: `validation_modeB/nio_equilibrium/plot_nio.py` (follows Eder's
`plot_RIXS_XAS.ipynb`: summed cross-polarisations, `shift = 862.9`, XAS overlaid in white).

## Parameters

| key | value | note |
|---|---|---|
| `SO` | `11.507 0.083 0.102` | 2p core SOC, then valence SOC terms -- valence SOC is ON here |
| `SC1` | `1.0 0 0.1` | ligand p F0, F2 |
| `SC2` | `6.0 0 12.233 0 7.597` | valence d-d F0, F2, F4 (unreduced; >1 so the code's /49, /441 conversion fires) |
| `SC2EX` | `6.0 0 13.005 0 8.084` | same in the core-hole state |
| `FG` | `6.0 5.783 7.720 3.290` | **F0, G1, F2, G3** -- see the convention note below |
| `HFscale` | `0.8` | scales `SC2[1..]`, `SC2EX[1..]`, `FG[1..3]`, not the F0 monopoles |
| `CF` | `0 0 0.8 0.8 0.8` | crystal field on dxy, dxz, dyz |
| `tpd`, `sigpi`, `tpp` | `0.35`, `0.25`, `0.0` | hybridization |
| `MLCT` | `4.5` | with `EFFDEL = True`, so this is the *effective* delta |
| `MF` | `0.11` | mean field |
| `Coordination` | `Octahedral`, 1 site, 2 holes | Ni2+ d8 |
| `DIAG` | `4` | blocks here are < 2e3 so `Dense` is chosen anyway |
| `epsab`, `epsloss` | `0.4`, `0.01` | incident and loss broadening |
| `INCIDENT` | `-7.9 2.3 100` | start, end, number of points |

## FG convention -- the one real gotcha

`FG` is **`F0, G1, F2, G3`**, per upstream CTFAMultiplet
(https://github.com/khhsu0724/CTFAMultiplet/wiki/Input-Parameters). So this input carries
F2 = 7.720 and G1 = 5.783, which after `HFscale = 0.8` gives F2 = 6.176, G1 = 4.626, G3 = 2.632 --
matching Hartree-Fock for Ni2+ 2p3d (7.721, 5.787, 3.291) times 0.8. **This input is correct as
written.**

**When translating from the Fortran (`TD_Cluster_Xray.f90`), F2 and G1 must be exchanged.** The
Fortran declares `cf2` = F^2, `cf3` = G^1, `cf4` = G^3, i.e. `F0, F2, G1, G3`:

    Fortran: cf0 = 6.0   cf2 = 6.18 (F2)   cf3 = 4.635 (G1)   cf4 = 2.63 (G3)
    ours   : FG = 6.0   4.635   6.18   2.63

Getting this backwards is easy and quiet -- it changes the core-hole multiplet without any error
message. `validation_modeB/modeb_pumpon_v2_fixed/INPUT` has it backwards; see
`validation_modeB/atomic_limit_multiplet/` for how to check a parameter set against Quanty in the
atomic limit.

## Verifying a core-valence parameter set

`validation_modeB/atomic_limit_multiplet/quanty_crosscheck/q_pd.lua` builds the same 2p^5 3d^9
multiplet in Quanty. Switch off hybridization, CF, SOC and MF in the INPUT, run both, and compare
`exdos.txt` against Quanty's eigenvalues -- they agree here to 4.7e-05 eV with the degeneracy
sequence 5, 21, 15, 9, 3, 7. `validation_modeB/pdterms.py` computes the same multiplet analytically
from Wigner 3j/6j symbols as a third check.
