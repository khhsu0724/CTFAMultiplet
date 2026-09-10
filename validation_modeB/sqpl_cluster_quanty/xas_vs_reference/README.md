# The +1.9 eV XAS satellite = the legacy core-hole ligand-ligand hopping sign

`reference_XAS_X.dat` is the production Fortran run's own XAS (columns: incident energy,
dipole-weighted XAS, unweighted core DOS, incident energy again; Lorentzian `gammax = 0.1`),
copied from

    Ni-L-RIXS/new_shape_function/XX/16_points_new_shape_function_XX/XAS_X.dat

Both of our spectra are pole lists (`solver = 1`, `NEDOS = 40000`) broadened with the same
Lorentzian on the same grid. Windows are relative to each spectrum's own main L3 peak.

| | main L3 | w[-3..+1] | w[+1..+4] | **w[+4..+9]** | w[+14..+22] | r vs ref |
|---|---|---|---|---|---|---|
| reference | -4.633 | 0.4212 | 0.1294 | **0.1238** | 0.2505 | 1.0000 |
| ours, default | -4.633 | 0.5062 | 0.1614 | **0.0067** | 0.3091 | 0.8691 |
| ours, `EX_HOP_SIGN=-1` | -4.633 | 0.4212 | 0.1294 | **0.1238** | 0.2505 | **1.0000** |

`EX_HOP_SIGN=-1` gives the ligand-ligand hopping the opposite sign in the core-hole sector to
the valence sector, which is what the reference does (its line ~1555 `+4*tmato` against ~764
`-4*tmato`). It is not a gauge choice -- the p -> -p gauge on the ligand orbitals leaves that
term invariant, so the two sectors are genuinely different Hamiltonians, differing by up to
0.70 eV on the one-body levels.

An O-O hopping integral cannot depend on whether a core hole sits on the metal, so our default
is the physical convention and the reference's +1.9 eV satellite appears to be an artifact of
that sign rather than a charge-transfer feature. Stated as an argument, not a measurement: no
experimental Ni L3 XAS is compared here, and the published Fig. 7 maps inherit the same
convention.
