# RIXS was dropping the spin-flip channel (43-54% of the L-edge cross section)

Found by cross-checking equilibrium RIXS against Quanty 0.81 on the square-planar d8
cluster, after the XAS had already been verified exact (`../sqpl_cluster_quanty/`).

## The bug

All three RIXS paths de-excited the intermediate state back into **the same ground-sector
symmetry block the initial state came from**:

| path | file | what it did |
|---|---|---|
| equilibrium `RIXS()` | `src/photon.cpp` | `ContFracExpan(GS.hblks[g.first].ham, ...)` |
| Mode A `trRIXS_snapshot` | `src/tr_photon.hpp` | `ContFracExpan(GS.hblks[g.first].ham, ...)` |
| Mode B `trRIXS_2time` | `src/tr_photon.hpp` | final-state basis `{|l>}` built only from `GS.hblks[gs_blk_ind]` |

That is wrong whenever the ground sector is split into Sz blocks. The core hole carries
spin, and **core spin-orbit coupling mixes it while the intermediate state propagates**, so
absorbing an up valence hole and emitting a down one changes Sz by 1. Those final states
live in a different block and were silently never reached. That is the spin-flip (magnon)
RIXS channel -- the physics this repository is named for.

## Evidence

Quanty is set up in the hole language with our own one-body matrix (same construction as
`../sqpl_cluster_quanty/`), and the Kramers-Heisenberg sum is assembled from exact matrix
elements over the **complete** 231-state valence space:

    I(w_in, w_loss) = sum_i sum_f | sum_n <f|D_out^dag|n><n|D_in|i> / (E_n - E_i - w_in + i*eps_ab) |^2
                      * L(w_loss - (E_f - E_i), eps_loss)

with `eps_ab = eps_loss = 0.1` (the code's own defaults, confirmed from its log).

| channel | | scale ours/Quanty | Pearson r | max&#124;diff&#124;/peak |
|---|---|---|---|---|
| X in, Y out | before | **0.567** | 0.9306 | 28.3 % |
| X in, Y out | after  | **0.999999** | 0.9999972 | 0.37 % |
| X in, X out | before | **0.457** | 0.7737 | 56.9 % |
| X in, X out | after  | **0.999971** | 0.9999998 | 0.10 % |

The decisive number is the **scale**, which is not a fitted parameter -- it is the ratio of
absolute cross sections. Before the fix we had 57% (XY) and 46% (XX) of the true intensity.
The residual 0.1-0.37% after the fix is BiCGStab's tolerance plus the continued-fraction
truncation.

Individual features were misweighted by up to 10x, not merely scaled: on resonance the
loss = 3.66 eV final state came out at 0.11x its correct weight while the elastic line in
the XX channel came out 1.75x too strong.

Mode B was worse still. On the `example/PRX_d8_sqpl_Ledge` toy case it reported

    before:  Building 55 valence-eigenstate trajectories (de-excitation basis)...
    after :  Building 231 valence-eigenstate trajectories (de-excitation basis, 3 ground block(s))...

i.e. it summed over 55 of 231 final states -- missing 76% of the final-state space.

## Scope: which inputs were affected

Only those whose **ground sector is Sz-blocked**, which means valence spin-orbit coupling
off *and* no mean field. Concretely:

- **Affected**: `example/PRX_d8_sqpl_Ledge` and everything derived from it (`SO = 11.5 0`,
  no `MF`) -- 3 ground blocks of dimension 55/121/55, 1 core-hole block of 132. This is the
  PRX Fig. 7/8 parameter set, so **every tr-RIXS result in this repository built on it is
  affected** and has been re-run.
- **Not affected**: `example/NiO_octahedral_Ledge` (`SO = 11.507 0.083 0.102`, `MF = 0.11`)
  -- valence SOC and a mean field leave the ground sector as a **single** 168-dimensional
  block, so every final state was already included. Verified: the stored NiO reference is
  reproduced identically before and after (`RIXS_Ledge_X_Y` max|diff| 1.4e-4 = BiCGStab
  tolerance in both cases, totals 172.896567 in both).
- **Not affected**: anything with no core spin-orbit coupling. With `SO = 0 0` the fix is
  provably inert -- `basis_overlap` then finds no dipole overlap between a core-hole block
  and any other ground block. Verified bit-identical: **max|diff| = 0.000e+00**.

So no input written for the upstream code changes behaviour unless it was actually hitting
the bug.

## Cost

The fix adds one dipole projection and one continued fraction per extra ground block per
incident energy. For Mode B it also multiplies the number of valence trajectories and the
size of the `Krows` array by (total valence dim)/(one block's dim) -- 231/55 = 4.2x on the
toy case, 231/121 = 1.9x when starting from the Sz = 0 block. That is unavoidable: those
amplitudes are the missing physics.
