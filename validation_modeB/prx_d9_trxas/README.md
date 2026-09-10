# 3d9 (CuO) cluster, TM L3 edge: our code against PRX Fig. 2(b) and Fig. 3(d)-(f)

Set up from `Herr_Jost/tr_spectroscopy/Cu/Cu-K_X/TD_Cluster_Xray.f90`, which is the 3d9 cluster
configured for the **ligand K edge** (`edge = 2`, `ncore = 4`) -- the paper's Fig. 3(a)-(c).
Switched here to the **TM L edge**, for which that dataset contains no Fortran output, so the
comparison is against the published figures. Parameter translation is documented in
`example/PRX_d9_sqpl_Ledge/INPUT`.

Equilibrium check: ground state `d9 : 0.77881, d10L : 0.22119`, `Delta = 2.8`, L3 line at
+1.50 eV, L2 at +18.75 (core spin-orbit splitting 17.25 = 1.5 * 11.5).

## What matches the paper

| quantity | paper | ours |
|---|---|---|
| 3d_x2-y2 occupation, dominant long-time frequency | "about 6.8 eV" [Fig. 2(b2)] | **6.77 eV** |
| 2p_aa occupation | same frequency **plus two additional** components, out of phase [Fig. 2(b1)] | 6.77 eV **plus 2.66 and 2.90 eV** |
| L3 peak softening at `Delta t ~ 0` [Fig. 3(d)] | yes | yes |
| net L3 shift at long times | ~ -1.0 eV (read off Fig. 3(d), +-0.3) | **-1.16 eV** |
| transient dip depth | ~ -2.1 eV (read off, +-0.3) | -1.34 eV at `Delta t = +1.75` |
| main-peak intensity: sharp drop then oscillations [Fig. 3(e)] | yes | yes, same shape |
| X vs Y dichroism out of phase [Fig. 3(f)] | yes | yes, **correlation -0.994** |

The 6.77 eV hit is on a specific number the paper quotes, and the ligand orbitals reproducing that
frequency *plus two extra components* is the paper's own description of Fig. 2(b1) verbatim. So the
3d9 pump dynamics is independently confirmed.

Read-off values from the published figure carry roughly +-0.3 eV; the axis spans 930-933 eV in
about 2 cm.

## The mechanism, and why this cannot settle the 3d8 sign question

For d9 there is **one** hole. XAS conserves hole number in the hole picture -- it creates a 2p hole
and destroys a 3d hole -- so the core-excited sector is `2p^5 3d^10` with **no valence holes at
all**: exactly 2 levels, `-5.75` (x4, L3) and `+11.5` (x2, L2), pure core spin-orbit coupling.
`exdos.txt` confirms it.

Two consequences.

1. **The L3 line position is exactly `E_EX - E_GS(t)`**, and `E_EX` is pump-independent (no valence
   holes means no core-valence Coulomb and no hopping to modulate). So the permanent red shift must
   equal minus the permanently deposited energy. Verified: the pump deposits **+1.1847 eV** into
   the 3d9 ground state (`E_stat`: -1.5073 -> -0.3227) and the L3 line shifts by **-1.16 eV**. The
   agreement is exact to the bin width. This makes the d9 L edge an unusually clean probe -- there
   is no multiplet or hybridization structure in the final state to complicate it.

2. **The d9 L edge is exactly insensitive to the core-hole ligand-ligand hopping sign.** With no
   valence holes in the excited sector there is nothing for that hopping to act on. Verified:
   `EX_HOP_SIGN=-1` changes the XAS by `max|diff| = 0.000e+00` and leaves `eig.txt` and `exdos.txt`
   byte-identical. **So this calculation cannot adjudicate the 3d8 sign question**, and it was
   worth checking before running it.

What it *does* establish is the counter-consideration raised against the artifact reading of the
3d8 shift: our code is not unable to produce pump-induced L3 red shifts. Where the mechanism is
unambiguous -- pure energy deposition, no final-state complications -- it produces one of very
nearly the published size (-1.16 vs ~-1.0 eV). In the 3d8 trRIXS the core-excited state *does*
carry a valence hole, so its energy shifts under the pump too and the two shifts largely cancel,
which is why we get -0.20 eV there. That is consistent with the reference's extra ~1.8 eV coming
from the sign, though it does not prove it.

## Incidental: a paper-vs-code parameter mismatch

The paper's appendix [45] quotes `eps_p = 2.8` for the 3d9 cluster, which matches this Fortran
exactly. For the **3d8** cluster it quotes `eps_p = 6`, but that Fortran instead computes
`delxx = nd*aveu = 2 * 5.450003 = 10.900006` and uses it. Anyone reproducing the 3d8 results from
the paper's stated parameters will get a charge-transfer energy about 1.8x too small. Our
`example/PRX_d8_sqpl_Ledge/INPUT` uses `MLCT = 10.900004`, i.e. it follows the code, not the paper.

## Reproducing

    # equilibrium
    INPUT_equilibrium        -> ground state, L3/L2 positions, sign-insensitivity check
    # trXAS, dense delays and a narrow probe so the transient and the oscillations resolve
    INPUT_X, INPUT_Y         -> TRPROBESIGMA = 0.5, TRDELAY = -6 20 105
    python3 plot_d9_trxas.py <dir_X> <dir_Y>
    # occupation dynamics; TDDT = 0.005 because print_every = 20 is hardcoded, so the
    # output step is 20*TDDT and 6.8 eV needs adequate sampling
    INPUT_occupation_scan
