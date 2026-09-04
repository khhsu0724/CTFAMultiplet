# Next-Session Report: Two-Time-Correlator trRIXS/trXAS

**Read this first.** This is a concise orientation document written for a fresh session (possibly
using a different/higher-capability model) picking up this project. The full, chronological,
blow-by-blow debugging log from the session that produced this report is
[`SESSION_HANDOFF.md`](SESSION_HANDOFF.md) (600+ lines) — consult it for exact numbers, exact code
excerpts, and the full reasoning behind anything summarized below. Don't re-derive anything listed
under "Validated" or "Ruled out" without a specific new reason to doubt it.

## Reference papers (now saved in `reference_papers/`)

1. **`PhysRevX.15.011012.pdf`** — Jost, Lomeli, et al., *"Time-Resolved X-Ray Spectroscopy from the
   Atomic Orbital Ground State Up"*, Phys. Rev. X **15**, 011012 (2025). **This is the target
   paper** — the user (Eder Lomeli) is a co-author. It presents the trXAS/trRIXS results this
   project aims to reproduce, for both 3d⁹ (Cu-like) and 3d⁸ (Ni-like) cluster configurations.
   Fig. 2 (occupation dynamics), Figs. 3-6 (trXAS), Figs. 4-5,7-8 (trRIXS) are the key figures.
   **The specific end target is Fig. 7/8's "Hund's-exchange anti-Stokes" feature** — a peak at
   ≈−0.9 eV energy *gain* (anti-Stokes, i.e. the scattered photon has *more* energy than the
   incident one) appearing in the 3d⁸ cluster's trRIXS at long pump-probe delay, attributed to a
   pump-induced superposition of 3d⁸ configurations.
2. **`PhysRevB.99.104306.pdf`** — Chen, Wang, Jia, Moritz, Shvaika, Freericks, Devereaux,
   *"Theory for time-resolved resonant inelastic x-ray scattering"*, Phys. Rev. B **99**, 104306
   (2019). **This is the rigorous, general theory paper for tr-RIXS** — it derives the full
   four-time correlation function (their Eq. 9-11) from the Keldysh contour, then reduces it, with
   a phenomenological core-hole-lifetime decay `l(t1,t2) = exp(-|t2-t1|/τ_ch)`, to a two-time
   Kramers-Heisenberg cross section (their Eq. 18) that is structurally the *same* formula this
   project's "Mode B" (`trRIXS_2time`) implements and the same formula the paper-1 Fortran
   reference uses. **This paper was not available to cross-check against until this session** —
   see "Recommended next step" below for why this matters.

## The goal

Reproduce the trXAS/trRIXS results of paper 1 using this project's C++ exact-diagonalization
multiplet code, cross-validated directly against paper 1's own Fortran90 reference implementation
(`TD_Cluster_Xray.f90`, several variants under
`~/Documents/GradSchool/Devereaux Group/Herr_Jost/tr_spectroscopy/` on this machine). The specific,
concrete target is the 3d⁸/Ni/TM-L-edge anti-Stokes feature described above. This requires a
genuine **two-time correlator** propagator (paper 2's Eq. 18 / paper 1's Eq. 2) — a frozen-
instantaneous-Hamiltonian approximation is structurally incapable of it (see "Mode A vs Mode B"
below).

## What's validated (do not re-derive)

- **Hamiltonian and pump parameters** for the 3d⁸/TM-L-edge case, cross-checked byte-for-byte
  against 4+ variants of the Fortran source (`SO=11.5 0`, `CF=1.0`, `tpd=1.0`, `tpdzr=0.25`,
  `sigpi=0.225`, `TPPPI=0.25`, `FG=6.0 6.18 4.635 2.63`, `SC2=SC2EX=6.0 0 0.12857 0 0.025`,
  `MLCT=10.900004`, `PUMPFREQ=π`, `PUMPAMPX=5.0`, `PUMPSIGMA=1.224745` — this last value is
  `sqrt(Γ_p²)` with `Γ_p²=1.5 eV⁻²`, the pump width quoted in paper 1).
- **Occupation dynamics** (`scan_td_occupation` in `src/tdynamics.hpp`, dispatched via
  `TDOCCSCAN=True` in the INPUT `&PUMP` block) match paper 1's own Fig. 2(c1) closely — both
  against the Fortran reference numerically (~2-3% at spot checks) and, freshly this session,
  against paper 1's own plotted time window (Δt = −3 to 6 eV⁻¹) with a direct visual side-by-side
  (`validation_modeB/occupation_d8_prxfig2/occupation_vs_time_d8_comparison.png`) — same pre-pump
  plateaus, same oscillation onset, same long-time amplitude/period for all four tracked orbitals.
  **This confirms the core real-time propagator (`LanczosPropagate`/`build_trajectory` in
  `src/tdynamics.hpp`) is both correct and fast** (this code: ~2.6s for a 9 eV⁻¹ window; the
  full-scale Fortran: ~46s for a 75 eV⁻¹ window — roughly 2x faster per unit time-window even
  accounting for the different window sizes). **The propagator itself is not the problem** — see
  "What's failing" below; the issue is isolated to the two-time correlator's *assembly*, not the
  underlying single-state time evolution.
- **The pump/probe machinery must be genuinely two-time**, not a frozen-snapshot average — Mode A
  (below) cannot reproduce the anti-Stokes feature no matter how its parameters are tuned, because
  its internal oscillation period (~0.66 eV⁻¹) is faster than the core-hole lifetime
  (`τ_ch≈4.19 eV⁻¹`), and Mode A discards cross-time interference by squaring before averaging.

## Mode A vs Mode B (both exist in the code today)

- **Mode A** (`trXAS`/`trRIXS` in `src/tr_photon.hpp`, INPUT keys `TRXAS`/`TRRIXS`): treats each
  probe-window sample time as an independent frozen Hamiltonian, reuses the existing equilibrium
  solvers (`BiCGS`/`ContFracExpan`), and *incoherently averages the resulting intensities* across
  the probe window. Fast, validated, but structurally cannot capture the anti-Stokes feature.
- **Mode B** (`trRIXS_2time` in `src/tr_photon.hpp`, INPUT key `TRRIXS2`): builds the genuine
  two-time correlator `K(t1,t2) = <ψ_GS(t2)| D_out U_core(t2,t1) D_in |ψ_GS(t1)>` by real-time
  propagating the dipole-excited state, then Fourier-integrates over `t1,t2` with Gaussian probe
  windows and the core-hole-lifetime decay factor — this is paper 2's Eq. 18. **This is where all
  of this session's effort went, and where the current bug lives.** Not yet a clean production
  feature even once correct: `sim_window_mult`/`t2_max_window` are hardcoded literals in
  `main.cpp`'s dispatch (not real INPUT keys), there's no automatic `GSNEV`/`EXNEV` sufficiency
  check, and it writes one output file per (ground-state, matching-core-hole-block) group requiring
  external summation — all deferred, lower-priority cleanup once correctness is established.

## Validation methodology

Paper 1's own Fortran reference is `O(N_t²)` in the RIXS section and intractable locally at
production scale. This session built a *shrunk* toy copy (`nb=nt=1, timesteps=80`, saved at
`validation_modeB/fortran_toy_reference_v2/`) small enough to run in minutes, and compared its
output directly against the equivalent C++ Mode B run on the same toy-scale INPUT
(`validation_modeB/modeb_pumpon_v2_fixed/INPUT` and siblings). **Caution already learned the hard
way:** a shrunk window's `dt` must stay ≲0.05-0.08 or the Fortran's own per-step construction
doesn't converge — confirmed by watching a shrunk run itself converge as `timesteps` was increased.

## Bugs found and fixed this session (all real, all kept)

1. **Fortran reference crash**: a double-allocate of `allowed_tm` inside the `wsp` (final valence
   sector) loop crashes under `-fcheck=all` on the second match — silently truncated the
   *original* toy validation run to just 1 of 3 channels, undetected for a long time because the
   process still wrote a complete-looking output before dying. Fixed (in the `_v2` copy only) with
   a deallocate guard.
2. **Fortran reference file-naming collision**: the `tlabel` file-naming formula didn't depend on
   which `wsp` (core-hole block) was being written, so a second matching block silently overwrote
   the first's output files. Fixed with a `ppp`-dependent offset.
3. **Real C++ bug, kept**: `trRIXS_2time` computed each matching core-hole block (`exblk`) as an
   independent channel, incoherently summing their *intensities* — should coherently sum their
   *amplitudes* into one shared `Krows` before squaring once, since different core-hole blocks are
   different intermediate pathways to the *same* final state. Fixed in `src/tr_photon.hpp` (now
   takes `const std::vector<size_t>& ex_blk_inds` and accumulates coherently); confirmed this
   *doesn't* change the current toy validation's numbers (each channel only ever had one matching
   `exblk` there) but is still correct and necessary in general.

## Misdiagnoses ruled out this session (don't re-try without new evidence)

- **Missing per-delay probe-weight normalization** — falsified by direct Fortran source read: the
  reference's own two-time-correlator section never divides by any weight sum.
- **"3 channels" = 3 core-hole blocks for one ground state** — wrong; it's 3 *degenerate ground
  states* (the Sz=−1,0,+1 sublevels of a spin-triplet 3d⁸ ground state) for this toy system.
  Correctly identified which one (`chan0`, 2Sz=−2) matches the Fortran's own `ws=1`.
- **Lanczos/Krylov propagator convergence** — `TDM=50` vs `TDM=20` gave bit-for-bit identical Mode
  B output. Ruled out definitively; not a numerical accuracy issue.
- **Degenerate final-state coherent summing** — a real, general phenomenon (coherently summing a
  degenerate eigenspace is basis-dependent, since different diagonalization routines pick different
  arbitrary relative phases within it) but tested two ways (fully per-eigenvector incoherent, and
  grouped-by-degenerate-energy incoherent) and **both made an already-good feature measurably worse**
  without fixing the bad one. Reverted. See `SESSION_HANDOFF.md`'s "TESTED AND REJECTED" section for
  the full numbers and reasoning — this is a legitimate theoretical concern in general but is not
  (or not the dominant part of) this specific bug.

## What's failing right now

For the correctly-matched single channel, comparing the toy C++ Mode B output against the Fortran
reference across 5 delays at the 3d⁸ dd-excitation feature (`ab=-4.992 eV, loss≈2 eV`): the ratio
ranges from **0.7× to 575×** depending on delay — not a constant scale factor, not monotonic with
delay. A *different* feature in the same spectrum (`ab≈12.24 eV`) matches much better (~1.2×), so
the bug seems to specifically afflict the sharp, long-lived dd-excitation resonance rather than the
whole spectrum uniformly.

The most targeted diagnostic run this session did: dump the **raw** `K(t1,t2)` amplitude (built
once, shared across every incident energy — only the final `t1,t2` Fourier sum depends on incident
energy) and compare it directly, point-by-point, against the reference's own `tkern(t1,t2)`, at
both short (`t2-t1` up to ~0.4 eV⁻¹) and long (`t2-t1` up to ~6.3 eV⁻¹) propagation times. Result:
a real, persistent **magnitude mismatch of ~2.8-4.8×** (not 1.0×, not growing with propagation
time — ruling out an accumulating/compounding numerical error) with an **irregular phase mismatch**
(not constant, not a clean conjugate, not a linear drift). This proves a real bug exists in the
raw amplitude construction, upstream of the frequency-domain step — but a ~3-5× point-wise error is
too small on its own to obviously explain a 575× final blowup, so either (a) there's a *second*,
separate bug in the final `t1,t2`/frequency summation loop, or (b) the coherent double-sum needed to
build a sharp resonance genuinely amplifies a modest phase error this much in practice (physically
plausible, since sharp resonances require long-range phase coherence) — this wasn't distinguished
before the session ended.

## Recommended next step (new since the last handoff update)

**Cross-check the two-time-correlator formula itself against paper 2's rigorous derivation**, not
just against paper 1's own Fortran implementation. Everything done this session compared this
project's C++ against *one specific numerical implementation* (the Fortran reference) — a
structural/numerical match to that code doesn't rule out both codes sharing a subtly wrong
convention, and a mismatch doesn't distinguish "our code is wrong" from "the reference's own toy
shrinking broke something" without a third, independent, from-first-principles check. Paper 2's
Eq. 9-18 derives the tr-RIXS cross section generally (not tied to any specific code), including the
exact role of the core-hole lifetime `l(t1,t2)`, the exact time-ordering of the four dipole
operators, and the reduction from the full four-time correlator down to the effective two-time form
via the probe-shape-function approximation `g(t1,t)g(t2,t)g(t1',t)g(t2',t) ≈ g²(t1,t)g²(t2,t)`
(their Eq. 18's derivation, valid when core-hole lifetime ≪ probe duration — worth explicitly
checking this condition holds for the parameters in use). Read paper 2 §II.B-C in full (not just the
skim done this session), derive the *exact* expected `K(t1,t2)`/`csum` formula term-by-term from
Eq. 9-18, and compare that derivation — not just the existing code — against what
`trRIXS_2time` actually computes line-by-line. This is the one class of check that hasn't been
tried yet and could reveal a genuine, structural formula error (a missing factor, a wrong
time-ordering, a missing complex-conjugation) rather than another numerical-implementation-level
discrepancy.

Other concrete, not-yet-exhausted candidates, roughly in order of effort-to-payoff:

1. Line-by-line (not structural) comparison of the dipole matrix element construction
   (`basis_overlap`/`gen_dipole_state` in `src/photon.cpp`) against the Fortran's own vertex
   construction (~line 2400-2430 of `TD_Cluster_Xray.f90`) — checked structurally this session
   (plausibly matches), never at full numerical/line-by-line precision.
2. Whether `trRIXS_2time` needs an `E_gs`-referencing term analogous to Mode A's
   `instantaneous_energy` (`src/tr_photon.hpp`) — raised but not conclusively ruled in or out; peak
   *locations* already matching between codes is evidence against a gross missing offset, but a
   subtler version of this idea (e.g. affecting only the propagation phase, not the final
   `win`/`wout` axis) hasn't been fully checked.
3. Whether the toy validation window's small size (`nb=nt=1`) itself produces an artifact common to
   both codes that a larger window wouldn't — the original motivation for this session's "v2"
   diagnostic test, revisited multiple times but never cleanly isolated because other confounds
   (the Fortran crash, the degenerate-ground-state mix-up) kept surfacing first. Worth returning to
   *only after* the raw-amplitude bug above is resolved, since it's a much more expensive test to
   run (needs `nb,nt` big enough to matter, `O(N_t²)` scaling) and shouldn't be chased while a
   cheaper, more targeted lead is still open.

## Key file locations

- `NEXT_SESSION_REPORT.md` (this file) — start here.
- `SESSION_HANDOFF.md` — full chronological debugging log, 600+ lines, exact numbers and reasoning
  for everything summarized above.
- `reference_papers/PhysRevX.15.011012.pdf` — target paper (paper 1).
- `reference_papers/PhysRevB.99.104306.pdf` — tr-RIXS theory paper (paper 2), **not yet used** for
  a from-first-principles formula cross-check — see "Recommended next step."
- `src/tr_photon.hpp` — `trRIXS_2time` (Mode B) and `write_trRIXS_2time`; also Mode A's
  `trXAS`/`trRIXS`/`instantaneous_energy`.
- `src/main.cpp` — `TRRIXS2` INPUT parsing (~line 280s) and dispatch (~line 722-760).
- `src/tdynamics.hpp` — `LanczosPropagate`, `build_trajectory`, `scan_td_occupation` (all
  validated/trusted).
- `src/photon.cpp` — `basis_overlap`, `gen_dipole_state`, `proj_pvec` (dipole matrix element
  construction, shared by Mode A/B and the older equilibrium solvers).
- `validation_modeB/fortran_toy_reference_v2/` — the fixed toy Fortran reference (crash and
  file-naming bugs patched), its compiled binary, and its output.
- `validation_modeB/modeb_pumpon_v2_fixed/` — the matching C++ Mode B toy run (coherent exblk
  summation, the current best-performing version).
- `validation_modeB/occupation_d8_prxfig2/` — the fresh occupation-vs-time validation against
  paper 1's own Fig. 2(c1), including the side-by-side comparison PNG.

## Commit state

Everything described above is committed to `main` as of this report. Run `git log --oneline -20`
for the exact commit sequence and messages (each commit message explains *why*, not just *what*,
for the corresponding step above).
