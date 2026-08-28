# Session Handoff: Time-Resolved RIXS (trRIXS) Two-Time Correlator Implementation

**Purpose of this document**: this session is approaching its context limit. This file is a
complete technical handoff so a fresh Claude session (or the user) can pick up exactly where
this one left off, without needing to re-derive anything already established here.

## Overall goal

Reproduce the time-resolved RIXS/XAS results of **Phys. Rev. X 15, 011012 (2025)** ("Time-resolved
X-ray spectroscopy from the atomic limit") using this C++ exact-diagonalization code
(CTFAMultiplet-magnons), cross-validated against the paper's own Fortran90 reference
implementation (`TD_Cluster_Xray.f90`, several variants supplied by the user from
`~/Documents/GradSchool/Devereaux Group/Herr_Jost/tr_spectroscopy/`). The specific target: the
paper's Fig. 7/8 "Hund's-exchange anti-Stokes" feature at ≈0.9 eV energy gain in the d8 (Ni,
TM L-edge) cluster, appearing at long pump-probe delay, off the main absorption resonance.

## What is already solidly validated (do not re-derive)

1. **Hamiltonian + pump parameters** for the d8/TM-L-edge case (`CF=1.0`, `SO=11.5 0` i.e.
   dlso=0, `tpd=1.0` from `sqrt(3)*vpds/2`, `tpdzr=0.25`, `sigpi=0.225`, `TPPPI=0.25`,
   `FG=6.0 6.18 4.635 2.63`, `SC2=SC2EX=6.0 0 0.12857 0 0.025`, `MLCT=10.900004`,
   `PUMPFREQ=pi`, `PUMPAMPX=5.0`, `PUMPSIGMA=1.224745`) were derived directly from the Fortran
   source and cross-checked byte-for-byte across 4+ file variants (`NiO_final_final_Ledge`,
   `check_AS`, `new_shape_function`, and 3 files the user downloaded — files (1) and (3) are
   IDENTICAL in physics to `new_shape_function/XY`; file (2) is a different, untested
   d9/TM-L-edge configuration, not yet explored).
2. **Occupation dynamics match the Fortran exactly** (not just qualitatively vs. the paper's
   Fig. 2) — confirmed to ~2-3% at t=0,6,10 eV⁻¹ using the *unmodified, full-scale* reference
   run (`nb=4, nt=8, timesteps=320`). **Important pitfall already hit and resolved**: an
   earlier *shrunk* toy comparison (nb=1) appeared to show a huge, alarming mismatch — this
   turned out to be a **coarse-`timesteps`/`dt` convergence artifact**, not a physics bug (the
   Fortran's own per-step construction converges slowly with `dt`; confirmed by rerunning the
   Fortran itself at `timesteps=20→200→2000` and watching it converge toward my value). Lesson:
   never trust a shrunk-window toy comparison without first checking `dt`-convergence.
3. **The pump/probe machinery is genuinely two-time**, not a simple frozen-snapshot average.
   The existing "Mode A" implementation (`trXAS`/`trRIXS` in `src/tr_photon.hpp`, dispatched via
   `TRXAS`/`TRRIXS` INPUT keys) treats each probe-window sample time as an independent frozen
   Hamiltonian and averages the resulting *intensities* — this is a real approximation, and
   it's why Mode A never reproduced the 0.9 eV anti-Stokes feature no matter how the Hamiltonian
   parameters were tuned (confirmed via many `check_AS`/`new_shape_function` reruns, none of
   which are worth re-doing — see the "what NOT to re-try" section below).

## The Fortran reference: how to compile and run a cheap toy test

```bash
mpif90 -O3 -o code -ffree-line-length-none -fcheck=all -L/opt/homebrew/lib -larpack \
  -F/System/Library/Frameworks/Accelerate.framework -framework accelerate TD_Cluster_Xray.f90
```
(mpif90, ARPACK via homebrew, and Accelerate are all present on this machine — confirmed working.)

**The full-scale RIXS section (`tpivert` construction) is O(N_t²) and NOT tractable locally**
(N_t = (nb+nt)*timesteps ~ 3840 for the paper's real parameters). For a cheap local toy test,
shrink `timesteps`, `nb`, `nt` in a *copy* of the source, but:
- **Keep `dt = 2*pi/timesteps` no larger than ~0.05-0.08** (finer than this and the Fortran's
  own construction doesn't converge — confirmed empirically, see above). `timesteps=80` with
  `nb=1, nt=1` (giving dt≈0.0785, total window ±2π≈±6.28) is the configuration validated this
  session — runs in a few minutes total including the RIXS section.
- Also shrink `nfreq`/`ifreq` (frequency grid sizes, e.g. to 16 each) and replace the
  `allowed_tm` delay list with values that fit inside `[-nb*timesteps, nt*timesteps]` (e.g.
  `[-72,-64,...,72]` for `timesteps=80`).
- **Kill the process as soon as `Occ_Pulse.dat` is fully written** if you only need occupation
  data (it's written well before the expensive RIXS section starts) — use
  `until [ -s Occ_Pulse.dat ] || ! ps -p $PID; do sleep 3; done` then `kill -9 $PID`.

**A working toy reference set (source, compiled binary, and all output files) is saved at**
`validation_modeB/fortran_toy_reference/` in this repo (survives `/tmp` getting cleared, which
happened once already this session on a reboot — never rely on `/tmp` surviving across
sessions). The trSPEC files there are for **XY polarization only** (this specific Fortran
variant's `beta` loop is restricted to `beta=2,3`, i.e. Y,Z scattered only, not X) — pvin=X,
pvout=Y in any C++ comparison.

## The Mode B implementation (the new work this session)

New function `trRIXS_2time()` in `src/tr_photon.hpp` (near the end of the file), dispatched
from `main.cpp` via a new `TRRIXS2` boolean INPUT key (in `TRParam::trRIXS2`). This is
currently a **validation-only path**: it hardcodes some parameters in `main.cpp`'s dispatch
block (a `sim_window_mult` value chosen to match the toy Fortran's exact window, and
`t2_max_window=-1` for "no truncation") rather than exposing them as general INPUT keys — that
generalization is future work, not yet done.

### What it computes

```
csum(win, wloss, T) = sum_{t1<t2, t2-t1<=window} g(t1-T) g(t2-T) exp(-(t2-t1)/tau_ch)
                        K(t1,t2) exp(-i*win*t1 + i*(win-wloss)*t2) * dt^2
I(win, wloss, T) = |csum|^2
```
where `K(t1,t2) = sum_l <valence_l(t2)| D_out  U_core(t2,t1)  D_in |psi_GS(t1)>`, `l` ranging
over the **complete valence eigenstate basis** of the relevant block (not just the ground
state) — see "the bug that was found and fixed" below. `g(x)` is a normalized Gaussian probe
window matching the Fortran's own `gauss*exp(...)/sqrt(2*pi)` normalization exactly.

### Efficiency levers relative to the reference (all still valid, all implemented)

1. `K(t1,t2)` does not depend on delay `T` — built once, reused for every requested delay
   (the reference recomputes it once per delay, redundantly).
2. `t2` is truncated at a `tau_ch`-bounded window when desired (`t2_max_window` parameter) —
   the reference sums to the end of the simulated range regardless.
3. `U_core(t2,t1)` is applied via short-iterative-Lanczos propagation of a single state vector
   (reusing the existing `propagate_block_state`, unmodified, applied to the `EX` Hilbert space
   instead of `GS` — it was already fully generic), not full dense diagonalization plus an
   `ndim x ndim` matrix composition at every timestep like the reference does.

### Two bugs found and fixed this session, in order

1. **Equal-time (`t1=t2`) term.** Initially the K(t1,t2) loop started at `k=0` (i.e. included
   `t2=t1`); the reference's own loop (`do t2=t1+1,...`) never includes this term. Fixed by
   starting the loop at `k=1`. **This did NOT turn out to be the actual bug** — a before/after
   pump-off test gave bit-for-bit identical results, proving this term was already negligible.
   Kept the fix anyway (it's still correct to match the reference exactly), but kept looking.
2. **THE actual bug, confirmed and fixed**: the de-excitation ("l") side of `K(t1,t2)` was
   using only the single physical ground-state trajectory (`gs_traj`) instead of the reference's
   true structure — a **coherent sum over the complete valence eigenstate basis** (`l` ranges
   over `min(nev, binstates(wsp))` in the Fortran, which for this toy system is `nval=55`
   states, not "a few" — confirmed by direct inspection of the code). Fixed by building `nval`
   separate trajectories (each starting from a different static valence eigenvector, propagated
   forward under `GS`'s own real-time Hamiltonian) and — since inner products are linear in the
   bra — **pre-summing them into a single `combined_traj` per timestep** rather than doing `nval`
   separate inner products inside the expensive K(t1,t2) double loop (a real efficiency win: the
   summation cost is paid once per timestep during trajectory construction, not once per (t1,t2)
   pair). This fix was confirmed decisively via a pump-off equilibrium sanity check: before the
   fix, the spectrum showed a spurious dominant peak at loss=0 (elastic-line-like, ~1700x
   larger than the genuine loss=2 dd-excitation feature); after the fix, loss=2 became the
   correctly dominant peak, matching the reference's own equilibrium shape.

**Practical note for continuing this work**: this fix requires `GSNEV`/`EXNEV` set large enough
to cover the *entire* block dimension (not just a handful of low-lying states) — the toy system
needs `GSNEV=EXNEV=150` (actual block size is 55) to get the full basis; for the paper-scale
system this will need to be re-checked (larger blocks may need larger NEV, or a full
diagonalization).

## Current state: where the comparison stands right now

With both fixes applied, compared against `validation_modeB/fortran_toy_reference/` across 7
interior delays (see `validation_modeB/compare_modeb_full.py`, results saved in
`validation_modeB/modeb_pumpon_result/` and `modeb_pumpoff_result/`):

- **Peak location matches exactly** at 6 of 7 delays: both codes agree the dominant feature
  sits at `ab=-4.992, loss=2.00` — strong confirmation the core resonance structure and energy
  bookkeeping are now correct.
- **Peak magnitude is systematically too large, and — as of the latest re-analysis below — the
  *shape* of the magnitude vs. delay is wrong, not just a scale factor.** (An earlier pass at this
  session described this as a simple "ratio that grows with delay" and hypothesized a missing
  normalization; that hypothesis has since been **falsified** by direct source inspection — see
  "Immediate next step" below for the corrected picture and current best explanation.) Both `ref`
  and `mine` show a rise-then-fall shape in peak intensity vs. delay, but `ref` peaks near
  `T≈-1.88` while `mine` peaks later (`T≈0` to `+1.88`), and `ref` falls off much faster past its
  peak than `mine` does.

- **The last delay (+5.65) breaks down entirely** — the reference's peak has moved to a
  completely different (ab, loss) location there. This delay sits only ~0.63 eV⁻¹ from the
  simulation's hard edge (`t_hi=+6.28`), while `probe_sigma=5.09` is much wider than that
  margin — so the Gaussian probe window is barely attenuated by the time it hits the edge. Given
  the newer finding that even the *other 6* delays show a shape mismatch (not just this one), this
  is now understood as the most extreme case of a truncation effect that likely touches the whole
  delay range in this particular toy window, not an isolated one-off artifact.

## Mode A vs Mode B selection: current state, and the user's plan

Mechanically, the INPUT file's `&PUMP` block already has three booleans that select between
modes: `TRXAS`, `TRRIXS` (Mode A), `TRRIXS2` (Mode B) — parsed at `main.cpp` lines ~284-286,
dispatched in that priority order (`TRRIXS2` checked first) starting at line ~722. So switching
is *possible* today by editing the INPUT file. But it is **not yet a clean, general-purpose
toggle** — Mode B is currently a validation harness, not a finished feature:

1. `sim_window_mult` and `t2_max_window` are **hardcoded literals in `main.cpp`'s dispatch
   block** (not read from any INPUT key), tuned specifically to match the toy Fortran
   reference's exact bounded time window. They are almost certainly wrong for a different
   (e.g. full physical/paper-scale) run and require editing the source directly to change.
2. `GSNEV`/`EXNEV` must be set (via INPUT) large enough to cover the *entire* relevant block
   dimension for Mode B's de-excitation sum to be correct (see the bug-fix section above) —
   there is no automatic check, so an under-sized value silently gives a wrong answer with no
   warning.
3. Mode B currently writes **one set of output files per (ground-state, core-hole block)
   channel** (`trRIXS2_chan0_*`, `trRIXS2_chan1_*`, ...), unlike Mode A which writes one
   already-combined spectrum. Getting "the" final Mode B spectrum currently requires manually
   summing the channel files afterward (see `validation_modeB/compare_modeb_full.py` for how).
4. Mode B has the open magnitude-normalization bug described above, so its output isn't
   quantitatively trustworthy yet regardless of the above.

**The user has explicitly said**: they may want to implement a proper, general Mode A/Mode B
switch at some point (exposing the hardcoded parameters as real INPUT keys, auto-checking
NEV sizing, auto-summing channels) — but that is *not* the current priority. **For now, default
to Mode B** as the mode being worked on/fixed. Don't spend further effort tuning or falling back
to Mode A — it's confirmed structurally incapable of resolving the target feature (see point 3
in the validated-results section above).

## Immediate next step: normalization hypothesis FALSIFIED, new lead identified

**The "missing per-delay normalization" hypothesis queued at the end of the last session is
wrong — do not implement it.** Verified by reading the reference Fortran's own two-time-correlator
section directly (`Ni-L-RIXS/new_shape_function/XY/TD_Cluster_Xray.f90`, the `tm`/`omegai`/`omega`
loop around line 2783): `csum` is accumulated exactly like our `csum` (same `s1*s2*l12*tkern*phase`
structure, same `gauss = 1/probe_sigma` convention), and the result is written directly as
`dreal(csum4*dconjg(csum4))` — **there is no division by any per-delay weight sum anywhere in that
section.** (The unrelated equilibrium `ckern`/`spectra` block much later in the file does divide by
`twistx*twisty`, a k-point twist-average factor that has nothing to do with the two-time delay
scan — don't confuse the two.) So Mode A's `norm_weight[d]` pattern has no counterpart in the
reference for Mode B, and adding one would make our output *diverge* from the reference, not
converge. Confirmed by direct source read, not just architectural reasoning.

**Re-running `validation_modeB/compare_modeb_full.py` and reading the full loss-axis profiles
(not just the peak ratio) reveals the real pattern is more specific than "a ratio that grows":**
- `ref`'s peak intensity at `ab=-4.992,loss=2.00` **rises then falls with delay**, peaking around
  `T≈-1.88` (4.63e-5) and falling off almost 20x by `T≈+3.77` (8.49e-6) — a real physical decay
  shape, roughly symmetric-looking but stronger on the positive side.
- `mine`'s peak at the *same* `(ab,loss)` also rises then falls, but **peaks later** (around
  `T≈0` to `+1.88`, ~1.0e-3) and falls off much more gently.
- i.e. both codes show a peaked, non-monotonic delay dependence (ruling out a simple missing
  normalization *or* a pure constant/growing scale factor) — but the peak **location along the
  delay axis is shifted**, and the falloff **rate differs strongly** between the two.

**Leading explanation, not yet confirmed:** this toy validation window is small enough
(`probe_sigma=5.09`, half-window `t_hi-t_lo` only `2*6.28=12.57`) that the causal truncation
`t1<=t2<=t_hi` (both codes cap `t2` at the same simulated-range edge) removes a *different*
fraction of each delay's Gaussian-weighted (t1,t2) pairs depending on how close `T` sits to that
edge — this alone produces a rise-then-fall-with-delay shape in *both* codes, but a subtle
difference in how the truncation is realized (e.g. Krylov-propagation error accumulating over the
trajectory length being used for a given `t1`, or an off-by-one/asymmetry in how `t_lo`/`t_hi` are
placed relative to the delay grid) could shift/broaden it differently between them. **This was
previously (incorrectly) assumed to only affect the single edge-most delay (+5.65); the new data
suggests the whole delay range may be edge-affected to some degree** since `probe_sigma` is
comparable to the full window half-width.

## MAJOR FINDING: the reference toy Fortran run silently crashed after its FIRST channel —
## the "3-channel" comparison this whole session was built on compared against an incomplete reference

While setting up the v2 diagnostic test below (narrower probe, same window, to isolate edge
truncation), the rebuilt toy Fortran crashed:

```
At line 2712 of file TD_Cluster_Xray.f90
Fortran runtime error: Attempting to allocate already allocated variable 'allowed_tm'
```

Root cause: `allocate(allowed_tm(...))` sits inside the `do ppp=1,qstates` (final/core-hole block,
`wsp`) loop, guarded by `if(tspin.eq.0)`. But `tspin=0` is set unconditionally just above, and only
overwritten by `if(dlso.eq.0.0d0)tspin=sztot(ws)-sztot(wsp)` — which **never fires when SOC is on**
(`dlso != 0` for this system, `SO=11.5`). So `tspin.eq.0` is trivially true for **every** `ppp`, and
the code reallocates `allowed_tm` on every matching `ppp` without ever deallocating it. It crashes
on the *second* matching `wsp` — under `-fcheck=all` (the exact compile flags this whole project has
been using).

**This is not new to the v2 run.** Checking the exploratory run logs from earlier this session
(`/private/tmp/fortran_ref_build/run.log`, `fortran_ref_build3/run.log`) confirms the *original* toy
run that produced `validation_modeB/fortran_toy_reference/trSPEC_XY_*.dat` (the reference this whole
session's comparison was built on) **hit the exact same crash** — it just wasn't noticed because the
process still wrote a complete, plausible-looking set of 19 `trSPEC_XY_0010.dat`...`1810.dat` files
(one full tm-scan's worth) before dying, and nobody checked the exit code or scrolled to the end of
the log. Consistent with this: `fortran_toy_reference/` has no `XPS.dat` or later-stage output that
the full program would write after the RIXS section completes — it never got that far.

**What this means:** the reference only ever reflects the **first** matching (ground-state,
core-hole-block) `wsp` channel that Fortran's internal `qlist` ordering happens to try first
(dimension 55, per the log). The two other matching channels the C++ side computes and (in every
prior comparison this session) **incoherently summed in** (`chan0+chan1+chan2`, dims 55/121/55) were
**never computed by the reference at all** — the Fortran crashed before reaching them. Every
prior-session conclusion about magnitude and delay-shape mismatch (the "13×→98× growing ratio," the
"rise-then-fall shape mismatch") was drawn from comparing `mine = 3 channels summed` against
`ref = 1 channel only`. **This alone is enough to explain a large chunk of the magnitude
discrepancy**, and quite possibly the shape discrepancy too, without implicating any bug in Mode B's
physics/normalization at all.

**Fix applied (in `validation_modeB/fortran_toy_reference_v2/TD_Cluster_Xray.f90` only — the
original `fortran_toy_reference/` copy is left untouched as historical record):** added
`if (allocated(allowed_tm)) deallocate(allowed_tm)` immediately before the `allocate` call. This
lets the reference process all matching `wsp` channels and write all 3 sets of `trSPEC_XY_*.dat`
files (with different `tlabel` bases, since `tmom` differs by channel), so each can finally be
compared against its correct matching C++ channel individually, and only genuinely-matching
channels summed.

**Concrete next steps (v2 Fortran rerun in progress as of this writing, PID/log in
`validation_modeB/fortran_toy_reference_v2/run.log`):**
1. Once it completes, inspect the *new* set of `trSPEC_XY_*.dat` tlabel groups it writes (there
   should now be up to 3 distinct 5-file groups, one per surviving `wsp`, each with its own
   `100*count+10*tmom` offset — read off the actual `tmom` per group by checking which
   count-bases appear, analogous to how `tmom=1` was reverse-engineered for the original 19-file
   run) — see "How the delay/channel labels were reverse-engineered" note below for the method.
2. Match each Fortran `wsp` group to its corresponding C++ `chan{0,1,2}` by comparing block
   dimension (55/121/55) and, if that's ambiguous, by comparing which channel's peak location lines
   up with which reference group's peak location.
3. Redo the comparison **per-channel** (not summed) first — this is the real apples-to-apples test.
   Only sum channels that Fortran actually confirms are separate, nonzero, matching contributions.
4. If per-channel comparison still shows the same rise-then-fall shape mismatch even after removing
   the summed-incomplete-reference confound, then (and only then) return to the edge-truncation
   question the v2 test was originally designed to isolate (narrower probe within the same window,
   already set up and ready to compare once step 1-3 are done).
5. This also means the original `fortran_toy_reference/` 19-file dataset and every ratio/shape
   number derived from it earlier in this document should be treated as **provisional and likely
   revised** once the fixed 3-channel comparison is in hand — don't cite the old 13×-98× table or
   the "peaks at T≈-1.88 vs T≈0..+1.88" shape claim as settled fact anymore.

## What NOT to re-try (already explored exhaustively, dead ends for this specific question)

- Re-deriving Hamiltonian parameters from the Fortran source — done, cross-validated 4+ times,
  not the remaining issue.
- Mode A (frozen-snapshot) parameter tuning (`TAUCH`, `EPSLOSS`, `TRPROBESIGMA`, incident-energy
  window scans, `check_AS`/`new_shape_function` CF value differences) — extensively explored
  earlier this session; Mode A is structurally the wrong tool for this feature (see point 3
  above), further tuning within Mode A will not recover the 0.9 eV feature.
- Assuming a shrunk-window toy Fortran run is directly comparable without checking its own
  `dt`-convergence first (see the coarse-`timesteps` pitfall above).
- The "d9 at TM L-edge" configuration (Downloads file (2)) — a genuinely different, untested
  system; not related to the current d8 investigation, don't conflate the two.

## Key file locations

- `src/tr_photon.hpp`: `trRIXS_2time()` (Mode B, new), `trRIXS()`/`trXAS()` (Mode A, existing),
  `write_trRIXS_2time()`.
- `src/main.cpp`: `TRRIXS2` INPUT key parsing (~line 285), dispatch block calling
  `trRIXS_2time()` once per (ground-state, matching core-hole block) channel (~line 726-750).
- `src/photon.hpp`: `TRParam::trRIXS2` field.
- `src/tdynamics.hpp`: `propagate_block_state()` (fully generic, works on `GS` or `EX` alike —
  no changes needed here for Mode B), `build_trajectory()`.
- `validation_modeB/`: all reference data and comparison scripts described above.

## Commit state

All of this session's work (Phase 1-2 propagator, Mode A trXAS/trRIXS, the `EMLO` parameter,
pump-profile output, and Mode B including the two bug fixes) is **committed** as commit
`287a7e3` ("Add time-resolved XAS/RIXS: real-time propagator, frozen-snapshot spectroscopy,
and two-time-correlator RIXS") on `main`, including this handoff doc and everything under
`validation_modeB/`. Any *new* changes made after reading this handoff (e.g. the growing-ratio
fix) will show up as fresh working-tree modifications on top of that commit — check `git log`/
`git status` to confirm you're starting from `287a7e3` before assuming anything is uncommitted.
