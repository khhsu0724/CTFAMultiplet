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

**This turned out to be a dead end in an informative way.** After fixing the crash, the fixed
Fortran run (a 3rd, then 4th/5th debug-instrumented attempt, all in
`validation_modeB/fortran_toy_reference_v2/`) revealed that `wsp` (final valence sector) only ever
takes **one** value in the `do ppp=1,qstates` loop before the code moves to the next polarization
pair — not because of a bug, but because of a real, hard `if(tspin.eq.0)` selection-rule gate:
`tspin = sztot(ws) - sztot(wsp)`, and debug prints confirmed `tspin = 0, -2, -4` for `ppp = 1, 2, 3`
respectively. Only `ppp=1` (Δ(valence Sz)=0) passes. **This is not a Fortran limitation — it is the
reference correctly enforcing that the electric dipole operator does not touch spin, so a valence→
core-hole→valence round trip must return the electron to a state of the same total Sz whenever the
valence sector itself has no spin-orbit coupling mixing Sz** (confirmed: `dlso` is this system's
*valence* SOC, printed as `Valence spin-orbit coupling: 0.00000000` — the `SO=11.5` in the INPUT is
the *core* 2p SOC, only relevant to the core-hole state, not to whether valence Sz is a good quantum
number).

## ROOT CAUSE FOUND: Mode B incoherently sums exblk (core-hole block) channels that should be summed coherently

The natural next question was: does the C++ side make the same mistake, treating `exblk_ind` as a
loop over *final valence* sectors and failing to filter by Sz? **No** — checking
`trRIXS_2time`'s own construction (`src/tr_photon.hpp`) shows `combined_traj` is built from
`gsblk = GS.hblks[gs_blk_ind]` alone, the *same single ground-state block* for every value of
`ex_blk_ind` in the `chan0/1/2` loop. So the C++ implementation already restricts the **final
valence state** to the same Sz sector as the initial ground state, for all 3 channels — exactly
matching Fortran's `ppp=1`-only conclusion. The bug is not there.

The actual difference: `exblk_ind` in the C++ loop enumerates **intermediate core-hole blocks**
(`EX.hblks`), and `chan0/1/2` correspond to 3 *different* core-hole symmetry blocks that are all
individually dipole-connected to the same ground state. Checking their actual computed intensities
confirms all three are genuinely nonzero and comparable in magnitude (at delay=0: chan0 max=7.8e-3,
chan1 max=6.9e-3, chan2 max=6.1e-3) — not near-zero noise.

The reference has no equivalent split at all: its entire core-hole manifold is **one undifferentiated
132-dimensional space**, and its `csum` vertex sum (`do k=1,gdegen(ws); do l=...; do i=1,min(nevn,
cstates(momgs)); csum=csum+dconjg(tfvert(wsp,l,i,t2))*tpivert(ws,k,i,t2,t1); end do...`) coherently
sums the amplitude over **every** core-hole eigenstate `i` in that manifold into one complex number
*before* `tkern(t1,t2)=csum` is stored — squaring (`|csum|²`) only happens once, much later, after
the full `t1,t2` Fourier sum. Physically this must be right: core SOC (`SO[0]=11.5`, present in the
core-hole Hamiltonian) genuinely mixes what would otherwise be separate Sz sub-blocks within the
core-hole manifold — the manifold is not physically separable into independent channels, so summing
per-block **intensities** (`|amplitude_block0|² + |amplitude_block1|² + |amplitude_block2|²`, what
`main.cpp`'s `trRIXS2_chan{0,1,2}` dispatch currently does, with the 3 files meant to be summed
*after* squaring — see `validation_modeB/compare_modeb_full.py`) is not the same thing as summing per-
block **amplitudes and squaring once** (`|amplitude_block0 + amplitude_block1 + amplitude_block2|²`)
whenever cross-block interference terms are nonzero. Because the blocks are literally different
`EX.hblks` entries, their propagated states are also orthogonal at fixed time (no cross terms in an
inner product against a valence deexcitation vector *at one instant*) — but the interference that
matters here is between the SAME (t1,t2) pair's contributions from different blocks reaching the
SAME final valence eigenstate via DIFFERENT intermediate paths, which is exactly what coherent
amplitude summation captures and incoherent intensity summation destroys.

**This is the same category of bug already found and fixed once this session on the valence side**
(the `combined_traj` fix, coherently summing the de-excitation amplitude over all valence
eigenstates before the K-row inner product) — it was simply never applied to the intermediate
core-hole side. This is very likely the dominant remaining cause of both the magnitude and
delay-shape discrepancies chased through this entire session; the toy validation window's edge-
truncation questions (the original motivation for the whole v2 test) are almost certainly a much
smaller secondary effect by comparison, and should only be revisited after this fix is in place and
verified.

**Concrete fix (not yet implemented as of this writing):** restructure `trRIXS_2time` (and its
`main.cpp` TRRIXS2 dispatch caller) so that, for a fixed `(gs_blk_ind, gs_eig_ind)`, the `Krows`
complex-amplitude array is accumulated **coherently across every matching `ex_blk_ind`** — i.e. call
the K-row-building logic once per `exblk`, but add each `exblk`'s contribution into one shared
`Krows` (or an equivalent combined array) *before* the final `intensity[d][...] = std::norm(csum)`
step, rather than writing separate `trRIXS2_chan{N}_*.txt` files intended to be summed after the
fact. `write_trRIXS_2time` should then be called once with the single combined result. This changes
the loop structure in `main.cpp`'s TRRIXS2 dispatch block (currently `for (auto& exblk : EX.hblks) {
trRIXS_2time(...) }`, one independent call per channel) to something that passes all matching
`exblk`s into a single `trRIXS_2time` call, or accumulates their `Krows` externally before the final
delay/frequency sum.

**Implemented** (`trRIXS_2time` in `src/tr_photon.hpp` now takes `const std::vector<size_t>&
ex_blk_inds` and accumulates every matching core-hole block into one shared `Krows` before the final
`|.|^2`; `main.cpp`'s TRRIXS2 dispatch collects matching `exblk`s per ground state first, then makes
one `trRIXS_2time` call per ground state instead of per exblk). Compiles cleanly. **Keep this fix —
it is still correct and necessary in general** (multiple core-hole blocks legitimately reachable
from one ground state must be summed as amplitudes, not intensities) — but see below: for *this
specific toy system*, it turned out to make **no numerical difference**, because the premise that
motivated it (3 exblk matches per single ground state) was wrong for this system.

## CORRECTION: "3 channels" were never 3 exblk per ground state — they are 3 degenerate ground states

Rerunning the v2 toy test with the coherent-exblk fix in place still produced **3** separate output
files (`trRIXS2_chan0/1/2`), not 1. Checking `validation_modeB/modeb_pumpon_v2_fixed/run.log`
directly shows why: each of the 3 channels prints `accumulating core-hole block 0...` exactly **once**
— never block 1 or 2. So each channel's `ex_blk_inds` list only ever had one entry to begin with; the
coherent-sum fix had nothing to combine here. The channel count comes from the *outer* `for (auto& g
: gsi)` loop instead: `gsi.size() == 3`, with valence block dimensions 55, 121, 55 — the same three
`2Sz = -2, 0, +2` sectors seen throughout this investigation, all sitting at the ground-state energy.
**This is a spin-triplet (S=1) ground state, and "chan0/1/2" are its three degenerate Sz sublevels**,
not three core-hole intermediate states. Every "3-channel" comparison this entire session — the
original 13×-98× ratio table, the shape-mismatch analysis, the incomplete-Fortran-reference
discovery, the exblk coherent-sum investigation — was actually comparing (sum or individual pick of)
3 degenerate *initial states* against the reference's output, not 3 intermediate-state channels for
one initial state. The exblk-coherence question this section was originally about is a real, orthogonal
issue that simply doesn't happen to arise in this particular toy configuration.

**Why this matters:** the reference Fortran's `do qqq=1,qstates; ...; if(count2.eq.0)` structure
processes **only the first** `ws` satisfying `gsec(ws)==1` — i.e. it computes RIXS from *one single*
representative sublevel of the degenerate ground state, not a properly averaged (sum or mean) result
over all 2S+1=3 sublevels. Comparing our sum-over-3 (or any arbitrary one of our 3) against the
reference's single, specific, unidentified sublevel is not guaranteed to match even if both codes'
underlying physics is completely correct — the two quantities are conceptually different observables
unless the specific sublevel picked happens to equal (or be simply proportional to) an average.

**Tested: none of the 3 individual C++ channels cleanly matches the reference's single channel.**
Comparing `trRIXS2_chan{0,1,2}_t{delay}.txt` (in `validation_modeB/modeb_pumpon_v2_fixed/`)
individually against `trSPEC_XY_1010...1410.dat` (in `validation_modeB/fortran_toy_reference_v2/`)
across all 5 delays: peak locations disagree even at `delay=0` (reference's actual global peak moved
to `ab=12.240, loss=0.50`, a completely different resonance than the `ab=-4.992, loss=2.00` feature
this whole session has tracked — that feature is still present in the reference at `delay=0`, just
subdominant there, at `I=3.33e-4` vs. the new peak's `9.61e-4`), and the growing-with-delay ratio
pattern reappears in every one of the 3 individual-channel comparisons (roughly 3×-150× depending on
channel and delay) — i.e. **the core delay-shape/magnitude discrepancy this whole session has been
chasing is still unresolved**, and is now known to be independent of both bugs fixed so far (the
Fortran crash/tlabel-collision bugs, and the exblk-coherence question).

**Resolved: the reference processes `ws=1`, `2Sz=-2`, `binstates=55`** — confirmed via a targeted
debug print (`DEBUG which ws processed: ws= 1 2Sz= -2 binstates= 55`, consistent across every
polarization pair it reaches, in `validation_modeB/fortran_toy_reference_v2/run.log`). `GS.hblks` is
built in ascending-Sz order (`hilbert.cpp`'s `sz_Hash`/block-construction loop runs
`for (int sz = -max_2sz; sz <= max_2sz; sz+=2)`), and `main.cpp`'s `gsi` is populated by iterating
`GS.hblks` in index order, so **`gsi[0]` = block 0 = `2Sz=-2` = `chan0`** is the *correctly matched*
C++ channel — no ambiguity, confirmed independently by dimension (55) matching `binstates=55`.

## FINAL RESULT (as of this writing): a real, unconfounded discrepancy remains even for the correctly-matched channel

Comparing `chan0` (only) against the reference across all 5 delays, at the historically-tracked
`ab=-4.992, loss=2.00` feature:

| delay (eV⁻¹) | ref | mine (chan0) | ratio |
|---|---|---|---|
| −3.14 | 1.074e-03 | 6.495e-03 | 6.0× |
| −1.57 | 1.894e-03 | 1.342e-03 | 0.7× |
| 0.00 | 3.327e-04 | 7.748e-03 | 23× |
| +1.57 | 3.809e-04 | 8.007e-03 | 21× |
| +3.14 | 3.343e-05 | 1.924e-02 | 575× |

This is now a **clean, unconfounded result** — every previously-suspected confound has been
eliminated: not a missing normalization (falsified by direct source read), not an incomplete/crashed
reference (fixed, verified clean completion), not incoherent exblk summation (fixed, verified it
doesn't apply to this system), not a degenerate-ground-state mismatch (resolved, this is the single
correctly-matched channel). The ratio is not constant, not simply the old "13×→98× growing" pattern
either — it's noisier (dips below 1× at one delay) but still clearly blows up badly at the largest
tested delay.

**Also notable and likely diagnostic:** the reference's *dominant global peak* location itself moves
around with delay — at `delay=0` it jumps to a completely different resonance at `ab=12.240, loss=
0.50` (`I=9.61e-4`, larger than the `ab=-4.992` feature's `3.33e-4` there), and at `delay=+3.14` to
`ab=10.325, loss=0.50` (`I=1.39e-4`). Our `chan0` output stays dominated by the `ab=-4.992` feature
at every delay and never develops a comparably strong peak near `ab=10-12` — i.e. **the reference
appears to dynamically activate an entirely different spectral feature at certain delays that our
Mode B output does not reproduce at all**, not just get the magnitude of the same feature wrong. This
is a much sharper, more specific lead than anything found earlier in this document and should be the
starting point for the next investigation: what physical resonance sits at `ab≈10-12` (a different
incident energy channel entirely — check the INCIDENT grid in the INPUT, `-6.9065523978701657
21.813222178495966 16`, to see what physical transition `ab≈10-12` corresponds to), and why would it
be delay-dependent in the reference but not (or much less so) in ours? Candidate explanations not yet
tested: (a) this could be a *different* final valence eigenstate within the *same* `gsblk`/`2Sz=-2`
block that our `combined_traj` construction should already include (it sums over *all* `nval` valence
eigenstates in the block) — if so, the discrepancy would have to be in the K(t1,t2) amplitude
reaching that particular final state, not in which states are summed; (b) a genuine numerical
convergence issue in the Lanczos propagation (`hparam.td_m`) or timestep (`hparam.td_dt`) for this
higher-energy-loss channel specifically; (c) still-unruled-out edge-of-window truncation, though the
v2 test's delays were deliberately kept ≥3 probe-sigma from the boundary specifically to avoid this.

**Checked: the `ab≈12.24` feature is NOT missing from our output** — it's present in `chan0` at every
delay, at a comparable order of magnitude to the reference (e.g. at `delay=0`: ref=`9.61e-4`,
mine=`1.166e-3`, only a **1.2× ratio** — dramatically closer than the `ab=-4.992` feature's 23× at
the same delay). This reframes the diagnosis substantially: **it is not that we're missing a
transition/pathway. It's that our `ab=-4.992, loss≈2` (dd-excitation) feature specifically is
systematically too large by up to ~2 orders of magnitude, badly enough that it drowns out an
`ab≈12.24` feature that our code otherwise gets roughly right.** In the reference, at delays where
the dd-excitation feature is naturally weak (its own physical `T`-dependence, whatever drives that),
the `ab≈12.24` feature becomes the visible global max simply because the (correctly-sized) dd feature
dropped below it — an effect entirely masked in our output because our dd feature never drops enough
to be overtaken.

**Tested and ruled out: Lanczos/Krylov propagator convergence.** Reran the identical v2 `chan0` test
with `TDM` (Krylov subspace order) raised from 20 to 50 — near-exact diagonalization for the 55-dim
blocks involved. Result: **bit-for-bit identical** intensity at every one of the 5 delays (to the
digits printed). The propagator is already fully converged at `TDM=20`; this is not a numerical
accuracy issue and `TDM` is not the fix.

## DECISIVE TEST: raw K(t1,t2) amplitude compared directly between both codes — real error found, but far too small to explain the final discrepancy alone

Since `Krows` (our K(t1,t2)) is built once and shared across every incident energy, while only the
final delay/frequency loop depends on `win`/`wout`, the natural next test is to compare the **raw
two-time-correlator amplitude itself** — before any Fourier/frequency processing — directly against
the reference's `tkern(t1,t2)`. Both codes were instrumented to print `K(t1=0, t2=k·dt)` for
`k=1..5`, keeping everything else in the v2 `chan0` setup fixed. Result:

| k | ref (Re, Im) | mine (Re, Im) | \|mine\|/\|ref\| | phase diff (°) |
|---|---|---|---|---|
| 1 | (0.022965, −0.111863) | (−0.00609, −0.41457) | 3.63 | −12.4 |
| 2 | (−0.062462, −0.119060) | (0.13246, −0.34793) | 2.77 | +48.5 |
| 3 | (−0.099491, −0.019005) | (0.26449, −0.13095) | 2.91 | +142.9 |
| 4 | (−0.044650, 0.070359) | (0.30492, 0.11940) | 3.93 | −101.0 |
| 5 | (0.037773, 0.076761) | (0.25392, 0.27186) | 4.35 | −16.8 |

**A real discrepancy exists in the raw amplitude construction** — magnitudes are systematically too
large by a factor of roughly **2.8–4.3×** (not 1.0×, ruling out "amplitude construction is exactly
correct") — but this is drastically smaller than the up-to-575× discrepancy seen in the final
frequency-domain output. **The phase differences do not follow a clean pattern** (not a constant
offset, not a clean conjugate/mirror, not a smoothly-growing linear drift when unwrapped — differences
between consecutive k's phase-diffs are 61°, 94°, 116°, 84°, suggestive of *something* systematic but
too noisy to pin to a single constant frequency/energy offset from only 5 points).

**Conclusion: the bug is NOT purely in the final delay/frequency summation loop** (a modest 3-4×,
non-unity magnitude ratio proves the raw amplitude itself is measurably wrong) **but the final-output
discrepancy is far larger than this raw error alone can explain** — meaning the final coherent
(t1,t2) double-sum used to build a *sharp, long-lived resonance* (the dd-excitation feature
specifically) amplifies a modest point-wise phase/magnitude mismatch into a much larger error,
consistent with how narrow resonances are especially sensitive to phase coherence across a wide time
range. **Both a real amplitude-construction bug and further amplification in the coherent
(t1,t2)-sum are likely both part of the picture** — fixing the amplitude-construction bug is the
priority since it's the more tractable, better-localized piece.

**Where to look next:** the `~3-4×` magnitude error and irregular phase pattern point toward the
dipole matrix element / propagation pieces feeding `Krows`, not the frequency-domain loop. Concrete
candidates, cheapest-to-check first:
1. **Compare `blap_in`/`blap_out` (dipole overlap coefficients from `basis_overlap`) directly against
   Fortran's own dipole matrix elements** for this specific `(gs_blk_ind, ex_blk_ind)` pair — a
   normalization-convention mismatch here (e.g. a missing `sqrt` factor, a different Clebsch-Gordan/
   reduced-matrix-element convention, or an extra factor from `pm.pvin`/`pm.pvout` handling) would
   show up as exactly this kind of moderate, not-quite-constant magnitude/phase discrepancy once
   summed over several matrix elements.
2. **Check whether Fortran's `tpivert`/`tfvert` construction includes an energy-reference term ours
   doesn't** — an old code comment in the reference (`confirmed tpivert(tj,ti)~exp(-i(tj-ti)*E_CH(m)
   -i(ti-t0)*E_gs)`) explicitly ties the propagated amplitude's phase to *both* the core-hole
   eigenvalue *and* `E_gs` (ground energy) in a specific way. This echoes a concern the user raised
   much earlier this session about whether RIXS energies should be referenced to the ground state or
   the pump-excited state (`instantaneous_energy`, used explicitly in Mode A's `trRIXS`/`trXAS` to
   anchor the absorption resonance to `E_gs_t1`, an *instantaneous* expectation value under the
   possibly-pump-driven Hamiltonian) — worth checking carefully whether `trRIXS_2time` needs an
   analogous, currently-absent energy-reference term somewhere in the K(t1,t2) construction (not
   necessarily in the final `win`/`wout` phase factor, since peak *locations* already agree between
   codes, ruling out a gross missing offset there — but possibly in how the propagation phase relates
   `E_gs` to the core-hole eigenvalues internally).
3. Only after (1)-(2) are checked: revisit whether `combined_traj`'s treatment of near-degenerate
   final valence eigenstates (several of the eigenvalues in this block are 2-4 fold degenerate, per
   the printed spectrum — e.g. `11.9000044` appears with degeneracy 3 and 4 in adjacent slots) could
   introduce spurious coherent cross-terms; ruled out as a code-vs-code *difference* (Fortran's own
   `l` loop treats degenerate members identically, individually, coherently — same as ours) but not
   yet ruled out as a genuine subtlety in whether coherent summing over degenerate final states is
   the physically correct treatment in the first place.

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
