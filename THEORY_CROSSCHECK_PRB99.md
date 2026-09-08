# From-first-principles cross-check of Mode B against Phys. Rev. B 99, 104306

**Status: a structural formula error was found, fixed, and independently verified.**
The two-time correlator implemented by `trRIXS_2time` (Mode B) — and, identically, the one
implemented by paper 1's own Fortran reference `TD_Cluster_Xray.f90` — **summed the scattering
amplitude coherently over the final valence states before squaring**. Chen *et al.*'s derivation
requires the opposite: the final-state sum is **incoherent** (a sum of squared amplitudes). The
coherent form is not an observable at all — it depends on the arbitrary phase the diagonalization
routine assigns each eigenvector — which is why every numerical comparison against the Fortran
reference this project has ever run produced an irregular, irreproducible discrepancy.

---

## 1. The derivation, term by term

Notation follows PRB 99, 104306 (paper 2). Single-site cluster, so the site sums `m, n` collapse
and the `e^{iq·(r_m - r_n)}` factors are 1.

### 1.1 Their starting point

Eq. (10) defines the four-time correlation function (zero temperature, `|i0>` the initial state):

```
S(t1', t2', t2, t1) = < U(-inf,t1') D_ei^dag(t1') U(t1',t2') D_ef(t2')
                         U(t2',t2) D_ef^dag(t2) U(t2,t1) D_ei(t1) U(t1,-inf) >
```

and Eq. (11) is the full tr-RIXS cross section:

```
I(wi, wf, t) = INT dt2 INT^{t2} dt1 INT dt2' INT^{t2'} dt1'
                 e^{i*wi*(t1'-t1)} e^{-i*wf*(t2'-t2)}
                 l(t1,t2) l(t1',t2') g(t1,t) g(t2,t) g(t1',t) g(t2',t)
                 S(t1', t2', t2, t1)
```

with `l(t1,t2) = exp(-|t2-t1|/tau_ch)` the phenomenological core-hole lifetime and `g` the probe
envelope. The unprimed times are the ket side, the primed times the bra side.

### 1.2 Factorizing the four-time correlator

Write `|psi(t)> = U(t,-inf)|i0>` (the pump-driven initial-state trajectory), so
`<i0|U(-inf,t1') = <psi(t1')|`. Group the operators around the *middle* propagator `U(t2',t2)`,
which acts in the **final, core-hole-free valence sector**:

```
S = < X(t1',t2') | U(t2',t2) | X(t1,t2) >
|X(t1,t2)> = D_ef^dag(t2) U_core(t2,t1) D_ei(t1) |psi(t1)>
```

`|X>` is exactly what the code builds: excite (`D_in`), propagate in the core-hole manifold
(`U_core`), de-excite (`D_out`). Note `S` is **not** `conj(A') * A` yet — the middle propagator
prevents it.

Resolve that middle propagator with any complete orthonormal basis `{|l>}` of the valence sector,
taken at a reference time `t0`:

```
U(t2',t2) = U(t2',t0) U(t0,t2) = sum_l U(t2',t0) |l> <l| U(t2,t0)^dag
```

Defining the *time-evolved* valence basis `|val_l(t)> = U(t,t0)|l>`, this gives

```
S(t1',t2',t2,t1) = sum_l conj(a_l(t1',t2')) * a_l(t1,t2)
a_l(t1,t2) = <val_l(t2)| D_out U_core(t2,t1) D_in |psi_GS(t1)>            (*)
```

Line (*) is exactly the code's `a_l`, with `t0 = t_lo`. **So the code's per-`l` amplitude is
correct.** The error is in what is done with it next.

### 1.3 The quadruple integral collapses to a sum of moduli squared

Substituting into Eq. (11), the integrand factorizes into a bra part (primed times) and a ket part
(unprimed times) *for each `l` separately*, and the two parts are complex conjugates of one
another. Therefore

```
I(wi, wloss, T) = sum_l | A_l |^2
A_l = sum_{t1<t2} g(t1-T) g(t2-T) e^{-(t2-t1)/tau_ch} a_l(t1,t2)
                  e^{-i*wi*t1 + i*wout*t2} dt^2 ,      wout = wi - wloss
```

**Coherent over `(t1,t2)`; incoherent over the final valence state `l`.**

### 1.4 Their own equilibrium limit says the same thing

Paper 2 checks its Eq. (11) against equilibrium in §II.C. Setting `H(t)=H`, `g=1`, and doing the
time integrals gives their Eq. (15)-(16):

```
I_eq = sum_f |A_f|^2 delta(w - E_f + E_i)
A_f  = <f0| D_out  1/(wi - H + E_i + i*Gamma)  D_in |i0>
```

A **sum of squared amplitudes over final states `f`**, not the square of their sum. This is the
ordinary Kramers-Heisenberg formula, and it is the equilibrium limit of the expression derived in
§1.3 above — consistent, as it must be.

### 1.5 Why the coherent form cannot be right

`A_l` is linear in the bra: `A_l = <l|M>` for a single vector `|M>` (the weighted, phase-factored
sum of the back-propagated `|X(t1,t2)>`). Then

* correct:   `sum_l |<l|M>|^2 = <M|M>` — **independent of the basis `{|l>}`**;
* as coded:  `|sum_l <l|M>|^2 = |<u|M>|^2` with `|u> = sum_l |l>` — a projection onto a direction
  determined entirely by the arbitrary sign/phase the diagonalizer assigned to each eigenvector,
  and by the arbitrary choice of basis inside every degenerate eigenspace.

Flipping the sign of a single valence eigenvector — a change no physical result may notice —
changes `|sum_l A_l|^2` by an O(1) amount.

### 1.6 Other things this cross-check confirmed as already correct

Checked against Eq. (11) line by line and found to match the code as written:

* **Time ordering / causality**: `t1 < t2` on the ket side; `U_core` runs forward from `t1` to
  `t2`; the equal-time term is excluded. Matches.
* **Conjugation**: the bra side is the complex conjugate of the ket side; no separate conjugation
  is missing.
* **Core-hole lifetime**: one factor `l(t1,t2) = exp(-(t2-t1)/tau_ch)` per side. Since the code
  builds one amplitude and squares it, exactly one factor belongs in `A_l`. Matches.
* **Probe shape**: Eq. (11) carries `g(t1,t)g(t2,t)` per side, i.e. `g(t1-T) g(t2-T)` inside `A_l`.
  Matches. (Paper 2's `g(t1)g(t2)g(t1')g(t2') ~ g^2(t2,t)g^2(t2',t)` reduction of §II.C is an
  optional simplification valid when `tau_ch << sigma_probe`; the code implements the *unreduced*
  Eq. (11), which is the more accurate choice and is what the reference Fortran does too. For the
  parameters in use, `tau_ch = 4.19 eV^-1` and `sigma_pr = 1.0-5.1 eV^-1`, that condition does
  *not* hold, so not taking the shortcut is the right call.)
* **Energy referencing** (the open `instantaneous_energy`/`E_gs` question from the previous
  session): collecting the phases of `a_l` against `e^{-i*wi*t1 + i*wout*t2}` gives a stationary
  point at `wi = E_ch - E_gs` and `wloss = E_l - E_gs`. Incident energy is already referenced to
  the ground state and loss to the final-state excitation energy, matching the equilibrium
  solver's own axes. **No extra `E_gs` term is missing.** This item can be closed.
* **Coherent sum over core-hole blocks** (`ex_blk_inds`): correct and required — different
  intermediate core-hole states are summed inside `a_l` before squaring, per Eq. (10).
* **No per-delay normalization**: confirmed; Eq. (11) has none.

---

## 2. The fix

`trRIXS_2time` in `src/tr_photon.hpp` previously summed the valence trajectories into one
`combined_traj = sum_l |val_l(t)>` and took a single inner product, producing the scalar kernel
`K(t1,t2) = sum_l a_l(t1,t2)`, then `|sum over (t1,t2)|^2`. It now keeps the `l` index open
(`Krows[is1][j*nval + l]`), accumulates `A_l` per final state in the frequency loop, and returns
`sum_l |A_l|^2`.

Two clearly-marked diagnostic environment variables were added, both no-ops unless set:

* `TRRIXS2_COHERENT_FINAL` — restores the old (reference-Fortran) coherent convention, for
  diffing against that reference.
* `TRRIXS2_FLIP_VAL=<l>` — sign-flips valence eigenvector `l`, to exercise the basis-invariance
  test below.
* `TRRIXS2_WINDOW_MULT` — overrides the hardcoded `sim_window_mult=0.75` (which exists only to
  match the toy Fortran's bounded window), needed for the equilibrium-limit check.
* `RIXS_GS_INDEX` (in `src/photon.cpp`) — restricts the *equilibrium* RIXS solver to one
  degenerate ground state, so it can be compared like-for-like against a single Mode B channel.

A warning is now emitted if `GSNEV` is too small to give a complete valence eigenbasis, since the
final-state sum silently requires `sum_l |l><l| = 1`.

---

## 3. Numerical verification (toy d8 case, `validation_modeB/modeb_pumpon_v2_fixed/INPUT`, chan0)

### 3.1 Basis-invariance test — the decisive one

Same binary, same INPUT, one environment variable different (`TRRIXS2_FLIP_VAL=3`, i.e. the sign
of one valence eigenvector flipped). Max relative change in the output intensity:

| delay (eV^-1) | corrected `sum_l \|A_l\|^2` | as-was `\|sum_l A_l\|^2` |
|---|---|---|
| -3.14 | 9.2e-04 | **5.7e-01** |
|  0.00 | 2.4e-03 | **1.7e+00** |
| +3.14 | 4.3e-03 | **7.4e-01** |

The corrected result is invariant to within its own numerical noise. The as-implemented result
changes by 57-172%.

### 3.2 Run-to-run reproducibility

The `l`-summed kernel `sum_l a_l(t1,t2)` is not reproducible between runs of the same binary on
the same input — the diagonalizer does not always return the same basis inside the degenerate
valence eigenspaces. Two runs whose kernels differed by a factor of ~3 in magnitude
(`(0.0160,-0.0355)` vs `(0.0059,-0.1590)` at `t1=0, t2=dt`) nevertheless produced **the same
corrected intensity to 7e-4 - 2e-3**, while the coherent convention's intensity differed between
runs by **80-97%**. This is the exact signature of a basis-dependent quantity: run-to-run
variation is a rotation inside a degenerate eigenspace, which `sum_l |A_l|^2` cannot see and
`|sum_l A_l|^2` is entirely at the mercy of.

**This also retroactively explains the whole previous session's failed comparisons.** The
0.7x-575x, delay-dependent, feature-dependent, irregular-phase discrepancy against the Fortran
reference was never a physics bug to be localized: both codes were computing a quantity that is
not reproducible even against itself. The recorded raw-kernel mismatch (2.8-4.8x magnitude,
irregular phase, not growing with propagation time) is exactly what comparing two arbitrary
projections of the same vector looks like.

### 3.3 Positive check: the equilibrium limit

Independently of paper 2's algebra, the equilibrium RIXS path already in this codebase
(`RIXS()` in `src/photon.cpp`, `solver = 4`) evaluates the loss axis with a continued-fraction
expansion of `Im <v| (w + E_gs - H + i*eps)^-1 |v>`, which *is* `sum_f |<f|v>|^2` broadened —
the incoherent final-state sum of Eq. (15). Mode B was the only place in the repository doing it
the other way. That makes it a genuine, first-principles regression target: **with the pump off,
Mode B must reproduce the equilibrium Kramers-Heisenberg spectrum.**

Run with `PUMPAMPX = 0`, the same single ground state on both sides (`RIXS_GS_INDEX=0`),
`EPSAB = 1/tau_ch = 0.2387` matched to the two-time core-hole decay, and the equilibrium loss axis
Gaussian-broadened to the probe-limited resolution of the Mode B run:

**Incident-energy profile** (summed over loss, normalized), pump-off Mode B vs equilibrium KH:

```
  ab       -6.9  -5.0  -3.1  -1.2   0.8   2.7   4.6   6.5   8.4  10.3  12.2  14.2 ...
  EQ (KH)   0.05  1.00  0.31  0.09  0.01  0.00  0.00  0.00  0.00  0.01  0.23  0.56
  Mode B    0.47  1.00  0.31  0.09  0.03  0.01  0.01  0.01  0.03  0.16  0.35  0.15
```

Same resonance at `ab = -4.992`, same falloff. **The absorption resonance and its energy
referencing are confirmed correct from first principles**, independently of the Fortran.

**Discrimination between the two conventions**, Pearson correlation of the full 2D map against the
probe-broadened equilibrium spectrum:

| loss-axis blur | corrected `sum_l \|A_l\|^2` | as-was `\|sum_l A_l\|^2` |
|---|---|---|
| 0.6 eV | **+0.791** | +0.700 |
| 0.8 eV | **+0.794** | +0.700 |
| 1.0 eV | **+0.789** | +0.692 |
| 1.2 eV | **+0.780** | +0.682 |

The corrected form is consistently closer to the equilibrium limit. The residual gap is expected at
this scale: the toy window's probe (`sigma_pr = 1.0 eV^-1`) gives ~1 eV energy resolution against a
0.5 eV loss grid, and modelling the probe broadening as a simple Gaussian convolution of the
equilibrium spectrum is itself only approximate.

**Setup caveat learned here, worth recording:** an attempt to sharpen this test by widening the
probe to `sigma_pr = 3.0` with `TRRIXS2_WINDOW_MULT=2.865` (window `+/-12`) made the agreement much
*worse*, moving the apparent resonance to the incident-grid edge. That is a hard-truncation
artifact of the setup, not a code bug: window `+/-12` is only `4*sigma_pr`, where the probe
Gaussian is still at `e^-8 ~ 3e-4` rather than numerically zero. The toy setting
(`sigma_pr = 1.0`, window `+/-6.28` = `6.3*sigma_pr`) is properly converged. **Any Mode B run needs
its simulated window to be at least ~6 probe sigmas wide**, which the hardcoded
`sim_window_mult = 0.75` does not guarantee in general (it is `0.75*max(sigma_pr, tau_ch)`, and is
tau_ch-dominated here).

---

## 4. What this does and does not settle

**Settled:** the final-state sum. The formula `trRIXS_2time` now implements matches Eq. (11)
term by term, the result is basis-independent as a cross section must be, and with the pump off it
reproduces this code's own equilibrium Kramers-Heisenberg spectrum better than the old form did.

**Not settled by this work:**

1. **The paper-1 Fortran reference has the same error** (`do l=1,min(nev,binstates(wsp))`
   accumulating into one `csum`, squared once after the frequency sum, around line 2783 of
   `TD_Cluster_Xray.f90`). Agreement with that reference is therefore no longer a meaningful
   validation target for this quantity, and disagreement with it is not evidence of a bug here.
   Any future Mode B validation should go against the equilibrium limit, or against a corrected
   reference. **This is worth raising with the paper-1 authors** — it affects published Fig. 7/8
   intensities, though not necessarily peak positions.
2. **A real OpenMP nondeterminism in the propagation path.** Two identical *multithreaded* runs of
   the same binary on the same input differ by ~3% typically and up to 20% locally in the final
   intensity (2D Pearson r = 0.9989); run **single-threaded with `OMP_NUM_THREADS=1
   OMP_DYNAMIC=FALSE` and they are bit-for-bit identical**. The source is in
   `propagate_block_state`'s hot path. **A `sample` profile settles what is going on, and it is
   worth knowing: the fine-grained `#pragma omp parallel for` pragmas inside the propagation are
   pure overhead.** Of the samples inside `Dense::mvmult_cmplx`, essentially all are
   `__kmpc_fork_call` -> `__kmp_join_barrier` -> `pthread_cond_wait` -- thread fork/join around
   loops of only ~55 elements. Measured end-to-end on the production Hamiltonian,
   `OMP_NUM_THREADS=1` runs **3.75x faster** than `OMP_NUM_THREADS=4` (851 ms vs 3191 ms), and is
   bit-reproducible. So Mode B should be run **single-threaded, with concurrency across separate
   processes**; the right place for OpenMP here is a coarse grain (per delay, per incident
   energy), not inside a 55-element vector loop. (`calc_HYB`, which rebuilds the whole-Hilbert
   hybridization every timestep, looks like an obvious suspect but the profile puts it at only
   ~4% -- `LanczosPropagate` is 77%.) Note what this does to the two conventions: across the same
   pair of runs the corrected observable moves by ~3%, while the basis-dependent kernel
   `sum_l a_l(0,dt)` moves by O(1) -- `(-0.148,-0.334)` vs `(0.0068,-0.160)`, against the
   single-threaded, reproducible `(0.4297,-0.1444)`. Any run whose numbers must be exactly
   reproducible should be single-threaded until this is fixed -- which is also the faster option.
3. **The equilibrium-limit agreement is good, not exact** (r ~ 0.79). This is now the sharpest
   available lead and the right regression test to sharpen: rerun it at production resolution
   (a wider probe with a window of at least ~6 probe sigmas, and a finer loss grid) and see
   whether the residual gap closes. If it does not, there is a second, smaller problem still in
   Mode B — and unlike every comparison in `SESSION_HANDOFF.md`, this test does not depend on the
   Fortran reference at all.
4. **Whether this fixes the target physics** — the PRX Fig. 7/8 Hund's-exchange anti-Stokes
   feature at ~-0.9 eV — is untested at production scale. The toy validation window
   (`sigma_pr = 1.0 eV^-1`) has ~1 eV energy resolution, far too coarse to resolve it.
