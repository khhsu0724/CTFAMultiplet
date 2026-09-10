# CTFAMultiplet-magnons: current status

Read this first. It supersedes `NEXT_SESSION_REPORT.md` and `SESSION_HANDOFF.md`, which are kept
as the reasoning trail for how we got here.

Last updated: 2026-09-10, branch `fix/core-valence-multipoles-and-final-state-sum`.

---

## 1. What is validated, and against what

**Quanty 0.81 is the validation authority**, not the legacy Fortran (section 3 explains why).
Quanty lives at `~/Downloads/Quanty`; the recipe for mirroring this code's hole-language
Hamiltonian in it is in `validation_modeB/sqpl_cluster_quanty/README.md`.

| what | against | agreement | where |
|---|---|---|---|
| d^8 and 2p^5 3d^9 atomic multiplets | Quanty | all 12 levels, degeneracies identical | `validation_modeB/atomic_limit_multiplet/` |
| atomic L-edge XAS | Quanty | 11 poles to 4e-05 eV, weights to 2e-06, L3/L2 = 2.228806 both | same |
| square-planar cluster valence spectrum | Quanty | 14 levels to 0.002 eV (our Lanczos truncation) | `validation_modeB/sqpl_cluster_quanty/` |
| square-planar cluster L-edge XAS | Quanty | 120 poles to 4e-05 eV, every weight 1.0000, sum rule 1e-7 | same |
| equilibrium RIXS (KH, `solver = 4`) | Quanty | absolute scale 0.999999, r = 0.9999997 | `validation_modeB/rixs_spinflip_channel/` |
| pumped one-body Hamiltonian | legacy Fortran | 1e-14 at every pump amplitude | `validation_modeB/pump_dynamics_crosscheck/` |
| pre-pump occupations | legacy Fortran | 2e-07 | same |
| NiO equilibrium XAS/RIXS | user's stored reference | XAS bit-identical; eig/dos identical; RIXS 1.4e-4 (solver tol) | `validation_modeB/nio_equilibrium/` |

Not independently validated: **the pumped (time-dependent) RIXS itself**. Quanty cannot do it.
Its equilibrium limit is exact, its Hamiltonian and pump are exact, and its propagator is
self-converged -- but the tr-RIXS assembly has no external check. Treat new tr-RIXS physics
accordingly.

## 2. Bugs fixed in this code

All three are backward compatible; the NiO reference is reproduced identically after all of them.

1. **Inter-ligand hopping sign** (`SquarePlanar::get_tmat`). The two surviving O_x-O_y elements
   had opposite signs; Slater-Koster plus the Bloch sum requires them equal. Only bit when
   `TPPPI != 0`. Changed the valence eigenvalues by up to 0.03 eV and the pumped Hamiltonian by
   4.5e-02.
2. **RIXS dropped the spin-flip channel** -- all three paths (equilibrium `RIXS()`, Mode A
   `trRIXS_snapshot`, Mode B `trRIXS_2time`) de-excited back into the initial state's own Sz
   block. Core spin-orbit coupling mixes the core-hole spin while the intermediate state
   propagates, so Sz can change by 1. **43-54% of the L-edge cross section.** Only bites
   Sz-blocked ground sectors (valence SOC off and no `MF`); with `SO = 0 0` the fix is
   bit-identical.
3. **Mode B final-state sum was coherent** (earlier session) -- `|sum_l A_l|^2` is basis
   dependent and not an observable. Now incoherent, per PRB 99, 104306 Eq. (11)/(15).

Retracted, do not resurrect: an `FG` k-channel "fix" (the upstream `F0,G1,F2,G3` order was
already right) and an "exchange-sign fix" (the upstream sign is correct, Quanty-confirmed to
4.7e-05 eV). Both are documented as retractions in place.

## 3. Why the legacy Fortran is not a validation target

Four independent defects, in decreasing confidence:

1. **Core-hole ligand-ligand hopping sign.** Its valence sector is `xcef -2TX -2TY -4TO` and its
   core-hole sector `xcef +2TX +2TY +4TO`. The p -> -p gauge on the ligands maps the first to
   `+2TX +2TY -4TO` (verified by explicit conjugation, identical spectrum), so the two sectors are
   **not** gauge-equivalent -- they differ by the sign of the ligand-ligand hopping, worth 0.70 eV
   on the one-body levels. An O-O hopping cannot depend on whether a core hole sits on the metal.
   **This single sign reproduces its +1.9 eV XAS satellite exactly (r = 1.0000) and its ~2 eV
   post-pump red shift to 0.02 eV** -- see `validation_modeB/prx_fig8_antistokes/FIG7_RESOLUTION.md`.
   Both published features therefore look like artifacts of it. Use `EX_HOP_SIGN=-1` to reproduce
   its convention.
2. **Coherent final-state summation** -- not an observable (see 2.3 above). Explains neither the
   satellite nor the red shift.
3. **First-order time propagation.** Its occupation error against our converged result is 7.87e-2
   at `dt = 0.0785` and 1.97e-2 at `dt = 0.0196` -- ratio 4.01 for 4x finer dt, converging toward
   us. Our own propagation moves by 9e-3 over a 20x `TDDT` range.
4. **`checkint` off-by-one.** Its `ncore == 6` blocks test `u < norbspc` where every other block
   uses `u <= norbspc`, leaving the last ligand orbital interacting with the TM core hole,
   asymmetrically in index order. Small (~1 eV on 1 of 6 orbitals) and not quantified.

What the Fortran does get right, and where earlier sessions were wrong to suspect it: its
interaction tensor is **strictly on-site** (`checkint` enforces it explicitly), so there is no
inter-site monopole; and its pump and equilibrium valence Hamiltonian match ours exactly.

## 4. Known open question

**The pump-induced total-intensity decay.** The reference loses ~2.8x of its trRIXS intensity
across the pump and saturates; ours is flat with a slight rise. Measured on the converged delays
only (6, 7, 8), delay 6 -> 8: reference 0.800, ours 0.986, ours + legacy core-hole sign 0.944,
ours + coherent final states 1.157. Not truncation, not the core-hole sign, not the coherent sum
alone. Note our behaviour is what the paper's own text describes ("During the pump, there is a
slight increase in intensity... these changes are not pronounced").

## 5. Gotchas that have each cost real time

- **`solver = 1` for pole comparisons**, never a small `epsab`. `ContFracExpan` has no
  reorthogonalization, so shrinking `epsab` splits true poles into ghost peaks a few meV apart.
- **`epsab` is a Lorentzian HWHM; Quanty's `Gamma` is a FWHM.** Factor 2.
- **Mode B writes one file set per degenerate ground state.** The observable is their sum; use
  `util/sum_trrixs2_channels.py`. Reading `chan0` alone drops 2 of 3 channels for the PRX d8 set.
- **Probe-window truncation.** `sim_window_mult` now defaults to `n_sigma_probe` (4) and
  `trRIXS_2time` warns, before the expensive trajectory build, about any delay under 3 probe-sigma
  from a window edge. `TRRIXS2_WINDOW_MULT=0.75 TRRIXS2_T2_WINDOW=-1` reproduces the old
  (under-converged) validation runs bit-identically.
- **Mode B cost is `O(n_is * tau_ch/dt)`** now that `t2_max_window` defaults to `8*tau_ch`. It was
  `O(n_is^2)`; a converged 15-delay production run was ~32 h per channel under the old default.
- **RIXS is not bit-reproducible across thread counts** (~1e-4): `CG_tol = 1e-8` is tested against
  the *squared* BiCGStab residual. XAS and all eigen/DOS output are bit-reproducible.
- **Build**: the CMake target is `CTFAMultiplet` and produces `./main`; `build/main-test` needs a
  hand relink. See the memory note, or:
  `g++ -std=c++1y -O3 -ffast-math -march=native -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -I../src -c ../src/<f>.cpp`
  then link naming `libopenblas.dylib`, `libarpack.dylib`, `libomp.dylib` directly.
  `main.cpp` takes ~2 min; `tr_photon.hpp` is header-only so touching it rebuilds `main.o`.
- Don't relink `build/main-test` while runs are in flight.

## 6. Parameter sets

- `example/NiO_octahedral_Ledge/` -- the user's validated NiO set. Valence SOC + `MF = 0.11`, so a
  single 168-dim ground block. Use this as the regression target.
- `example/PRX_d8_sqpl_Ledge/` -- the PRX Fig. 7/8 d8 set, `FG = 6.0 4.635 6.18 2.63` (the
  corrected `F0,G1,F2,G3` ordering; the Fortran declares `F0,F2,G1,G3`, so F2 and G1 swap on
  translation). `SO = 11.5 0` and no `MF`, so 3 Sz ground blocks -- the configuration that exposed
  the spin-flip bug.
