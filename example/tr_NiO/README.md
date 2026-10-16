# tr_NiO — time-resolved XAS and RIXS on a NiO cluster

Worked example for the pump-probe machinery: a single-site octahedral NiO cluster (d8, two
valence holes) driven by an optical pump and probed at the Ni L3 edge.

The pump enters through a **Peierls substitution on the hybridization only** — the hopping
phase picks up `k -> k + A(t)`, with

    A(t) = PUMPAMPX * exp(-(t-PUMPT0)^2 / 2 PUMPSIGMA^2) * sin(PUMPFREQ * (t-PUMPT0))

All times are natural units, hbar = 1, so time is in eV^-1 (1 eV^-1 = 0.658 fs).

## The four inputs

| file | what it does | rough cost |
|---|---|---|
| `INPUT_equilibrium` | static XAS + RIXS, no pump. The reference the time-resolved runs reduce to. | seconds |
| `INPUT_occupation_scan` | orbital hole occupations vs time under the pump (`TDOCCSCAN`). The cheapest way to see whether a pump setting does anything. | seconds |
| `INPUT_trXAS` | time-resolved XAS over a delay series (`TRXAS`). | ~2 min |
| `INPUT_trRIXS2` | time-resolved RIXS via the two-time correlator (`TRRIXS2`). | ~1 h |

Run each by copying it to `INPUT` in a working directory and running the binary there.

## Input keys this adds

**Pump** — all default to zero/off, so an input without them behaves exactly as before.

| key | meaning |
|---|---|
| `TIMEDEP` | master switch for the real-time machinery |
| `PUMPFREQ` | pump angular frequency in eV (**not** cycles — `pi` means omega = pi eV) |
| `PUMPAMPX/Y/Z` | Peierls-phase amplitude, dimensionless radians. Convert to field/fluence with `util/pump_units.py` |
| `PUMPSIGMA` | Gaussian envelope width in eV^-1 |
| `PUMPT0` | pulse centre |
| `TDDT`, `TDM` | propagator step and Krylov dimension |
| `TDTMIN`, `TDTMAX` | time window for `TDOCCSCAN` |
| `MF` | molecular field on the d spin (eV), splitting the Ms sublevels |
| `TPPPI` | ligand-ligand pi bonding (square planar); defaults to 0, reproducing the previous `tpp`-only form |

**Spectroscopy**

| key | meaning |
|---|---|
| `TRXAS` | time-resolved XAS |
| `TRRIXS` | time-resolved RIXS, frozen-Hamiltonian snapshot (Mode A) |
| `TRRIXS2` | time-resolved RIXS, two-time correlator (Mode B) — PRB **99**, 104306 |
| `TRDELAY` | pump-probe delays: `start end count` |
| `TRPROBESIGMA` | probe width in eV^-1 |
| `TRDTPROBE` | probe sampling step (defaults to `TDDT`) |
| `TAUCH` | core-hole lifetime in eV^-1. A **physical** lifetime, deliberately independent of `epsab`/`epsloss`, which are numerical broadenings |
| `EMLO` | lower bound of the loss axis. Defaults to `-2`, the value previously hardcoded, so existing inputs are unaffected. Set it negative to resolve anti-Stokes (energy-gain) features |

## Things that will bite you

**`TRPROBESIGMA` sets the energy resolution, not `epsloss`.** In Mode B the loss-axis width
is the Fourier transform of the probe envelope, FWHM ~ `2.25/TRPROBESIGMA`. `epsloss` does
nothing there. Resolving 100 meV features needs `TRPROBESIGMA ~ 20`, which is expensive.

**Include a delay that precedes the pulse.** Every time-resolved path starts its trajectory
from the static ground state at `t_lo = min(TRDELAY) - 4*TRPROBESIGMA` and propagates
forward. If `t_lo` falls after the pulse, the pump is never applied and you silently get the
equilibrium answer for every delay and every amplitude. The code warns when it detects this,
but the fix is yours: keep an early delay in the list, as `INPUT_trRIXS2` does with `-5.0`.

**Mode B writes one file set per degenerate ground state.** `trRIXS2_chan<n>_*` is one
channel, not the cross section; the observable is their sum. Use
`util/sum_trrixs2_channels.py`.

**Mode B needs a complete eigenbasis.** The final-state resolution of identity is only
complete if `GSNEV`/`EXNEV` cover the whole block, which is why they are large here. The code
warns if they are not.

**The shipped pulse is sub-cycle.** `PUMPFREQ = 0.25` with `PUMPSIGMA = 1.2247` gives
`omega*sigma = 0.306`, i.e. an envelope ~0.115 of an optical cycle, with a bandwidth 3.3x the
centre frequency. It is a convenient numerical kick, not a realizable 5 um pulse.
`util/pump_units.py` converts an amplitude to a peak field and fluence and reports this; see
`util/PUMP_UNITS.md` for the derivation.
