# Calibrating the pump amplitude against experiment, via the 3d9 L3 shift

The 3d9 L3 trXAS shift is **exactly** minus the permanently deposited energy: the core-excited
sector is `2p^5 3d^10` with no valence holes, so `E_EX` is pump-independent and the line sits at
`E_EX - E_GS(t)` (see `README.md`). That makes the shift a direct, assumption-free readout of how
hard the pump is driving the cluster -- and therefore a way to calibrate `PUMPAMPX` against a
measured shift.

Experimental reference point (unpublished; details withheld): a Ni L-edge trXAS red shift of
order 100 meV under an optical pump, measured on a square-planar, nominally 3d9 Ni(1+) system
for which the 3d9 cluster here is a reasonable model.

## At the PRX pump frequency, omega = pi eV

| `PUMPAMPX` | Fortran `pump_ampx` | deposited energy | L3 red shift |
|---|---|---|---|
| 0.25 | 0.5 | 0.0003 eV | 0 meV |
| 0.5 | 1.0 | 0.0043 eV | 4 meV |
| 0.75 | 1.5 | 0.0226 eV | 23 meV |
| **1.0** | **2.0** | **0.0704 eV** | **70 meV**  <- matches experiment |
| 1.5 | 3.0 | 0.2440 eV | 244 meV |
| 2.0 | 4.0 | 0.2538 eV | 254 meV |
| 3.0 | 6.0 | 0.2044 eV | 204 meV |
| **5.0** | **10.0** | **1.1882 eV** | **1188 meV**  <- the value used in the paper |

## At the experiment's pump frequency, omega = 4.75 eV

| `PUMPAMPX` | deposited energy | L3 red shift |
|---|---|---|
| 2.0 | 0.0084 eV | 8 meV |
| **~2.8** | -- | **~70 meV** (interpolated) |
| 3.0 | 0.1002 eV | 100 meV |
| 4.0 | 0.1265 eV | 126 meV |
| 5.0 | 0.1156 eV | 116 meV |
| 6.0 | 0.4457 eV | 446 meV |
| 10.0 | 2.8127 eV | 2813 meV |

4.75 eV couples much more weakly than `pi` at the same amplitude: it is off-resonance for this
cluster (`Delta = 2.8`, relevant excitations ~1.5-3 eV) and the pulse is short
(`PUMPSIGMA = 1.2247` -> bandwidth ~0.8 eV), so there is little spectral overlap.

## Two things to take from this

1. **The paper's pump is roughly 5x stronger than the experiment corresponds to.** At
   `PUMPAMPX = 5` the cluster absorbs 1.19 eV per hole and the L3 line moves 1.19 eV -- about 15x
   the measured ~100 meV. That is by design: the paper is studying multilevel Rabi oscillations,
   which require strong driving. But it means **neither** our 3d8 trRIXS shift (-0.20 eV) nor the
   legacy code's (-2.03 eV) is in the experimentally relevant regime, and the measurement does not
   discriminate between them. It argues that a 2 eV shift is not physical at accessible fluence;
   it argues equally against 1.19 eV.
2. **The response is strongly non-monotonic in amplitude.** At `omega = pi`, going 1.5 -> 2.0 ->
   3.0 gives 0.244 -> 0.254 -> 0.204 eV, then 5.0 jumps to 1.188. "Turn the pump down a bit" is
   not a smooth knob here. Any amplitude change needs its deposited energy re-measured, which is
   what `TDOCCSCAN = True` and the `E_stat` column are for.

## Caveat

Matching deposited energy per hole in a one-site cluster is a model-level statement, not a rigorous
fluence match to a sample: there is no volume, no screening environment, and no relaxation channel
here. Use it to put the simulation in the right regime, not to claim a quantitative fluence.
