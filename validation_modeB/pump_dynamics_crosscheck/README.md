# Pump dynamics: the Hamiltonians and pump are identical; the legacy propagator is first order in dt

## What matches exactly

**The pumped Hamiltonian.** Our `SquarePlanar::get_tmat` and the Fortran's
`xcef - 2*tmatx*sin((kx+Ax)/2) - 2*tmaty*sin((ky+Ay)/2) - 4*tmato*sin(..)*sin(..)` give the same
one-body spectrum to **1e-14 at every pump amplitude**, once the amplitude convention is applied
(the Fortran's `vec_px` = 2x our `pump_shift`, because its Peierls argument is `(k+A)/2` and ours
is `k/2 + A`; hence `PUMPAMPX = 5.0` against its `pump_ampx = 10.0`):

| our pump_shift A | max &#124;eig(ours) - eig(Fortran at 2A)&#124; | same, before the tpp fix |
|---|---|---|
| 0.00 | 8.9e-15 | 4.5e-02 |
| 1.00 | 3.6e-15 | 2.3e-02 |
| 3.00 | 1.2e-14 | 4.5e-02 |
| 5.00 | 1.1e-14 | 1.2e-02 |

(The right-hand column is why the inter-ligand sign fix mattered: before it the two codes were
solving different Hamiltonians under the pump as well as at rest.)

**The pump envelope.** The Fortran's `pump_ampx*sin(2*pi*pump_freq*t)*exp(-0.5*t^2/pump_sigma)`
with `pump_freq = 0.5`, `pump_sigma = 1.5` is angular frequency `pi` and a Gaussian `exp(-t^2/3)`,
i.e. our `PUMPFREQ = pi` and `PUMPSIGMA = sqrt(1.5) = 1.2247448713915890`. Exact.

**The initial state.** Pre-pump occupations agree to **2e-07**.

## What does not: the legacy time propagation is first order in dt

Occupations over the pulse, ours (converged) against the two Fortran runs, which differ only in
`timesteps` and therefore in `dt = 2*pi/timesteps`:

| Fortran run | dt | max &#124;occupation difference&#124; vs ours |
|---|---|---|
| toy (`timesteps = 80`) | 0.078540 | 7.87e-02 |
| production (`timesteps = 320`) | 0.019635 | 1.97e-02 |

**Ratio 4.01 for a 4x finer dt** -- clean first-order convergence, toward our result. Our own
propagation is converged: over a 20x range of `TDDT` (0.078540 -> 0.003927) our occupations move by
only 9.1e-03, and 1.2e-04 over the last 2.5x.

| our TDDT | max &#124;diff&#124; vs our finest | max &#124;diff&#124; vs the toy Fortran |
|---|---|---|
| 0.078540 | 9.06e-03 | 7.35e-02 |
| 0.039270 | 2.26e-03 | 6.83e-02 |
| 0.009817 | 1.20e-04 | 6.66e-02 |
| 0.003927 | -- | 6.65e-02 |

Note our discrepancy against the Fortran *saturates* at 6.6e-02 rather than vanishing as our dt
shrinks: the two schemes differ at O(dt) in their within-step treatment of H(t), so matching dt
does not match trajectories. The reference is the one that moves when its dt changes.

The differences are confined to the pump-driven orbitals -- `dx2` 7.7e-02, `pxx` 7.9e-02,
`pyy` 4.1e-02 -- while every orbital the x-polarized pump does not touch (`dxy`, `dxz`, `dyz`,
`pxy`, `pxz`, `pyx`, `pyz`) agrees to 1e-16 or better, and both codes conserve the hole number.

## Scope of this finding

At the production `dt` the reference's pumped-state error is ~2% in occupation. That is a real
convergence limitation but it is *not* obviously large enough to explain the outstanding ~2.8x
pump-induced intensity decay in its trRIXS (see `../prx_fig8_antistokes/FIG7_RESOLUTION.md`),
which remains open. It is recorded here as a separate, independently measured fact.

## Reproducing

`scan_td_occupation` hardcodes `print_every = 20`, so to land the output on the Fortran's own grid
choose `TDDT = ref_dt/20`. `INPUT_occupation_scan` uses `TDDT = 0.00392699081698724`, giving output
every 0.0785398 over `[-2*pi, 2*pi]`. The Fortran's `Occ_Pulse.dat` contains **three identical
concatenated blocks** (161 rows each), one per degenerate ground state -- they agree to 4e-14, as
they must, since the Peierls pump couples to charge and not to spin.
