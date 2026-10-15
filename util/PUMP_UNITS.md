# From PUMPAMPX to a field and a fluence

Derivation behind `util/pump_units.py`. Every formula here was checked numerically against
the code's own waveform; the check values are quoted inline.

## 1. The Peierls substitution

A uniform vector potential couples to a tight-binding model by multiplying each hopping by a
phase (minimal coupling, `p -> p - qA`):

    t_ij  ->  t_ij * exp[ i (q/hbar) A . (r_i - r_j) ]

Equivalently it shifts the crystal momentum, `k -> k - (e/hbar) A` for electrons. Writing
**A** as uniform over the cluster is the dipole approximation, and it is very safe here: the
pump wavelength at 0.25 eV is 5 um against a 2.09 A bond, a ratio of 2.4e4.

## 2. Which length enters — the bond, not the lattice constant

`cluster.cpp` sets `kx = PI/2 + pump_shift(t,0)` and then builds `sx = 2*sin(kx)`. That
structure factor identifies the length. For a TM site with two ligands at +-d x_hat, the
sigma-oriented d-p hybridization in momentum space is

    V(k) = t_pd ( e^{i k_x d} - e^{-i k_x d} ) = 2i t_pd sin(k_x d)

— the relative minus sign is the odd parity of p_x through the TM site — which matches the
code's `2*sin(kx)` with

    kx  ==  k_x * d ,      d = TM-ligand bond length

and `PI/2` is the hole-language zone corner `k_x = pi/a = pi/(2d)`. For NiO, a = 4.177 A so
**d = a/2 = 2.0885 A**. Using the lattice constant here instead of the bond costs a factor 2
in every field and a factor 4 in every fluence.

So `PUMPAMPX` is the amplitude of the dimensionless Peierls phase

    phi(t) = (e d / hbar) A_phys(t) = A0 * exp(-(t-t0)^2 / 2 sigma^2) * sin(omega (t-t0))

in radians, and `A_phys = hbar * phi / (e d)`.

## 3. Phase to electric field

In the gauge the Peierls substitution assumes (zero scalar potential), `E = -dA_phys/dt`:

    E = -(hbar / (e d)) dphi/dt

In the code's natural units (hbar = 1, energies in eV, times in eV^-1) the constants cancel.
With `t_SI = hbar_eV * t_code`, so `dphi/dt_SI = (dphi/dt_code)/hbar_eV`, and
`hbar_SI = e * hbar_eV`:

    E = (e hbar_eV / (e d)) * (dphi/dt_code) / hbar_eV = (dphi/dt_code) / d

and since 1 eV per elementary charge is exactly 1 volt,

    **E [V/m] = |dphi/dt| [eV] / d [m]**

## 4. Many-cycle versus sub-cycle — do not use E = omega*A here

If the envelope is slow compared with the carrier, `|dphi/dt| ~ phi_0 * omega`, giving the
familiar pair

    E_peak = phi_0 * omega[eV] / d          phi_0 = d * E_peak / omega[eV]

This pulse is not in that limit: `omega*sigma = 0.25 * sqrt(1.5) = 0.306`, i.e. the envelope
is 0.115 of an optical cycle. Two consequences, both with closed forms.

**The peak phase is far below A0.** Maximising `exp(-t^2/2 sigma^2) sin(omega t)` with
`sin(omega t) ~ omega t` gives the maximum at `t = sigma`, so

    phi_peak ~ A0 * omega * sigma * e^(-1/2)

For A0 = 25.3: 4.70 analytic against **4.63 numeric**, peaking at t = -1.21 eV^-1 against
sigma = 1.22.

**The field is set by the carrier, not the envelope.** At t = 0 the derivative is all
carrier:

    |dphi/dt|_peak ~ A0 * omega

6.325 eV analytic, **6.325 eV numeric** (exact). The ratio to the many-cycle estimate is

    E_peak / (omega * phi_peak) = e^(1/2) / (omega*sigma) = 5.38

against 5.47 numerically. So `E = omega*A` underestimates this pulse's field by a factor 5.5.

## 5. Field to fluence

Fluence is the time-integrated energy flux, and in vacuum `|S| = eps0 c E^2`:

    F = eps0 * c * integral E(t)^2 dt_SI = eps0 * c * hbar_eV * integral E^2 dt_code

For a many-cycle Gaussian `E = E_0 cos(omega t) exp(-t^2/2 sigma_t^2)`, using
`<cos^2> = 1/2` and `integral exp(-t^2/sigma_t^2) dt = sigma_t sqrt(pi)`:

    F = (1/2) eps0 c E_0^2 sigma_t sqrt(pi),      sigma_t = tau_FWHM / 2.3548

inverted for the planning direction,

    E_0 = sqrt( 2F / (eps0 c sigma_t sqrt(pi)) ),     phi_0 = d E_0 / omega[eV]

For the simulated sub-cycle pulse the integral is done numerically instead.

## 6. The chain, A0 = 25.3

| step | value |
|---|---|
| peak Peierls phase | 4.63 rad = 1.47 pi |
| peak \|dphi/dt\| | 6.325 eV |
| peak vector potential | 1.458e-05 V.s/m |
| peak field | 3.029e+10 V/m = 3.03 V/Ang = 303 MV/cm |
| fluence | 2569 J/m^2 = **257 mJ/cm^2** |

and A0 = 5 gives 0.91 rad, 0.60 V/Ang (60 MV/cm), 10 mJ/cm^2.

## 7. Sanity checks that the waveform is admissible

A propagating pulse carries no DC field, i.e. `integral E dt = -[A(+inf) - A(-inf)] = 0`.
Numerically -1.3e-19 V.s/m, zero to machine precision, because `phi` vanishes at both ends.
A Gaussian vector potential with no carrier would fail this and imply a net momentum kick.

## 8. What this conversion does not include

* **In-crystal, not incident.** This is the field driving the bond. At normal incidence the
  transmitted amplitude is `2/(1+n)` and the flux inside is `n eps0 c E^2`, so the incident
  fluence is larger by `1/T = (1+n)^2/(4n)`, about 1.18 for n ~ 2.3 at 5 um. Modest, but the
  sign is always "the experiment needs more".
* **No propagation.** One cluster, no absorption depth or pulse reshaping. At 0.25 eV NiO is
  far below its ~4 eV gap, so linear absorption is weak and the field does penetrate, but
  nothing here models that.
* **No local-field or screening correction** beyond what the cluster Hamiltonian contains.
* **Peierls only.** The code couples **A** to `tpd`/`tpp` and nothing else: no intra-atomic
  dipole coupling, no diamagnetic A^2 term. That is the standard treatment for a
  tight-binding model and is the right leading term for a below-gap drive, but it is an
  approximation.

Treat the numbers as order-of-magnitude with a well-defined leading term.
