#!/usr/bin/env python3
"""Convert the code's PUMPAMPX into an experimental peak field and fluence, and back.

Full derivation with numerical checks: util/PUMP_UNITS.md

PUMPAMPX is NOT a field. It is the amplitude of a dimensionless Peierls phase, in radians:
`cluster.cpp` sets `kx = PI/2 + pump_shift(t,0)`, and that `kx` is k_x*a/2, so the length
that converts phase to vector potential is the TM-ligand bond d, not the lattice constant.

    phi(t) = (e*d/hbar) * A_phys(t) = A0 * exp(-(t-t0)^2/2 sigma^2) * sin(omega (t-t0))
    E(t)   = -dA_phys/dt   ->   E[V/m] = (dphi/dt)[eV] / d[m]
    many-cycle limit: phi_peak = d * E_peak / omega[eV]

Caveats. This is the field driving the bond INSIDE the crystal; the incident fluence is
larger by reflection and refraction (n(5 um) ~ 2.3 for NiO, so a factor of a few). The
fluence uses the free-space Poynting flux c*eps0*E^2. There is one cluster and no
propagation or absorption depth. Treat everything here as order-of-magnitude.

Usage:  pump_units.py [A0 ...]        report the simulated pulses
        pump_units.py --realistic     what a feasible 5 um pump reaches
"""
import sys
import numpy as np

HBAR = 6.582119569e-16      # eV s
EPS0 = 8.8541878128e-12
C = 2.99792458e8
D_NIO = 2.088e-10           # Ni-O bond length (half the 4.177 A cubic lattice constant)


def simulated(A0, omega=0.25, sigma=np.sqrt(1.5), d=D_NIO):
    """Peak phase, peak in-crystal field and fluence of the code's own waveform."""
    t = np.linspace(-8 * sigma, 8 * sigma, 400001)
    phi = A0 * np.exp(-t**2 / (2 * sigma**2)) * np.sin(omega * t)
    E = np.gradient(phi, t) / d                      # V/m
    F = C * EPS0 * np.sum(E**2) * (t[1] - t[0]) * HBAR   # J/m^2
    return abs(phi).max(), abs(E).max(), F / 10.0    # rad, V/m, mJ/cm^2


def realistic(F_mJcm2, tau_fs, omega=0.25, d=D_NIO, n=1.0):
    """Peak field and Peierls phase of a many-cycle Gaussian pulse of this fluence."""
    sig_t = tau_fs * 1e-15 / 2.3548
    E0 = np.sqrt(2 * F_mJcm2 * 10.0 / (n * C * EPS0 * sig_t * np.sqrt(np.pi)))
    return E0, d * E0 / omega


def main(argv):
    om, sig = 0.25, np.sqrt(1.5)
    print(f"pulse as simulated: omega = {om} eV (period {2*np.pi*HBAR/om*1e15:.1f} fs), "
          f"sigma = {sig:.4f} eV^-1 = {sig*HBAR*1e15:.3f} fs")
    print(f"  envelope FWHM {2.3548*sig*HBAR*1e15:.2f} fs = "
          f"{2.3548*sig*om/(2*np.pi):.3f} optical cycles -> SUB-CYCLE, "
          f"bandwidth 1/sigma = {1/sig:.2f} eV = {1/sig/om:.1f}x the centre frequency")
    if "--realistic" in argv:
        print("\nfeasible 5 um pumps (many-cycle):")
        print(f"  {'fluence':>12s} {'FWHM':>7s} {'E_peak':>12s} {'phi_peak':>10s}")
        for F in (0.5, 1, 3, 10, 30):
            for tau in (50, 100):
                E0, p = realistic(F, tau, om)
                print(f"  {F:8.1f} mJ/cm2 {tau:4.0f} fs {E0/1e8:9.2f} MV/cm {p:10.3f} rad")
        return
    _, p1 = realistic(1.0, 50, om)
    for A0 in [float(a) for a in argv[1:]] or [5.0, 25.3]:
        pk, E, F = simulated(A0, om, sig)
        print(f"\nPUMPAMPX = {A0}")
        print(f"  peak Peierls phase          {pk:7.2f} rad = {pk/np.pi:.2f} pi   "
              f"(NOT A0 -- the sin and the Gaussian barely overlap at omega*sigma = {om*sig:.3f})")
        print(f"  peak in-crystal field       {E:.3e} V/m = {E/1e10:.2f} V/Ang = {E/1e8:.0f} MV/cm")
        print(f"  fluence of this waveform    {F:8.1f} mJ/cm2")
        print(f"  a 50 fs pulse reaching the same phase needs {(pk/p1)**2:8.1f} mJ/cm2")


if __name__ == "__main__":
    main(sys.argv)
