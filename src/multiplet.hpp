#include "hilbert.hpp"
#include "diagonalize.hpp"
#include "gaunt.hpp"
#ifndef MULTIPLET
#define MULTIPLET

double calc_U(double* gaunt1, double* gaunt2, const double* SC, int size);
// target selects which Block matrix (ham / ham_static / ham_hyb) terms are filled into --
// see Hilbert::HTarget. Every function below defaults to HTarget::HAM, i.e. the existing
// equilibrium Hamiltonian, so every pre-existing call site is unaffected. calc_ham(...,
// nohyb=true, target=HAM_STATIC) is how the real-time propagator (src/tdynamics.hpp)
// builds the time-independent piece of H(t) into a block's ham_static, separately from
// (and without disturbing) the ham used to find the equilibrium ground state.
void calc_ham(Hilbert& hilbs, const HParam& hparam, bool nohyb = false,
				Hilbert::HTarget target = Hilbert::HTarget::HAM);
void calc_coulomb(Hilbert& hilbs, const std::vector<double*>& SC,
					Hilbert::HTarget target = Hilbert::HTarget::HAM);
vecd CFmat(int l, const double* CF);
void calc_CF(Hilbert& hilbs, const double* CF, Hilbert::HTarget target = Hilbert::HTarget::HAM);
void calc_SO(Hilbert& hilbs, const double lambda, int l_in,
				Hilbert::HTarget target = Hilbert::HTarget::HAM);
void calc_CV(Hilbert& hilbs, const double* FG, Hilbert::HTarget target = Hilbert::HTarget::HAM);
void calc_MF(Hilbert& hilbs, const double h_field, Hilbert::HTarget target = Hilbert::HTarget::HAM);
// t: real-time propagation time (natural units, eV^-1), forwarded to Cluster::get_tmat_real(t).
// target: HAM for the (time-independent) equilibrium hybridization, HAM_HYB when rebuilding
// only the pump-modulated piece at time t during real-time propagation.
void calc_HYB(Hilbert& hilbs, const HParam& hparam, bool nohyb, double t = 0.0,
				Hilbert::HTarget target = Hilbert::HTarget::HAM);


#endif