#include "tdynamics.hpp"
#include "photon.hpp"
#include <cstdlib>
#ifndef TR_PHOTON
#define TR_PHOTON

// Time-resolved XAS/RIXS (Phase 3-4). Reuses the equilibrium spectral solvers
// (ContFracExpan/BiCGS, from diagonalize.hpp) and dipole machinery (gen_dipole_state/
// basis_overlap, from photon.cpp) essentially unchanged: at each probe-window sample time
// t1, the (possibly pump-modulated) core-hole Hamiltonian is treated as *frozen* at its
// instantaneous value -- valid whenever the core-hole lifetime (HParam::tau_ch) is short
// compared to the pump width (HParam::pump_sigma), which holds for realistic femtosecond
// pump pulses and sub-femtosecond core-hole lifetimes. This avoids ever needing an explicit
// real-time sweep over the core-hole dwell time, which is what makes this tractable (see
// tau_ch's documentation in atoms.hpp for the cost this avoids).
//
// All times are natural units (hbar=1), eV^-1, matching tdynamics.hpp.

inline double gaussian_weight(double t, double t0, double sigma) {
	double dt = t - t0;
	return std::exp(-dt*dt/(2*sigma*sigma));
}

// Rebuilds ham_hyb across every block of hilbs at the given instant t (no propagation --
// a pure snapshot), for use as a frozen instantaneous Hamiltonian. Mirrors the clear+refill
// sequence in propagate_block_state/validate_propagator exactly.
inline void snapshot_ham_hyb(Hilbert& hilbs, const HParam& hparam, double t) {
	for (auto& blk : hilbs.hblks) clear_ham_hyb(blk);
	if (hilbs.HYB_on) calc_HYB(hilbs, hparam, /*nohyb=*/false, t, Hilbert::HTarget::HAM_HYB);
	return;
}

// Instantaneous energy expectation <psi|H(t)|psi> under the (possibly pump-modulated)
// snapshot Hamiltonian at time t. Needed as the CFE/BiCGS energy anchor in place of the
// static equilibrium gs_en once the pump has acted: psi(t1) is then a genuine superposition
// with real amplitude on excited configurations, sitting at a HIGHER energy than the static
// ground state -- anchoring to the static gs_en instead systematically hides any resulting
// anti-Stokes (energy-gain) RIXS feature, since it under-estimates how far "up" the initial
// state actually sits.
inline double instantaneous_energy(Hilbert& hilbs, const HParam& hparam, size_t blk_ind,
									double t, const vecc& psi) {
	snapshot_ham_hyb(hilbs, hparam, t);
	CombinedHam<double> Ht(hilbs.hblks[blk_ind].ham_static, hilbs.hblks[blk_ind].ham_hyb);
	vecc Hpsi(psi.size(),0);
	Ht.mvmult_cmplx(psi, Hpsi);
	return ed::inner_product(psi,Hpsi).real();
}

// Propagates psi0 forward from t_lo to t_hi at step dt, recording the state at every grid
// point (index i <-> t = t_lo + i*dt). For the small test systems this code targets, caching
// the dense trajectory is simplest and cheapest to get right; a rolling/windowed cache would
// be needed to scale to much larger Hilbert spaces (not needed here).
inline std::vector<vecc> build_trajectory(Hilbert& hilbs, size_t blk_ind, const vecc& psi0,
											double t_lo, double t_hi, double dt, int m,
											const HParam& hparam) {
	int nsteps = (int)std::round((t_hi-t_lo)/dt);
	std::vector<vecc> traj(nsteps+1);
	vecc psi(psi0);
	traj[0] = psi;
	double t = t_lo;
	for (int s = 0; s < nsteps; ++s) {
		propagate_block_state(hilbs, blk_ind, psi, t, dt, hparam, m);
		t += dt;
		traj[s+1] = psi;
	}
	return traj;
}

inline void write_trXAS(const vecd& aben, const std::vector<vecd>& intensity,
							const vecd& delays, std::string prefix) {
	std::ofstream idx(prefix + "_delays.txt");
	idx << std::setw(15) << "delay (eV^-1)" << "  file" << std::endl;
	for (size_t d = 0; d < delays.size(); ++d) {
		std::string fname = prefix + "_t" + std::to_string(delays[d]) + ".txt";
		idx << std::setw(15) << delays[d] << "  " << fname << std::endl;
		std::ofstream fout(fname);
		fout << std::setw(15) << "peaks (eV)" << std::setw(15) << "intensity" << std::endl;
		for (size_t i = 0; i < aben.size(); ++i)
			fout << std::setw(15) << aben[i] << std::setw(15) << intensity[d][i] << std::endl;
	}
	std::cout << "Wrote " << delays.size() << " trXAS delay files, index: "
				<< prefix + "_delays.txt" << std::endl;
	return;
}

// Time-resolved XAS. hparam_gs/hparam_ex must be correctly-scoped copies of the run's
// HParam (ground-state vs. core-hole Coulomb/SOC values respectively -- see main.cpp, where
// process_hilbert_space leaves hparam mutated to core-hole values after building EX and
// never restores them). GS/EX must already be built+diagonalized via the ordinary
// equilibrium pipeline (calc_ham+Block::diagonalize), exactly as for validate_propagator.
inline void trXAS(Hilbert& GS, Hilbert& EX, const HParam& hparam_gs, const HParam& hparam_ex,
					const PM& pm, TRParam& trp) {
	std::cout << std::endl << "=== Time-resolved XAS ===" << std::endl;
	trp.set_delay_points();
	if (trp.delay_points.empty()) {std::cout << "No TRDELAY points set" << std::endl; return;}
	double dt_probe = (trp.dt_probe > 0) ? trp.dt_probe : hparam_gs.td_dt;
	double t_lo = *std::min_element(trp.delay_points.begin(),trp.delay_points.end())
					- trp.n_sigma_probe*trp.probe_sigma;
	double t_hi = *std::max_element(trp.delay_points.begin(),trp.delay_points.end())
					+ trp.n_sigma_probe*trp.probe_sigma;
	int n_is = (int)std::round((t_hi-t_lo)/dt_probe);

	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.nev; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	std::vector<bindex> gsi;
	for (size_t i = 0; i < GS.hblks.size(); ++i)
		for (size_t j = 0; j < GS.hblks[i].nev; ++j)
			if (std::abs(GS.hblks[i].eig[j]-gs_en) < TOL) gsi.push_back({i,j});
	std::cout << "Ground state(s): " << gsi.size() << ", energy: " << gs_en << std::endl;

	init_time_dependent_ham(GS, hparam_gs);
	init_time_dependent_ham(EX, hparam_ex);

	double nedos = pm.nedos;
	vecd xas_aben(nedos,0);
	for (int i = 0; i < nedos; ++i)
		xas_aben[i] = pm.ab_range[0] + (pm.ab_range[1]-pm.ab_range[0])/nedos*i;
	double Gamma_ch = 1.0/hparam_gs.tau_ch;

	// Total probe weight per delay, for normalizing the weighted average (so the pump-off,
	// wide-probe equilibrium limit recovers the same overall magnitude as static XAS()).
	vecd norm_weight(trp.delay_points.size(),0);
	for (int is = 0; is <= n_is; ++is) {
		double t1 = t_lo + is*dt_probe;
		for (size_t d = 0; d < trp.delay_points.size(); ++d)
			norm_weight[d] += gaussian_weight(t1, trp.delay_points[d], trp.probe_sigma);
	}

	std::vector<vecd> xas_int(trp.delay_points.size(), vecd(nedos,0));
	std::vector<blapIndex> blap;

	for (auto& g : gsi) {
		size_t N = GS.hblks[g.first].size;
		vecc psi0(N,0);
		for (size_t i = 0; i < N; ++i) psi0[i] = dcomp(GS.hblks[g.first].eigvec[g.second*N+i],0);
		std::cout << "Building ground-state trajectory over [" << t_lo << "," << t_hi << "]..." << std::endl;
		auto traj = build_trajectory(GS, g.first, psi0, t_lo, t_hi, dt_probe, hparam_gs.td_m, hparam_gs);

		for (auto& exblk : EX.hblks) {
			size_t exblk_ind = &exblk - &EX.hblks[0];
			if (!GS.SO_on && !EX.SO_on && GS.hblks[g.first].get_sz() != exblk.get_sz()) continue;
			basis_overlap(GS, EX, bindex(g.first,exblk_ind), blap, pm);
			if (blap.empty()) continue;

			for (int is = 0; is <= n_is; ++is) {
				double t1 = t_lo + is*dt_probe;
				vecc dipole_vec = gen_dipole_state(GS, EX, pm, bindex(g.first,exblk_ind), traj[is], blap);
				if (ed::norm(dipole_vec) < TOL) continue;

				double E_gs_t1 = instantaneous_energy(GS, hparam_gs, g.first, t1, traj[is]);
				snapshot_ham_hyb(EX, hparam_ex, t1);
				CombinedHam<double> Ht_EX(exblk.ham_static, exblk.ham_hyb);

				vecd shape(nedos,0);
				ContFracExpan(&Ht_EX, dipole_vec, E_gs_t1, xas_aben, shape, Gamma_ch, pm.niterCFE);

				for (size_t d = 0; d < trp.delay_points.size(); ++d) {
					double w = gaussian_weight(t1, trp.delay_points[d], trp.probe_sigma);
					if (w < 1e-6) continue;
					for (int i = 0; i < nedos; ++i) xas_int[d][i] += w*shape[i];
				}
			}
		}
	}

	for (size_t d = 0; d < trp.delay_points.size(); ++d) {
		double denom = (gsi.size()>0 ? gsi.size() : 1) * (norm_weight[d] > TOL ? norm_weight[d] : 1);
		for (auto& v : xas_int[d]) v /= denom;
	}
	write_trXAS(xas_aben, xas_int, trp.delay_points, "trXAS_"+pm.edge+"edge_"+pol_str(pm.pvin));
	std::cout << "=== End trXAS ===" << std::endl << std::endl;
	return;
}

inline void write_trRIXS(const vecd& ab, const vecd& em, const std::vector<vecd>& peaks_per_delay,
							const vecd& delays, std::string prefix) {
	std::ofstream idx(prefix + "_delays.txt");
	idx << std::setw(15) << "delay (eV^-1)" << "  file" << std::endl;
	for (size_t d = 0; d < delays.size(); ++d) {
		std::string fname = prefix + "_t" + std::to_string(delays[d]) + ".txt";
		idx << std::setw(15) << delays[d] << "  " << fname << std::endl;
		std::ofstream fout(fname);
		fout << std::setw(18) << "absorption (eV)" << std::setw(18) << "energy loss (eV)"
				<< std::setw(18) << "intensity" << std::endl;
		for (size_t y = 0; y < em.size(); ++y) {
			for (size_t x = 0; x < ab.size(); ++x) {
				double val = peaks_per_delay[d][y*ab.size()+x];
				if (val < TOL) continue;
				fout << std::setw(18) << ab[x] << std::setw(18) << em[y]
						<< std::setw(18) << val << std::endl;
			}
		}
	}
	std::cout << "Wrote " << delays.size() << " trRIXS delay files, index: "
				<< prefix + "_delays.txt" << std::endl;
	return;
}

// Time-resolved RIXS, "Mode A" (frozen-instantaneous-Hamiltonian): reuses the existing
// equilibrium RIXS() solver=4 inner loop (BiCGS for the incident-energy/core-hole-dwell
// axis, ContFracExpan for the loss-energy/final-state axis) essentially unchanged, fed a
// pump-snapshotted frozen core-hole Hamiltonian and the actual propagated ground state at
// each probe-window sample, instead of the static equilibrium values -- valid whenever
// tau_ch << pump_sigma (see atoms.hpp's tau_ch documentation). The loss-energy axis uses
// the *static* equilibrium GS Hamiltonian (GS.hblks[*].ham) rather than a pump-aware
// snapshot -- a documented simplification (the plan's "hybrid_s_handoff" for probe samples
// landing directly on the pump is not implemented here; most probe samples in a realistic
// scan are far enough from the pump center that this makes no difference).
inline void trRIXS(Hilbert& GS, Hilbert& EX, const HParam& hparam_gs, const HParam& hparam_ex,
					const PM& pm, TRParam& trp) {
	std::cout << std::endl << "=== Time-resolved RIXS ===" << std::endl;
	trp.set_delay_points();
	if (trp.delay_points.empty()) {std::cout << "No TRDELAY points set" << std::endl; return;}
	if (pm.inc_e_points.empty()) {std::cout << "No INCIDENT points set" << std::endl; return;}
	double dt_probe = (trp.dt_probe > 0) ? trp.dt_probe : hparam_gs.td_dt;
	double t_lo = *std::min_element(trp.delay_points.begin(),trp.delay_points.end())
					- trp.n_sigma_probe*trp.probe_sigma;
	double t_hi = *std::max_element(trp.delay_points.begin(),trp.delay_points.end())
					+ trp.n_sigma_probe*trp.probe_sigma;
	int n_is = (int)std::round((t_hi-t_lo)/dt_probe);

	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.nev; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	std::vector<bindex> gsi;
	for (size_t i = 0; i < GS.hblks.size(); ++i)
		for (size_t j = 0; j < GS.hblks[i].nev; ++j)
			if (std::abs(GS.hblks[i].eig[j]-gs_en) < TOL) gsi.push_back({i,j});
	std::cout << "Ground state(s): " << gsi.size() << ", energy: " << gs_en << std::endl;

	init_time_dependent_ham(GS, hparam_gs);
	init_time_dependent_ham(EX, hparam_ex);

	double nedos = pm.nedos;
	double Gamma_ch = 1.0/hparam_gs.tau_ch;
	size_t n_ab = pm.inc_e_points.size();
	vecd rixs_em(nedos,0);
	for (int i = 0; i < nedos; ++i) rixs_em[i] = pm.em_lo + (pm.em_energy-pm.em_lo)/nedos*i;

	vecd norm_weight(trp.delay_points.size(),0);
	for (int is = 0; is <= n_is; ++is) {
		double t1 = t_lo + is*dt_probe;
		for (size_t d = 0; d < trp.delay_points.size(); ++d)
			norm_weight[d] += gaussian_weight(t1, trp.delay_points[d], trp.probe_sigma);
	}

	std::vector<vecd> rixs_peaks(trp.delay_points.size(), vecd(nedos*n_ab,0));
	std::vector<blapIndex> blap_in, blap_out;

	for (auto& g : gsi) {
		size_t N = GS.hblks[g.first].size;
		vecc psi0(N,0);
		for (size_t i = 0; i < N; ++i) psi0[i] = dcomp(GS.hblks[g.first].eigvec[g.second*N+i],0);
		std::cout << "Building ground-state trajectory over [" << t_lo << "," << t_hi << "]..." << std::endl;
		auto traj = build_trajectory(GS, g.first, psi0, t_lo, t_hi, dt_probe, hparam_gs.td_m, hparam_gs);

		for (auto& exblk : EX.hblks) {
			size_t exblk_ind = &exblk - &EX.hblks[0];
			if (!GS.SO_on && !EX.SO_on && GS.hblks[g.first].get_sz() != exblk.get_sz()) continue;
			basis_overlap(GS, EX, bindex(g.first,exblk_ind), blap_in, pm);
			if (blap_in.empty()) continue;
			if (pm.pvin != pm.pvout) basis_overlap(GS, EX, bindex(g.first,exblk_ind), blap_out, pm, true);
			else blap_out = blap_in;

			for (int is = 0; is <= n_is; ++is) {
				double t1 = t_lo + is*dt_probe;
				bool any_delay_near = false;
				for (size_t d = 0; d < trp.delay_points.size(); ++d)
					if (gaussian_weight(t1,trp.delay_points[d],trp.probe_sigma) >= 1e-6) {any_delay_near = true; break;}
				if (!any_delay_near) continue;

				vecc dipole_vec = gen_dipole_state(GS, EX, pm, bindex(g.first,exblk_ind), traj[is], blap_in);
				if (ed::norm(dipole_vec) < TOL) continue;

				double E_gs_t1 = instantaneous_energy(GS, hparam_gs, g.first, t1, traj[is]);
				snapshot_ham_hyb(EX, hparam_ex, t1);
				CombinedHam<double> Ht_EX(exblk.ham_static, exblk.ham_hyb);

				vecd shape(nedos*n_ab,0);
				for (size_t ix = 0; ix < n_ab; ++ix) {
					double ab_en = pm.inc_e_points[ix];
					dcomp z(E_gs_t1+ab_en, -Gamma_ch);
					vecc midvec = BiCGS(&Ht_EX, dipole_vec, z, pm.CG_tol);
					midvec = gen_dipole_state(GS, EX, pm, bindex(g.first,exblk_ind), midvec, blap_out, false);
					if (ed::norm(midvec) < TOL) continue;
					vecd col(nedos,0);
					ContFracExpan(GS.hblks[g.first].ham, midvec, E_gs_t1, rixs_em, col, pm.eps_loss, pm.niterCFE);
					for (int iy = 0; iy < nedos; ++iy) shape[iy*n_ab+ix] = col[iy];
				}

				for (size_t d = 0; d < trp.delay_points.size(); ++d) {
					double w = gaussian_weight(t1, trp.delay_points[d], trp.probe_sigma);
					if (w < 1e-6) continue;
					for (size_t k = 0; k < shape.size(); ++k) rixs_peaks[d][k] += w*shape[k];
				}
			}
		}
	}

	// NOTE: unlike XAS() (which unconditionally divides by gsi.size() for every solver),
	// equilibrium RIXS()'s solver=4 path does NOT divide by gsi.size() at all -- that
	// normalization only appears in the solver=1/2/3 branches. Match solver=4's actual
	// behavior here (probe-weight normalization only) so the equilibrium limit agrees
	// exactly with the existing RIXS_*.txt solver=4 output.
	for (size_t d = 0; d < trp.delay_points.size(); ++d) {
		double denom = (norm_weight[d] > TOL ? norm_weight[d] : 1);
		for (auto& v : rixs_peaks[d]) v /= denom;
	}
	write_trRIXS(pm.inc_e_points, rixs_em, rixs_peaks, trp.delay_points,
					"trRIXS_"+pm.edge+"edge_"+pol_str(pm.pvin)+"_"+pol_str(pm.pvout));
	std::cout << "=== End trRIXS ===" << std::endl << std::endl;
	return;
}

// ============================================================================
// "Mode B": genuine two-time correlator trRIXS, matching the reference Fortran's actual
// time-domain Kramers-Heisenberg formula (not Mode A's frozen-instantaneous-snapshot
// approximation). Mode A incoherently averages already-squared per-snapshot intensities
// across the probe window, which destroys the cross-time interference between different
// absorption times that produces correctly-resolved spectral features when that
// interference matters (i.e. when tau_ch is not much shorter than the internal dynamics
// timescale). Mode B instead accumulates a coherent complex amplitude over a causally-
// ordered, tau_ch-bounded (t1,t2) window and squares only once at the end:
//
//   csum(win,wloss,T) = sum_{t1<=t2<=t1+W} g(t1-T) g(t2-T) exp(-(t2-t1)/tau_ch)
//                          K(t1,t2) exp(-i*win*t1 + i*(win-wloss)*t2) * dt^2
//   I(win,wloss,T) = |csum|^2
//
// where K(t1,t2) = <psi_GS(t2)| D_out  U_core(t2,t1)  D_in |psi_GS(t1)>, computed by
// explicitly propagating the dipole-excited state through EX's own (pump-modulated)
// Hamiltonian from t1 to t2 -- reusing propagate_block_state exactly as-is, since it
// already takes an arbitrary Hilbert& (GS or EX). g(x) is a normalized Gaussian probe
// window, 1/(probe_sigma*sqrt(2*pi)) * exp(-x^2/(2*probe_sigma^2)), matching the
// reference's own "gauss*exp(...)/sqrt(2*pi)" normalization exactly.
//
// Efficiency vs. the reference: (1) K(t1,t2) does not depend on the delay T, so it's
// built ONCE and reused across every requested delay (the reference recomputes it once
// per delay, redundantly); (2) the t2 sum is truncated at a tau_ch-bounded window (the
// reference sums to the end of the simulated range even though the exp(-(t2-t1)/tau_ch)
// factor makes distant terms negligible); (3) U_core(t2,t1) is applied via short-
// iterative-Lanczos propagation of a single state vector, not full dense diagonalization
// plus an ndim x ndim matrix composition at every timestep.
//
// Multiple matching core-hole blocks (ex_blk_inds) for the same ground state are accumulated
// COHERENTLY into one shared Krows amplitude array before the final |csum|^2 -- summing per-block
// *intensities* instead (as an earlier version of this function did, one independent call and
// output file per exblk) silently discards interference between amplitude contributions that
// reach the same final valence state via different intermediate core-hole blocks. Confirmed against
// the reference Fortran: its entire core-hole manifold is one undifferentiated space, and its own
// vertex sum (the "i" loop inside csum, before tkern(t1,t2)=csum) sums the amplitude over every
// core-hole eigenstate coherently, squaring only once, much later, after the full (t1,t2) sum --
// never per-block. Summing multiple degenerate ground states coherently too (not yet implemented)
// would require the same treatment if that case is ever exercised.
// ============================================================================

inline void write_trRIXS_2time(const vecd& ab, const vecd& em, const std::vector<vecd>& peaks_per_delay,
								const vecd& delays, std::string prefix) {
	std::ofstream idx(prefix + "_delays.txt");
	idx << std::setw(15) << "delay (eV^-1)" << "  file" << std::endl;
	for (size_t d = 0; d < delays.size(); ++d) {
		std::string fname = prefix + "_t" + std::to_string(delays[d]) + ".txt";
		idx << std::setw(15) << delays[d] << "  " << fname << std::endl;
		std::ofstream fout(fname);
		fout << std::setw(18) << "absorption (eV)" << std::setw(18) << "energy loss (eV)"
				<< std::setw(18) << "intensity" << std::endl;
		for (size_t y = 0; y < em.size(); ++y) {
			for (size_t x = 0; x < ab.size(); ++x) {
				fout << std::setw(18) << ab[x] << std::setw(18) << em[y]
						<< std::setw(18) << peaks_per_delay[d][y*ab.size()+x] << std::endl;
			}
		}
	}
	std::cout << "Wrote " << delays.size() << " trRIXS(2-time) delay files, index: "
				<< prefix + "_delays.txt" << std::endl;
	return;
}

inline void trRIXS_2time(Hilbert& GS, Hilbert& EX, const HParam& hparam_gs, const HParam& hparam_ex,
							const PM& pm, size_t gs_blk_ind, size_t gs_eig_ind,
							const std::vector<size_t>& ex_blk_inds,
							const vecd& ab_grid, const vecd& loss_grid, const vecd& delays,
							double probe_sigma, double tau_ch, double dt, int m,
							double sim_window_mult = 4.0, double t2_max_window = -1,
							std::string prefix = "trRIXS2",
							bool coherent_final_states = false) {
	std::cout << std::endl << "=== Time-resolved RIXS (two-time correlator, Mode B) ===" << std::endl;

	init_time_dependent_ham(GS, hparam_gs);
	init_time_dependent_ham(EX, hparam_ex);

	Block<double>& gsblk = GS.hblks[gs_blk_ind];
	size_t N = gsblk.size;
	vecc psi0(N,0);
	for (size_t i = 0; i < N; ++i) psi0[i] = dcomp(gsblk.eigvec[gs_eig_ind*N+i],0);

	// sim_window_mult extends [t_lo,t_hi] beyond the requested delays to capture the probe's
	// own Gaussian tails -- independent of t2_max_window (below), which bounds how far t2 is
	// allowed to run past t1 for the tau_ch decay. Conflating these two was a real bug in an
	// earlier version of this function: for validating directly against a reference run whose
	// own (t1,t2) sum is NOT tau_ch-truncated, t2_max_window must default to "no truncation"
	// (covering the full remaining simulated range), independent of how tight sim_window_mult
	// is set for matching that reference's own bounded time window.
	double window = sim_window_mult * std::max(probe_sigma, tau_ch);
	double t_lo = *std::min_element(delays.begin(),delays.end()) - window;
	double t_hi = *std::max_element(delays.begin(),delays.end()) + window;
	int n_is = (int)std::round((t_hi-t_lo)/dt);
	std::cout << "Building ground-state trajectory over [" << t_lo << "," << t_hi << "]..." << std::endl;
	auto gs_traj = build_trajectory(GS, gs_blk_ind, psi0, t_lo, t_hi, dt, m, hparam_gs);

	// De-excitation side. PRB 99, 104306 Eq. (11) contains, between the two de-excitation
	// dipole operators, a propagator U(t2',t2) acting in the FINAL (core-hole-free) valence
	// sector. Resolving it with a complete valence basis {|l>} taken at the reference time t_lo,
	//     U(t2',t2) = sum_l U(t2',t_lo)|l><l|U(t2,t_lo)^dag ,
	// factorizes the four-time correlator of Eq. (10) as
	//     S(t1',t2',t2,t1) = sum_l conj(a_l(t1',t2')) a_l(t1,t2),
	//     a_l(t1,t2) = <val_l(t2)| D_out U_core(t2,t1) D_in |psi_GS(t1)>,   |val_l(t)> = U(t,t_lo)|l>,
	// so that the (t1,t2,t1',t2') quadruple integral collapses to a sum of moduli squared,
	//     I = sum_l |A_l|^2,   A_l = sum_{t1<t2} (weights) a_l(t1,t2) e^{-i*win*t1 + i*wout*t2},
	// i.e. INCOHERENT over the final valence state l, coherent only over (t1,t2). This is the
	// nonequilibrium form of the paper's own equilibrium limit, their Eq. (15),
	// I_eq = sum_f |A_f|^2 delta(w - E_f + E_i) -- a sum of squared amplitudes over final
	// states, not the square of their sum.
	//
	// An earlier version of this function summed a_l over l into one scalar kernel and squared
	// once (|sum_l A_l|^2), mirroring the paper-1 reference Fortran, whose own kernel loop
	// ("do l=1,min(nev,binstates(wsp)) ... csum=csum+dconjg(tfvert(wsp,l,i,t2))*tpivert(...)",
	// then a single dreal(csum4*dconjg(csum4)) after the frequency sum) does the same. That is
	// not an observable: |sum_l <l|M>|^2 = |<u|M>|^2 with |u> = sum_l |l> is a projection onto a
	// direction fixed entirely by the arbitrary sign/phase LAPACK happens to assign each
	// eigenvector, whereas sum_l |<l|M>|^2 = <M|M> is basis independent. Set the
	// TRRIXS2_FLIP_VAL environment variable to a valence index to exercise exactly that
	// invariance (see below). Measured on the toy d8 validation case: flipping one valence
	// eigenvector moves the coherent path's intensity by 57-172%, and the corrected path by
	// <=0.4% -- the latter being this code's own run-to-run numerical noise (two identical runs
	// differ by up to ~1e-4 here), not a real dependence.
	size_t nval = gsblk.nev;
	if (nval < N) std::cerr << "WARNING: trRIXS_2time needs a COMPLETE valence eigenbasis for the "
							<< "final-state sum (sum_l |l><l| = 1); block dim " << N << " but only "
							<< nval << " eigenvectors available -- raise GSNEV." << std::endl;
	// Diagnostic only: flip the sign of one valence eigenvector. A correct cross section cannot
	// depend on this; the coherent-final-state path visibly does.
	int flip_val = -1;
	if (const char* fv = std::getenv("TRRIXS2_FLIP_VAL")) {
		flip_val = std::atoi(fv);
		std::cout << "DIAGNOSTIC: sign-flipping valence eigenvector " << flip_val << std::endl;
	}
	// Vmat[is] holds every |val_l(t_is)> contiguously (component i of state l at [l*N+i]), so the
	// K-row loop below can project one de-excited vector onto the whole final-state basis at once.
	std::cout << "Building " << nval << " valence-eigenstate trajectories (de-excitation basis)..." << std::endl;
	std::vector<vecc> Vmat(n_is+1, vecc(nval*N, dcomp(0,0)));
	for (size_t l = 0; l < nval; ++l) {
		vecc l0(N,0);
		double sgn = ((int)l == flip_val) ? -1.0 : 1.0;
		for (size_t i = 0; i < N; ++i) l0[i] = dcomp(sgn*gsblk.eigvec[l*N+i],0);
		auto lt = build_trajectory(GS, gs_blk_ind, l0, t_lo, t_hi, dt, m, hparam_gs);
		for (int is = 0; is <= n_is; ++is)
			for (size_t i = 0; i < N; ++i) Vmat[is][l*N+i] = lt[is][i];
		if (l % 20 == 0) std::cout << "  valence trajectory " << l << "/" << nval << std::endl;
	}

	int W_steps = (t2_max_window > 0) ? (int)std::round(t2_max_window/dt) : n_is;
	std::cout << "t2 window: " << W_steps << " steps ("
				<< (t2_max_window > 0 ? t2_max_window : (t_hi-t_lo)) << " eV^-1)" << std::endl;

	// a_l(t1,t2), causally ordered STRICTLY t2>t1 (matching the reference's own t2 loop, which
	// starts at t1+1 and never includes the equal-time t2=t1 term -- that term carries no
	// energy resolution at all and Fourier-transforms into a spurious peak at zero energy
	// loss if included), tau_ch-window-truncated, and built ONCE, reused for every requested
	// delay (the reference rebuilds this once per delay, redundantly).
	// Krows[is1][j*nval + l] holds a_l(t1, t1+(j+1)*dt): one amplitude per FINAL valence state l
	// (kept separate, per Eq. (11) -- see the derivation above), summed coherently over every
	// matching core-hole block in ex_blk_inds, since those are genuinely different intermediate
	// pathways to the same final state and interfere.
	std::vector<std::vector<dcomp>> Krows(n_is+1);
	for (int is1 = 0; is1 <= n_is; ++is1)
		Krows[is1].assign((size_t)std::min(W_steps, n_is-is1)*nval, dcomp(0,0));

	for (size_t ex_blk_ind : ex_blk_inds) {
		std::vector<blapIndex> blap_in, blap_out;
		basis_overlap(GS, EX, bindex(gs_blk_ind,ex_blk_ind), blap_in, pm);
		if (pm.pvin != pm.pvout) basis_overlap(GS, EX, bindex(gs_blk_ind,ex_blk_ind), blap_out, pm, true);
		else blap_out = blap_in;

		std::cout << "  accumulating core-hole block " << ex_blk_ind << "..." << std::endl;
		for (int is1 = 0; is1 <= n_is; ++is1) {
			int kmax = (int)(Krows[is1].size()/nval);
			vecc xi = gen_dipole_state(GS, EX, pm, bindex(gs_blk_ind,ex_blk_ind), gs_traj[is1], blap_in, true);
			if (ed::norm(xi) < TOL) continue;
			double t = t_lo + is1*dt;
			for (int k = 1; k <= kmax; ++k) {
				propagate_block_state(EX, ex_blk_ind, xi, t, dt, hparam_ex, m);
				t += dt;
				vecc deexc = gen_dipole_state(GS, EX, pm, bindex(gs_blk_ind,ex_blk_ind), xi, blap_out, false);
				const vecc& V = Vmat[is1+k];
				dcomp* dst = &Krows[is1][(size_t)(k-1)*nval];
				for (size_t l = 0; l < nval; ++l) {
					dcomp acc(0,0);
					const dcomp* vl = &V[l*N];
					for (size_t i = 0; i < N; ++i) acc += std::conj(vl[i]) * deexc[i];
					dst[l] += acc;
				}
			}
			if (is1 % 20 == 0) std::cout << "  K-row " << is1 << "/" << n_is << std::endl;
		}
	}

	// DEBUG: dump sum_l a_l(t1=0, t2=k*dt) for k=1..5 -- the l-summed combination is what the
	// reference stores as tkern(t1,t2), so this stays directly comparable to its own debug dump
	// even though the l components are now kept separate for the (correct) incoherent final sum.
	{
		int is0 = (int)std::round((0.0-t_lo)/dt);
		if (is0 >= 0 && is0 <= n_is) {
			for (int k = 1; k <= 5 && (size_t)k*nval <= Krows[is0].size(); ++k) {
				dcomp tk(0,0);
				for (size_t l = 0; l < nval; ++l) tk += Krows[is0][(size_t)(k-1)*nval+l];
				std::cout << "DEBUG Krows t1=0 t2=" << k << " val=" << tk << std::endl;
			}
		}
	}

	size_t n_ab = ab_grid.size(), n_loss = loss_grid.size();
	std::vector<vecd> intensity(delays.size(), vecd(n_loss*n_ab, 0));
	double gnorm = 1.0/(probe_sigma*std::sqrt(2*PI));
	if (coherent_final_states)
		std::cout << "DIAGNOSTIC: summing final valence states COHERENTLY (reproduces the "
					<< "reference Fortran's convention; not the Eq. (11) cross section)." << std::endl;

	// The (t1,t2) sum factorizes: with t2 = t1 + (j+1)*dt, the phase is
	//     e^{-i*win*t1 + i*wout*t2} = e^{+i*win*(j+1)*dt} * e^{-i*wloss*t2},
	// so the win-dependent part depends only on the OFFSET j and the wloss-dependent part only
	// on t2. That lets the whole (win, wloss) grid be filled from ONE pass over the (t1,t2)
	// pairs per delay, instead of redoing the double sum for every (win, wloss).
	//
	// The loop is ordered t2-outer so that everything hot stays in cache: for a fixed t2, the
	// partial sums Cblk[win][l] are only n_ab*nval complex (tens of kB) and the Krows row being
	// read is nval complex, while the full accumulator A[win][wloss][l] is touched once per t2.
	// An earlier, more obvious ordering (win-outer, re-reading all of Krows for every incident
	// energy) is arithmetically identical but ran ~2 orders of magnitude slower at production
	// size -- it streamed the multi-GB Krows array n_ab times and scattered its writes over it.
	int kmax_all = 0;
	for (int is1 = 0; is1 <= n_is; ++is1) kmax_all = std::max(kmax_all, (int)(Krows[is1].size()/nval));
	vecd l12tab(kmax_all);
	for (int j = 0; j < kmax_all; ++j) l12tab[j] = std::exp(-(j+1)*dt/tau_ch);
	if (coherent_final_states)
		std::cout << "DIAGNOSTIC: summing final valence states COHERENTLY (reproduces the "
					<< "reference Fortran's convention; not the Eq. (11) cross section)." << std::endl;

	#pragma omp parallel for schedule(dynamic)
	for (int d = 0; d < (int)delays.size(); ++d) {
		double T = delays[d];
		vecd gtab(n_is+1);
		for (int is = 0; is <= n_is; ++is) {
			double tt = t_lo + is*dt - T;
			gtab[is] = gnorm*std::exp(-tt*tt/(2*probe_sigma*probe_sigma));
		}
		// ph[ix][j] carries the incident phase, the core-hole decay and the dt^2 measure.
		std::vector<dcomp> ph((size_t)n_ab*kmax_all);
		for (size_t ix = 0; ix < n_ab; ++ix)
			for (int j = 0; j < kmax_all; ++j)
				ph[ix*kmax_all+j] = std::exp(dcomp(0,ab_grid[ix]*(j+1)*dt))*(l12tab[j]*dt*dt);
		std::vector<dcomp> A((size_t)n_ab*n_loss*nval, dcomp(0,0));
		std::vector<dcomp> Cblk((size_t)n_ab*nval);

		for (int it2 = 1; it2 <= n_is; ++it2) {
			double g2 = gtab[it2];
			if (g2 < 1e-10) continue;
			std::fill(Cblk.begin(), Cblk.end(), dcomp(0,0));
			bool any = false;
			for (int is1 = 0; is1 < it2; ++is1) {
				int j = it2-1-is1;
				if (j >= (int)(Krows[is1].size()/nval)) continue;
				double g1 = gtab[is1];
				if (g1 < 1e-10) continue;
				double base = g1*g2;
				const dcomp* src = &Krows[is1][(size_t)j*nval];
				for (size_t ix = 0; ix < n_ab; ++ix) {
					dcomp w = base*ph[ix*kmax_all+j];
					dcomp* dst = &Cblk[ix*nval];
					for (size_t l = 0; l < nval; ++l) dst[l] += w*src[l];
				}
				any = true;
			}
			if (!any) continue;
			double t2 = t_lo + it2*dt;
			for (size_t iy = 0; iy < n_loss; ++iy) {
				dcomp e = std::exp(dcomp(0,-loss_grid[iy]*t2));
				for (size_t ix = 0; ix < n_ab; ++ix) {
					dcomp* a = &A[(ix*n_loss+iy)*nval];
					const dcomp* c = &Cblk[ix*nval];
					for (size_t l = 0; l < nval; ++l) a[l] += e*c[l];
				}
			}
		}
		// Eq. (11)/(15): incoherent over the final valence state, coherent over (t1,t2).
		for (size_t ix = 0; ix < n_ab; ++ix) {
			for (size_t iy = 0; iy < n_loss; ++iy) {
				const dcomp* a = &A[(ix*n_loss+iy)*nval];
				double tot = 0;
				if (coherent_final_states) {
					dcomp csum(0,0);
					for (size_t l = 0; l < nval; ++l) csum += a[l];
					tot = std::norm(csum);
				} else {
					for (size_t l = 0; l < nval; ++l) tot += std::norm(a[l]);
				}
				intensity[d][iy*n_ab+ix] = tot;
			}
		}
		std::cout << "  delay " << d+1 << "/" << delays.size() << " done" << std::endl;
	}

	write_trRIXS_2time(ab_grid, loss_grid, intensity, delays, prefix);
	std::cout << "=== End trRIXS (Mode B) ===" << std::endl << std::endl;
	return;
}

#endif
