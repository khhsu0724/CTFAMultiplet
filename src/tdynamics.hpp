#include "multiplet.hpp"
#ifndef TDYNAMICS
#define TDYNAMICS

// Real-time propagation of many-body state vectors under a (possibly pump-modulated,
// time-dependent) Hamiltonian, via a short-iterative-Lanczos (Krylov subspace) exponential
// integrator. This is the time-domain analog of the existing solver=4 (Lanczos/ContFracExpan)
// machinery: it never forms or diagonalizes the full N-dimensional Hamiltonian, only ever
// running an m-step (m~10-30) Lanczos recursion per timestep.
//
// All times are natural units (hbar=1), i.e. eV^-1 (1 eV^-1 = 0.658 fs), matching the energy
// units (eV) used everywhere else in this code.

// Propagates psi <- exp(-i*H*dt)*psi via an m-step Lanczos run on ham (any Matrix<T>, e.g.
// a Block's ham/ham_static, or a CombinedHam view). Reconstructs the propagator from the
// small (m x m) real symmetric tridiagonal matrix built by Lanczos(): diagonalize it once
// (cheap dense diagonalization, size m not N), then exponentiate its eigenvalues.
template <typename T>
void LanczosPropagate(Matrix<T>* ham, vecc& psi, double dt, int m = 20) {
	double vnorm = ed::norm(psi);
	if (vnorm < TOL) return; // nothing to propagate
	vecc alpha(m,0), betha(m,0);
	std::vector<vecc> V;
	V.reserve(m);
	int m_eff = Lanczos(ham, psi, alpha, betha, m, &V);
	// Lanczos()'s pre-existing breakdown convention returns niter_CFE=i (not i+1) when
	// betha[i+1] vanishes at iteration i -- i.e. it can report m_eff=0 even though
	// alpha[0]/V[0] were already validly computed (this happens whenever psi is already
	// very close to an exact eigenvector of ham, e.g. propagating the ground state with
	// the pump off). Handle that degenerate case directly: it's just a pure phase.
	if (m_eff <= 0) {
		dcomp phase = std::exp(dcomp(0,-alpha[0].real()*dt));
		#pragma omp parallel for
		for (int n = 0; n < (int)psi.size(); ++n) psi[n] = vnorm * phase * V[0][n];
		return;
	}

	// Real symmetric tridiagonal T_m from the Lanczos recursion. H is real-symmetric, so
	// alpha/betha should be real to numerical precision -- this also acts as an implicit
	// Hermiticity sanity check on whatever Hamiltonian is passed in.
	std::vector<double> Tm(m_eff*m_eff, 0.0);
	for (int i = 0; i < m_eff; ++i) Tm[i*m_eff+i] = alpha[i].real();
	for (int i = 0; i < m_eff-1; ++i) Tm[i*m_eff+(i+1)] = Tm[(i+1)*m_eff+i] = betha[i+1].real();

	std::vector<double> eigval(m_eff,0.0), eigvec(m_eff*m_eff,0.0);
	ed_dsyevr(Tm.data(), eigvec.data(), eigval.data(), (size_t)m_eff);

	// exp(-i*T_m*dt) applied to e_1: c_i = sum_k U[i,k] * exp(-i*lambda_k*dt) * U[0,k],
	// using the same eigvec[k*n+i] convention as everywhere else in this codebase
	// (eigvec[k*n+i] = component i of eigenvector k).
	vecc c(m_eff, 0);
	for (int i = 0; i < m_eff; ++i)
		for (int k = 0; k < m_eff; ++k)
			c[i] += eigvec[k*m_eff+i] * std::exp(dcomp(0,-eigval[k]*dt)) * eigvec[k*m_eff+0];

	int N = ham->get_mat_dim();
	std::fill(psi.begin(), psi.end(), dcomp(0));
	for (int i = 0; i < m_eff; ++i) {
		dcomp coeff = vnorm * c[i];
		#pragma omp parallel for
		for (int n = 0; n < N; ++n) psi[n] += coeff * V[i][n];
	}
	return;
}

// Wipes a block's ham_hyb before refilling it for the next timestep. Dense::clear_mat()
// fully deallocates its internal array (unlike EZSparse::clear_mat(), which just empties
// its index/value vectors but stays usable) so a re-malloc is required in that case before
// the next round of fill_mat calls.
inline void clear_ham_hyb(Block<double>& blk) {
	blk.ham_hyb->clear_mat();
	blk.ham_hyb->malloc((int)blk.size);
	return;
}

// One-time setup for a Hilbert space that will be propagated in real time. Builds the
// time-independent piece of H(t) (Coulomb+SOC+CF+CV+MF) into every block's ham_static, and
// allocates (but does not yet fill) every block's ham_hyb. This deliberately never touches
// Block::ham -- the equilibrium ground state used to seed propagation must be found first,
// via the ordinary unmodified calc_ham(hilbs,hparam)/Block::diagonalize() pipeline, before
// calling this function.
inline void init_time_dependent_ham(Hilbert& hilbs, const HParam& hparam) {
	calc_ham(hilbs, hparam, /*nohyb=*/true, Hilbert::HTarget::HAM_STATIC);
	for (auto& blk : hilbs.hblks) blk.malloc_ham_dyn();
	return;
}

// A(t) = A0 * exp(-(t-t0)^2/(2*sigma^2)) * sin(freq*(t-t0)), computed directly from HParam
// (no Cluster instance needed -- matches Cluster::pump_shift()'s formula exactly, duplicated
// here so the pump profile can be written before/independent of any Cluster's own state).
inline double pump_shift_direct(const HParam& hparam, double t, int dir) {
	if (hparam.pump_sigma <= 0) return 0.0;
	double dtc = t - hparam.pump_t0;
	double A0 = (dir == 0) ? hparam.pump_amp_x : (dir == 1) ? hparam.pump_amp_y : hparam.pump_amp_z;
	return A0 * std::exp(-dtc*dtc/(2*hparam.pump_sigma*hparam.pump_sigma)) * std::sin(hparam.pump_freq*dtc);
}

// Writes the pump's Peierls-shift profile A_x(t), A_y(t), A_z(t) over [t_lo,t_hi] at step dt
// -- a plain reference file for plotting, independent of whatever propagation/spectroscopy
// mode is actually run.
inline void write_pump_profile(const HParam& hparam, double t_lo, double t_hi, double dt,
								std::string fname = "pump_profile.txt") {
	std::ofstream fout(fname);
	fout << std::setw(15) << "time (eV^-1)" << std::setw(15) << "A_x(t)"
			<< std::setw(15) << "A_y(t)" << std::setw(15) << "A_z(t)" << std::endl;
	int nsteps = (dt > 0) ? (int)std::round((t_hi-t_lo)/dt) : 0;
	for (int i = 0; i <= nsteps; ++i) {
		double t = t_lo + i*dt;
		fout << std::setw(15) << t << std::setw(15) << pump_shift_direct(hparam,t,0)
				<< std::setw(15) << pump_shift_direct(hparam,t,1)
				<< std::setw(15) << pump_shift_direct(hparam,t,2) << std::endl;
	}
	std::cout << "Wrote pump profile: " << fname << std::endl;
	return;
}

// Advances psi (a state vector living entirely within block blk_ind -- states never leave
// their Sz/Jz block under the pump, since the hybridization term is spin-diagonal, see
// calc_HYB) from t_now to t_now+dt under H(t) = ham_static + ham_hyb(t), with ham_hyb
// evaluated at the midpoint t_now+dt/2. calc_HYB is not scoped to a single block internally,
// so every block's ham_hyb is rebuilt on each call (cheap: it's a one-body term, not the
// two-body Coulomb term) -- a possible future optimization is scoping calc_HYB to only the
// block(s) actually being propagated, not required for correctness here.
inline void propagate_block_state(Hilbert& hilbs, size_t blk_ind, vecc& psi, double t_now,
									double dt, const HParam& hparam, int m = 20) {
	for (auto& blk : hilbs.hblks) {
		if (!blk.ham_static || !blk.ham_hyb)
			throw std::runtime_error("propagate_block_state: call init_time_dependent_ham first");
		clear_ham_hyb(blk);
	}
	// Guard exactly like calc_ham does -- calc_HYB itself does not check HYB_on (that check
	// lives in calc_ham), so calling it unconditionally here would add a hybridization term
	// that was never part of the original equilibrium Hamiltonian whenever HYB=false.
	if (hilbs.HYB_on)
		calc_HYB(hilbs, hparam, /*nohyb=*/false, t_now + dt/2.0, Hilbert::HTarget::HAM_HYB);
	Block<double>& blk = hilbs.hblks[blk_ind];
	CombinedHam<double> Ht(blk.ham_static, blk.ham_hyb);
	LanczosPropagate(&Ht, psi, dt, m);
	ed::norm_vec(psi); // correct residual numerical drift; propagation should already
						// conserve norm to the Lanczos subspace's accuracy (see validation
						// checklist -- track pre-renormalization norm as a diagnostic)
	return;
}

// Orbital occupation numbers <n_c> (real-orbital basis, same order as Cluster::orb_names)
// for an arbitrary state vector psi living in Hilbert block blk_ind -- NOT restricted to
// an already-diagonalized eigenvector. Mirrors photon.cpp's occupation() exactly (same
// annihilation-operator/Hash/Fsign machinery, same seph2real_mat rotation), generalized to
// read amplitudes directly from psi instead of Block::eigvec, and using complex arithmetic
// throughout (std::norm for |amplitude|^2, since psi is complex during propagation).
inline vecd occupation_raw(Hilbert& hilbs, size_t blk_ind, const vecc& psi) {
	int nvo = hilbs.cluster->vo_persite, nco = hilbs.cluster->co_persite;
	vecc U = hilbs.cluster->get_seph2real_mat();
	vecd occ_totsu(nvo,0), occ_totsd(nvo,0);
	Hilbert minus_1vh(hilbs,-1);
	Block<double>& blk = hilbs.hblks[blk_ind];
	for (int c = 0; c < nvo; ++c) {
		vecc wvfncsu(minus_1vh.hsize,0);
		vecc wvfncsd(minus_1vh.hsize,0);
		for (int ci = 0; ci < nvo; ++ci) {
			if (std::abs(U[c*nvo+ci]) < TOL) continue;
			for (size_t j = 0; j < blk.size; ++j) {
				dcomp sj = psi[j];
				if (std::abs(sj) < TOL) continue;
				ulli op_state = hilbs.Hashback(bindex(blk_ind,j));
				ulli opsu = BIG1 << (ci+nco), opsd = BIG1 << (ci+2*nco+nvo);
				if (opsu & op_state) {
					bindex ind = minus_1vh.Hash(op_state-opsu);
					double fsgn = hilbs.Fsign(&opsu,op_state,1);
					wvfncsu[ind.second] += fsgn*sj*U[c*nvo+ci];
				}
				if (opsd & op_state) {
					bindex ind = minus_1vh.Hash(op_state-opsd);
					double fsgn = hilbs.Fsign(&opsd,op_state,1);
					wvfncsd[ind.second] += fsgn*sj*U[c*nvo+ci];
				}
			}
		}
		for (auto& w : wvfncsu) occ_totsu[c] += std::norm(w);
		for (auto& w : wvfncsd) occ_totsd[c] += std::norm(w);
	}
	vecd occ_tot(nvo,0);
	for (int i = 0; i < nvo; ++i) occ_tot[i] = occ_totsu[i] + occ_totsd[i];
	return occ_tot;
}

// Propagates the ground state from t_min to t_max (pump centered at hparam.pump_t0, which
// should lie inside this window), writing orbital occupations vs. time to fname and to
// stdout every print_every steps -- the direct analog of PRX Fig. 2(b1)/(c1). Requires hilbs
// already built+diagonalized via the ordinary equilibrium pipeline, as with validate_propagator.
inline void scan_td_occupation(Hilbert& hilbs, const HParam& hparam, double t_min, double t_max,
								std::string fname = "occupation_vs_time.txt", int print_every = 20) {
	std::cout << std::endl << "=== Time-dependent occupation scan ===" << std::endl;

	double gs_en = hilbs.hblks[0].eig[0];
	size_t gs_blk = 0, gs_ind = 0;
	for (size_t b = 0; b < hilbs.hblks.size(); ++b)
		for (size_t i = 0; i < hilbs.hblks[b].nev; ++i)
			if (hilbs.hblks[b].eig[i] < gs_en) {gs_en = hilbs.hblks[b].eig[i]; gs_blk = b; gs_ind = i;}

	Block<double>& blk = hilbs.hblks[gs_blk];
	size_t N = blk.size;
	vecc psi(N,0);
	for (size_t i = 0; i < N; ++i) psi[i] = dcomp(blk.eigvec[gs_ind*N+i],0);

	init_time_dependent_ham(hilbs, hparam);

	double dt = hparam.td_dt;
	int m = hparam.td_m;
	int nsteps = (int)std::round((t_max - t_min)/dt);
	auto& names = hilbs.cluster->get_orb_names();
	int nvo = hilbs.cluster->vo_persite;

	std::ofstream fout(fname);
	fout << std::setw(12) << "t";
	for (auto& nm : names) fout << std::setw(12) << nm;
	// Energy absorbed from the pump. E_inst uses H(t) (pump-modulated hybridization at t);
	// E_stat uses the A=0 Hamiltonian, so after the pulse E_stat is conserved and its rise
	// above the static ground-state energy IS the energy deposited. That deposited energy is
	// what shifts the absorption resonance: the t1-side resonance sits at w_in = E_ch - E_n,
	// so a pumped state sitting dE higher red-shifts the whole incident-energy axis by dE.
	fout << std::setw(14) << "E_inst" << std::setw(14) << "E_stat" << std::endl;
	std::cout << std::setw(10) << "t";
	for (auto& nm : names) std::cout << std::setw(10) << nm;
	std::cout << std::endl;

	// t_min is very likely far enough before the pump (centered at hparam.pump_t0) that psi
	// is still (to a very good approximation) the equilibrium ground state there.
	double t = t_min;
	for (int s = 0; s <= nsteps; ++s) {
		if (s % print_every == 0 || s == nsteps) {
			vecd occ = occupation_raw(hilbs, gs_blk, psi);
			fout << std::setw(12) << t;
			for (int c = 0; c < nvo; ++c) fout << std::setw(12) << occ[c];
			// <psi|H(t_snap)|psi>; snapshot_ham_hyb/instantaneous_energy live in tr_photon.hpp
			// (which includes this file), so the same clear+refill+mvmult is inlined here.
			auto energy_at = [&](double t_snap) {
				for (auto& blk : hilbs.hblks) clear_ham_hyb(blk);
				if (hilbs.HYB_on)
					calc_HYB(hilbs, hparam, /*nohyb=*/false, t_snap, Hilbert::HTarget::HAM_HYB);
				CombinedHam<double> Ht(hilbs.hblks[gs_blk].ham_static, hilbs.hblks[gs_blk].ham_hyb);
				vecc Hpsi(psi.size(),0);
				Ht.mvmult_cmplx(psi, Hpsi);
				return ed::inner_product(psi,Hpsi).real();
			};
			double e_inst = energy_at(t);
			double e_stat = energy_at(hparam.pump_t0 + 1.0e4*std::max(hparam.pump_sigma,1.0));
			fout << std::setw(14) << e_inst << std::setw(14) << e_stat << std::endl;
			if (s % (print_every*5) == 0 || s == nsteps) {
				std::cout << std::setw(10) << t;
				for (int c = 0; c < nvo; ++c) std::cout << std::setw(10) << occ[c];
				std::cout << std::endl;
			}
		}
		if (s == nsteps) break;
		for (auto& b : hilbs.hblks) clear_ham_hyb(b);
		if (hilbs.HYB_on)
			calc_HYB(hilbs, hparam, /*nohyb=*/false, t + dt/2.0, Hilbert::HTarget::HAM_HYB);
		CombinedHam<double> Ht(blk.ham_static, blk.ham_hyb);
		LanczosPropagate(&Ht, psi, dt, m);
		ed::norm_vec(psi);
		t += dt;
	}
	std::cout << "Wrote " << fname << std::endl;
	std::cout << "=== End occupation scan ===" << std::endl << std::endl;
	return;
}

// Self-contained diagnostic harness (Phase 1-2 validation checklist items 2-5). Requires
// hilbs to already be built+diagonalized via the ordinary equilibrium pipeline
// (calc_ham(hilbs,hparam) + Block::diagonalize()) -- exactly what main.cpp already does
// for GS before this is called. Propagates the true ground state for a fixed number of
// steps and reports norm drift (pre-renormalization), overlap with the initial state
// (should stay a pure phase if the pump is off), and instantaneous energy.
inline void validate_propagator(Hilbert& hilbs, const HParam& hparam, int nsteps = 200) {
	std::cout << std::endl << "=== Time-dependent propagator validation ===" << std::endl;

	double gs_en = hilbs.hblks[0].eig[0];
	size_t gs_blk = 0, gs_ind = 0;
	for (size_t b = 0; b < hilbs.hblks.size(); ++b)
		for (size_t i = 0; i < hilbs.hblks[b].nev; ++i)
			if (hilbs.hblks[b].eig[i] < gs_en) {gs_en = hilbs.hblks[b].eig[i]; gs_blk = b; gs_ind = i;}
	std::cout << "Ground state block: " << gs_blk << ", energy: " << gs_en << std::endl;

	Block<double>& blk = hilbs.hblks[gs_blk];
	size_t N = blk.size;
	vecc psi0(N,0);
	for (size_t i = 0; i < N; ++i) psi0[i] = dcomp(blk.eigvec[gs_ind*N+i],0);

	init_time_dependent_ham(hilbs, hparam);

	vecc psi(psi0);
	double t = 0.0, dt = hparam.td_dt;
	int m = hparam.td_m;
	std::cout << "dt=" << dt << ", Krylov dim m=" << m << ", steps=" << nsteps << std::endl;
	std::cout << std::setw(10) << "t" << std::setw(16) << "norm(pre-renorm)"
				<< std::setw(16) << "|<psi0|psi(t)>|" << std::setw(16) << "phase"
				<< std::setw(16) << "phase (expected)" << std::endl;
	for (int s = 0; s < nsteps; ++s) {
		// Mirror propagate_block_state but keep the pre-renormalization norm for diagnostics.
		for (auto& b : hilbs.hblks) clear_ham_hyb(b);
		if (hilbs.HYB_on)
			calc_HYB(hilbs, hparam, /*nohyb=*/false, t + dt/2.0, Hilbert::HTarget::HAM_HYB);
		CombinedHam<double> Ht(blk.ham_static, blk.ham_hyb);
		LanczosPropagate(&Ht, psi, dt, m);
		double norm = ed::norm(psi);
		ed::norm_vec(psi);
		t += dt;
		if (s % std::max(1,nsteps/10) == 0 || s == nsteps-1) {
			dcomp overlap(0,0);
			for (size_t i = 0; i < N; ++i) overlap += std::conj(psi0[i])*psi[i];
			std::cout << std::setw(10) << t << std::setw(16) << norm
						<< std::setw(16) << std::abs(overlap) << std::setw(16) << std::arg(overlap)
						<< std::setw(16) << -gs_en*t << std::endl;
		}
	}
	std::cout << "=== End validation ===" << std::endl << std::endl;
	return;
}

#endif
