#include <stdlib.h>
#include "multiplet.hpp"

using namespace std;
// Includes multiplet interaction

double calc_U(double* gaunt1, double* gaunt2, const double* SC, int size) {
	double u = 0;
	for (int i = 0; i < size; ++i) u += gaunt1[i] * gaunt2[i] * SC[i];
	return u;
}

void calc_ham(Hilbert& hilbs, const HParam& hparam, bool nohyb, Hilbert::HTarget target) {
	// Assemble Hamiltonian of the hilbert space
	for (auto& blk : hilbs.hblks) {
		if (target == Hilbert::HTarget::HAM_STATIC) {
			blk.malloc_ham_static();
		} else {
			if (hilbs.num_ch == 0) blk.malloc_ham(hparam.gs_diag_option);
			if (hilbs.num_ch == 1) blk.malloc_ham(hparam.ex_diag_option);
		}
	}
	calc_coulomb(hilbs,hparam.SC,target);
	if (hilbs.SO_on) {
		calc_SO(hilbs,hparam.SO[0],1,target);
		calc_SO(hilbs,hparam.SO[1],2,target);
	}
	if (hilbs.CF_on) calc_CF(hilbs,&hparam.CF[0],target);
	if (hilbs.CV_on) calc_CV(hilbs,&hparam.FG[0],target);
	if (hilbs.MF_on) calc_MF(hilbs,hparam.MF,target);
	// nohyb=true (only ever used to build the static-only piece for real-time
	// propagation, target=HAM_STATIC) skips the hybridization term entirely --
	// it's rebuilt separately, per-timestep, into ham_hyb (see calc_HYB, tdynamics.hpp).
	if (hilbs.HYB_on && !nohyb) calc_HYB(hilbs,hparam,nohyb,0.0,target);
	return;
}

void calc_coulomb(Hilbert& hilbs, const vector<double*>& SC, Hilbert::HTarget target) {
	//Calculate Coulomb Matrix Element
	for (int i = 0; i < hilbs.atlist.size(); ++i) {
		int l = hilbs.atlist[i].l;
		if (l <= 0 || ed::is_zero_arr(SC[l-1],l*2+1)) continue;
		vector<int> ml_arr(l*2+1,0);
		for (int ml = -l; ml <= l; ++ml) ml_arr[ml+l] = ml;
		for (int msum = -l*2; msum <= l*2; ++msum) {
			vector<pair<int,int>> mpair_p, mpair_a;
			int m = -l;
			while(m < msum-m) {
				if (msum-m <= l) {
					mpair_p.push_back({m,msum-m});
					mpair_a.push_back({m,msum-m});
					mpair_a.push_back({msum-m,m});
				}
				++m;
			}
			if(m == msum-m) mpair_a.push_back({m,msum-m});
			// Antiparallel spin pairs
			for (auto m12 : mpair_a) {
				struct QN qn12[2] = {{m12.first,-0.5,i},{m12.second,0.5,i}}; // psi_kl, lhs 
				for (auto m34 : mpair_a) {
					struct QN qn34[2] = {{m34.first,-0.5,i},{m34.second,0.5,i}}; // psi_ij,rhs
					vpulli entries = hilbs.match(2,qn12,qn34);
					for (auto e : entries) {
						double matelem = calc_U(gaunt(l,m12.first,l,m34.first),
										gaunt(l,m34.second,l,m12.second),SC[l-1],2*l+1);
						if (e.first != e.second) matelem *= hilbs.Fsign(qn12,e.first,2)*hilbs.Fsign(qn34,e.second,2);
						hilbs.fill_hblk(matelem,e.first,e.second,target);
					}
				}
			}
			// Parallel spin pair
			for (auto spin : {-0.5,0.5}) {
				for (auto m12 : mpair_p) { 
					struct QN qn12[2] = {{m12.first,spin,i},{m12.second,spin,i}};
					for (auto m34 : mpair_p) {
						struct QN qn34[2] = {{m34.first,spin,i},{m34.second,spin,i}};
						vpulli entries = hilbs.match(2,qn12,qn34);
						for (auto e : entries) {
							double matelem = calc_U(gaunt(l,m12.first,l,m34.first),
											gaunt(l,m34.second,l,m12.second),SC[l-1],2*l+1) -
										 	calc_U(gaunt(l,m12.second,l,m34.first),
										 	gaunt(l,m34.second,l,m12.first),SC[l-1],2*l+1);
							if (e.first != e.second) matelem *= hilbs.Fsign(qn12,e.first,2)*hilbs.Fsign(qn34,e.second,2);
							hilbs.fill_hblk(matelem,e.first,e.second,target);
						}
					}
				}
			}
		}
		// !!!!!!Shift energy for particle hole symmetry
		// double eshift = hilbs.pheshift(ed::trace(mat,hilbs.hmat_size),2);
		// for (int i = 0; i < hilbs.hmat_size; ++i) mat[i+hilbs.hmat_size*i] -= eshift;
	}
	return;
}

vecd CFmat(int l, const double* CF) {
	// Crystal field present when the 2 state has same spin, hole language
	vecd cfmat;
	try {
		if (l == 2) {
			cfmat = vecd(25,0);
			cfmat[0*5+0] = cfmat[4*5+4] = (CF[0] + CF[2])/2;
			cfmat[0*5+4] = cfmat[4*5+0] = (CF[0] - CF[2])/2;
			cfmat[1*5+1] = cfmat[3*5+3] = (CF[3] + CF[4])/2;
			cfmat[1*5+3] = cfmat[3*5+1] = (CF[4] - CF[3])/2;
			cfmat[2*5+2] = CF[1];
		} else throw invalid_argument("non-d and non-octahedral crystal field not coded yet");
	} catch(const exception &ex) {cout << ex.what() << "\n";}
	return cfmat;
}

void calc_CF(Hilbert& hilbs, const double* CF, Hilbert::HTarget target) {
	// Calculate Crystal Field, Charge Transfer Matrix Element
	for (int i = 0; i < hilbs.atlist.size(); ++i) {
		int l = hilbs.atlist[i].l;
		vecd cfmat;
		if (l == 2 && !hilbs.atlist[i].is_lig) cfmat = CFmat(l, CF);
		else continue;
		for (int j = 0; j < (l*2+1)*(l*2+1); ++j) {
			if (cfmat[j] == 0) continue;
			int ml1 = j%(l*2+1)-l, ml2 = j/(l*2+1)-l;
			for (auto spin : {-0.5,0.5}) {
				QN qn1(ml1,spin,i), qn2(ml2,spin,i);
				vpulli entries = hilbs.match(1,&qn1,&qn2);
				for (auto& e : entries) {
					if (e.first == e.second) hilbs.fill_hblk(cfmat[j],e.first,e.second,target);
					else hilbs.fill_hblk(cfmat[j]*hilbs.Fsign(&qn1,e.first,1) // Negative for hole language
						 *hilbs.Fsign(&qn2,e.second,1),e.first,e.second,target);
				}
			}
		}
	}
}

void calc_SO(Hilbert& hilbs, const double lambda, int l_in, Hilbert::HTarget target) {
	// Calculate Spin Orbit Coupling Matrix Element
	if (lambda == 0) return;
	for (int i = 0; i < hilbs.atlist.size(); ++i) {
		int l = hilbs.atlist[i].l;
		// Only TM 2p/3d core orbitals get spin orbit coupling for now
		if (l != l_in || hilbs.atlist[i].is_lig) continue;
		for (int ml = -l; ml <= l; ++ml) {
			// longitudinal phonon
			for (auto spin : {-0.5,0.5}) {
				QN qn(ml,spin,i);
				vpulli entries = hilbs.match(1,&qn,&qn);
				// No sign traversing issue due to symmetry, spin = 1/2 factored in equation
				double matelem = -lambda * spin * ml;
				for (auto& e : entries) hilbs.fill_hblk(matelem,e.first,e.second,target);
			}
			// transverse phonon
			if (ml < l) {
				QN qn1(ml+1,-0.5,i), qn2(ml,0.5,i);
				vpulli entries = hilbs.match(1,&qn1,&qn2);
				double matelem = -lambda/2 * sqrt((l-ml)*(l+ml+1));
				for (auto& e : entries) {
					hilbs.fill_hblk(matelem*hilbs.Fsign(&qn1,e.first,1)*hilbs.Fsign(&qn2,e.second,1),e.first,e.second,target);
					hilbs.fill_hblk(matelem*hilbs.Fsign(&qn2,e.second,1)*hilbs.Fsign(&qn1,e.first,1),e.second,e.first,target);
				}
			}
		}
	}
}

void calc_MF(Hilbert& hilbs, const double h_field, Hilbert::HTarget target) {
	// Mean-field molecular/exchange field standing in for interatomic superexchange
	// with neighboring magnetic sites (Nag et al., PRL 124, 067202 (2020)). Splits
	// the Ms sublevels of the valence d-shell by h_field per unit Ms, which is what
	// turns the SOC/core-valence-exchange-enabled Delta-Ms=1,2 spin-flip transitions
	// (already present via calc_SO/calc_CV) into resolvable inelastic RIXS peaks
	// instead of leaving them degenerate with the elastic line. Diagonal in Sz, so
	// it does not break the Sz block structure.
	if (h_field == 0) return;
	for (int i = 0; i < hilbs.atlist.size(); ++i) {
		auto& at = hilbs.atlist[i];
		if (at.is_lig || !at.is_val || at.l != 2) continue; // TM valence d-shell only
		for (int ml = -at.l; ml <= at.l; ++ml) {
			for (auto spin : {-0.5,0.5}) {
				QN qn(ml,spin,i);
				vpulli entries = hilbs.match(1,&qn,&qn);
				double matelem = -h_field * spin;
				for (auto& e : entries) hilbs.fill_hblk(matelem,e.first,e.second,target);
			}
		}
	}
	return;
}

void calc_CV(Hilbert& hilbs, const double* FG, Hilbert::HTarget target) {
	// Calculate core valence interactions
	if (!hilbs.is_ex) return;
	for (int ci = 0; ci < hilbs.val_ati; ci++) {
		int vi = hilbs.atlist[ci].vind;
		int cl = hilbs.atlist[ci].l, vl = hilbs.atlist[vi].l;
		for (int vml = -vl; vml <= vl; ++vml) {
		for (int vmr = -vl; vmr <= vl; ++vmr) {
		for (int cml = -cl; cml <= cl; ++cml) {
		for (int cmr = -cl; cmr <= cl; ++cmr) {
			if (vml+cml != vmr+cmr) continue;
			vector<pair<double,double>> spairs = {{0.5,0.5},{-0.5,0.5},{0.5,-0.5},{-0.5,-0.5}};
			for (auto & s : spairs) {
				struct QN qnl[2] = {{vml,s.first,vi},{cml,s.second,ci}};
				struct QN qnr[2] = {{vmr,s.first,vi},{cmr,s.second,ci}};
				vpulli entries = hilbs.match(2,qnl,qnr);
				// FG is {F^0, G^1, F^2, G^3} -- the upstream CTFAMultiplet convention, see
				// https://github.com/khhsu0724/CTFAMultiplet/wiki/Input-Parameters ("F0,G1,F2,G3").
				// calc_U contracts its array by multipole order k, and that ordering is already
				// the right one: the Coulomb channel reads FG[0]=F^0 at k=0 and FG[2]=F^2 at k=2,
				// the exchange channel reads FG[1]=G^1 at k=1 and FG[3]=G^3 at k=3, and the
				// mismatched slots are killed by the Gaunt parity rules (for l_v=2, l_c=1 the
				// Coulomb product is nonzero only for even k and the exchange only for odd k).
				// The Udp printout in main.cpp, F^0 - G^1/15 - 3G^3/70, encodes the same ordering.
				double meF = calc_U(gaunt(vl,vml,vl,vmr),gaunt(cl,cmr,cl,cml),FG,2*cl+1);
				for (auto& e : entries)
					hilbs.fill_hblk(meF*hilbs.Fsign(qnl,e.first,2)*hilbs.Fsign(qnr,e.second,2),e.first,e.second,target);
				qnl[1] = {vml,s.second,vi}, qnl[0] = {cml,s.first,ci};
				qnr[0] = {vmr,s.first,vi}, qnr[1] = {cmr,s.second,ci};
				entries = hilbs.match(2,qnl,qnr);
				double meG = calc_U(gaunt(vl,vmr,cl,cml),gaunt(vl,vml,cl,cmr),FG,cl+vl+1);
				// Core-valence EXCHANGE sign. Cowan eq. 12.6 gives the exchange coefficient as
				// g_k = -(-1)^S <l||C^k||l'>^2 {6j}, i.e. the parallel-spin (triplet) combination
				// is pulled DOWN -- Hund's rule. Getting this sign backwards leaves every term
				// energy numerically correct but interchanges the singlet and triplet of each L,
				// which is exactly what the atomic-limit test showed (2p^5 3d^9 with no CF, no
				// hybridization and no core SOC gave 1F at 4.1138 and 3F at 7.8862, the exact
				// analytic values but swapped). EXCH_SIGN_LEGACY=1 restores the old sign.
				static const bool exch_legacy = (std::getenv("EXCH_SIGN_LEGACY") != nullptr);
				if (exch_legacy) { if (s.first == s.second) meG *= -1; }
				else             { if (s.first != s.second) meG *= -1; }
				for (auto& e : entries)
					hilbs.fill_hblk(meG*hilbs.Fsign(qnl,e.first,2)*hilbs.Fsign(qnr,e.second,2),e.first,e.second,target);
			}
		}}}}
	}
	return;
}

void calc_HYB(Hilbert& hilbs, const HParam& hparam, bool nohyb, double t, Hilbert::HTarget target) {
	// Charge Transfer and Hybridization
	if (hilbs.cluster->no_HYB) return;
	else hilbs.cluster->set_hyb_params(hparam);
	int nvo = hilbs.cluster->vo_persite * hilbs.tot_site_num();
	int nco = hilbs.cluster->co_persite * hilbs.tot_site_num();
	// HYBRIDIZATION is momentum dependent, need to rewrite code here
	vecd hybmat = ed::make_blk_mat(hilbs.cluster->get_tmat_real(t),hilbs.tot_site_num());
	// Diagnostic: the paper-1 reference Fortran builds the core-hole sector's hopping with the
	// OPPOSITE overall sign to the valence sector (its line ~1552 vs ~761: +2*tmatx/tmaty/tmato
	// against -2*tmatx/tmaty/tmato, otherwise byte-identical code). Flipping all three together
	// is a p -> -p gauge transformation only if the p-p term (tmato, the tpp channel) vanishes;
	// with tpp_pi = 0.25 it does not, so the two conventions are genuinely different core-hole
	// Hamiltonians. EX_HOP_SIGN=-1 reproduces the reference's convention for the core-hole
	// (is_ex) space so the two can be compared.
	if (hilbs.is_ex) {
		if (const char* sg = std::getenv("EX_HOP_SIGN")) {
			if (std::atof(sg) < 0) for (auto& v : hybmat) v = -v;
		}
	}
	// Get Hybridization Information, there should be num_orb x num_orb matrix providing hybdrization information
	// Loop through hyb matrix, TODO: Loop through each site
	for (int i = 0; i < nvo; ++i) {
		for (int j = 0; j < nvo; ++j) {
			if (abs(hybmat[i*nvo+j]) < TOL) continue;
			// Fill in spin down matrix element
			ulli lhssd = BIG1 << i, rhssd = BIG1 << j, incsd = (rhssd|lhssd);
			ulli lhssd_shift = (incsd-lhssd) << nco, rhssd_shift = (incsd-rhssd) << nco;	
			ulli lhssu = BIG1 << (i+nvo), rhssu = BIG1 << (j+nvo), incsu = (rhssu|lhssu);
			ulli lhssu_shift = (incsu-lhssu) << (2*nco), rhssu_shift = (incsu-rhssu) << (2*nco);
			QN qnlu, qnru, qnld, qnrd;
			int atlist_ind = 0;
			for (auto& at : hilbs.atlist) {
				if ((i+nco) >= at.sind && (i+nco) <= at.eind) {
					int ml = at.val_l - i + at.sind - nco;
					qnlu = QN(ml,0.5,atlist_ind);
					qnld = QN(ml,-0.5,atlist_ind);
				}
				if ((j+nco) >= at.sind && (j+nco) <= at.eind) {
					int ml = at.val_l - j + at.sind - nco;
					qnru = QN(ml,0.5,atlist_ind);
					qnrd = QN(ml,-0.5,atlist_ind);
				}
				atlist_ind++;
			}
			// Fill in spin down matrix element
			vector<ulli> match_entries = hilbs.enum_hspace(incsd,0,ed::count_bits(incsd)-1,0);
			vpulli entries;
			for (auto &me : match_entries) {
				entries.emplace_back(pair<ulli,ulli>(me-lhssd_shift,me-rhssd_shift)); // We can probably fill in the matrix at this point
			}

			for (auto &e : entries) {
if (e.first == e.second) hilbs.fill_hblk(hybmat[i*nvo+j],e.first,e.second,target);
				else hilbs.fill_hblk(hybmat[i*nvo+j]*hilbs.Fsign(&qnld,e.first,1)*hilbs.Fsign(&qnrd,e.second,1),
									 e.first,e.second,target);
			}
			// Fill in spin up matrix element
			match_entries = hilbs.enum_hspace(incsu,0,ed::count_bits(incsu)-1,0);
			entries.clear();
			for (auto &me : match_entries) {
				entries.emplace_back(pair<ulli,ulli>(me-lhssu_shift,me-rhssu_shift));
			}
			for (auto &e : entries) {
if (e.first == e.second) hilbs.fill_hblk(hybmat[i*nvo+j],e.first,e.second,target);
				else hilbs.fill_hblk(hybmat[i*nvo+j]*hilbs.Fsign(&qnlu,e.first,1)*hilbs.Fsign(&qnru,e.second,1),
									 e.first,e.second,target);
			}
		}
	}
	return;
}
