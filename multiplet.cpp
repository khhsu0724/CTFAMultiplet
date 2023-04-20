#include <stdlib.h>
#include "multiplet.hpp"

using namespace std;
// Includes multiplet interaction

double calc_U(double* gaunt1, double* gaunt2, const double* SC, int size) {
	double u = 0;
	for (int i = 0; i < size; ++i) u += gaunt1[i] * gaunt2[i] * SC[i];
	return u;
}

void calc_ham(Hilbert& hilbs, const HParam& hparam) {
	// Assemble Hamiltonian of the hilbert space
	for (auto& blk : hilbs.hblks) blk.malloc_ham(hparam.diag_option);
	calc_coulomb(hilbs,hparam.SC); 
	if (hilbs.SO_on) {
		calc_SO(hilbs,hparam.SO[0],1);
		calc_SO(hilbs,hparam.SO[1],2);
	}
	if (hilbs.CF_on) calc_CF(hilbs,&hparam.CF[0]);
	if (hilbs.CV_on) calc_CV(hilbs,&hparam.FG[0]);
	if (hilbs.HYB_on) calc_HYB(hilbs,hparam);
	return;
}

void calc_coulomb(Hilbert& hilbs, const vector<double*>& SC) {
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
						hilbs.fill_hblk(matelem,e.first,e.second);
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
							hilbs.fill_hblk(matelem,e.first,e.second);
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
			cfmat[24] = cfmat[0] = (CF[0] + CF[2])/2;
			cfmat[20] = cfmat[4] = (CF[0] - CF[2])/2;
			cfmat[18] = cfmat[6] = (CF[3] + CF[4])/2;
			cfmat[16] = cfmat[8] = (CF[4] - CF[3])/2;
			cfmat[12] = CF[1];
		} else throw invalid_argument("non-d and non-octahedral crystal field not coded yet");
	} catch(const exception &ex) {cout << ex.what() << "\n";}
	return cfmat;
}

void calc_CF(Hilbert& hilbs, const double* CF) {
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
					if (e.first == e.second) hilbs.fill_hblk(cfmat[j],e.first,e.second);
					else hilbs.fill_hblk(cfmat[j]*hilbs.Fsign(&qn1,e.first,1) // Negative for hole language
						 *hilbs.Fsign(&qn2,e.second,1),e.first,e.second);
				}
			}
		}
	}
}

void calc_SO(Hilbert& hilbs, const double lambda, int l_in) {
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
				for (auto& e : entries) hilbs.fill_hblk(matelem,e.first,e.second);
			}
			// transverse phonon
			if (ml < l) {
				QN qn1(ml+1,-0.5,i), qn2(ml,0.5,i);
				vpulli entries = hilbs.match(1,&qn1,&qn2);
				double matelem = -lambda/2 * sqrt((l-ml)*(l+ml+1));
				for (auto& e : entries) {
					hilbs.fill_hblk(matelem*hilbs.Fsign(&qn1,e.first,1)*hilbs.Fsign(&qn2,e.second,1),e.first,e.second);
					hilbs.fill_hblk(matelem*hilbs.Fsign(&qn2,e.second,1)*hilbs.Fsign(&qn1,e.first,1),e.second,e.first);
				}
			}
		}
	}
}

void calc_CV(Hilbert& hilbs, const double* FG) {
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
				double meF = calc_U(gaunt(vl,vml,vl,vmr),gaunt(cl,cmr,cl,cml),FG,2*cl+1);
				for (auto& e : entries) 
					hilbs.fill_hblk(meF*hilbs.Fsign(qnl,e.first,2)*hilbs.Fsign(qnr,e.second,2),e.first,e.second);
				qnl[1] = {vml,s.second,vi}, qnl[0] = {cml,s.first,ci};
				qnr[0] = {vmr,s.first,vi}, qnr[1] = {cmr,s.second,ci};
				entries = hilbs.match(2,qnl,qnr);
				double meG = calc_U(gaunt(vl,vmr,cl,cml),gaunt(vl,vml,cl,cmr),FG,cl+vl+1);
				if (s.first == s.second) meG *= -1; // This is due to the ordering convention of the fermions in this code
				for (auto& e : entries)
					hilbs.fill_hblk(meG*hilbs.Fsign(qnl,e.first,2)*hilbs.Fsign(qnr,e.second,2),e.first,e.second);
			}
		}}}}
	}
	return;
}

void calc_HYB(Hilbert& hilbs, const HParam& hparam) {
	// Charge Transfer and Hybridization
	if (hilbs.cluster->no_HYB) return;
	else hilbs.cluster->set_hyb_params(hparam);
	int nvo = hilbs.cluster->vo_persite * hilbs.tot_site_num();
	int nco = hilbs.cluster->co_persite * hilbs.tot_site_num();
	// HYBRIDIZATION is momentum dependent, need to rewrite code here
	vecd hybmat = ed::make_blk_mat(hilbs.cluster->get_tmat_real(),hilbs.tot_site_num());
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
				if (e.first == e.second) hilbs.fill_hblk(hybmat[i*nvo+j],e.first,e.second);
				else hilbs.fill_hblk(hybmat[i*nvo+j]*hilbs.Fsign(&qnld,e.first,1)*hilbs.Fsign(&qnrd,e.second,1),
									 e.first,e.second);
			}
			// Fill in spin up matrix element
			match_entries = hilbs.enum_hspace(incsu,0,ed::count_bits(incsu)-1,0);
			entries.clear();
			for (auto &me : match_entries) {
				entries.emplace_back(pair<ulli,ulli>(me-lhssu_shift,me-rhssu_shift));
			}
			for (auto &e : entries) {
				if (e.first == e.second) hilbs.fill_hblk(hybmat[i*nvo+j],e.first,e.second);
				else hilbs.fill_hblk(hybmat[i*nvo+j]*hilbs.Fsign(&qnlu,e.first,1)*hilbs.Fsign(&qnru,e.second,1),
									 e.first,e.second);
			}
		}
	}
	return;
}
