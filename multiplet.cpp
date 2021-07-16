#include <iostream>
#include <fstream>
#include <cstdlib> 
#include <string>
#include <vector>
#include <utility>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include "diagonalize.h"
#include "gaunt.h"
#include "hilbert.h"

using namespace std;
// Includes multiplet interaction calculation for a specific orbital

double calc_U(double* gaunt1, double* gaunt2, double* SC, int size) {
	double u = 0;
	for (int i = 0; i < size; ++i) {
		u += gaunt1[i] * gaunt2[i] * SC[i];
	}
	return u;
}

void calc_coulomb(Hilbert& hspace, double* mat, double* SC) {
	//Calculate Coulomb Matrix Element
	if (hspace.occ_num < 2 || hspace.occ_num > hspace.orb_avail - 2) return;
	int n = hspace.hmat_size, m = n;
	int* ml_arr = new int[hspace.l*2+1]{0};
	for (int ml = -hspace.l; ml <= hspace.l; ++ml) {
		ml_arr[ml+hspace.l] = ml;
	}
	// Find all combination when m1+m2=m3+m4
	for (int msum = -hspace.l*2; msum <= hspace.l*2; ++msum) {
		vector<pair<int,int>> mpair_p; //parallel ml pairs
		vector<pair<int,int>> mpair_a; //anti-parallel ml pairs
		int m = -hspace.l;
		while(m < msum-m) {
			if (msum-m <= hspace.l) {
				mpair_p.push_back({m,msum-m});
				mpair_a.push_back({m,msum-m});
				mpair_a.push_back({msum-m,m});
			}
			++m;
		}
		if(m == msum-m) {
			//antiparallel spin, can occupy same ml
			mpair_a.push_back({m,msum-m});
		}

		for (auto m12 : mpair_a) {
			// Populate Hamiltonian that contains the antiparallel spin pair
			struct QN qn12[2] = {{m12.first,-1},{m12.second,1}}; // psi_kl, lhs 
			for (auto m34 : mpair_a) {
				struct QN qn34[2] = {{m34.first,-1},{m34.second,1}}; // psi_ij,rhs
				vector<pair<int,int>> entries = hspace.match_states(2, qn12, qn34);
				for (auto e : entries) {
					double matelem = calc_U(gaunt(hspace.l,m12.first,hspace.l,m34.first),
									gaunt(hspace.l,m34.second,hspace.l,m12.second),SC,2*hspace.l+1);
					if (hspace.sindex(e.first) != hspace.sindex(e.second)) { // Off Diagonal Term
						matelem *= hspace.Psign(qn12,qn34,e.first,e.second,2,2);
					}
					hspace.fill_in_mat(mat,matelem,hspace.sindex(e.first),hspace.sindex(e.second));
					// mat[hspace.sindex(e.first)+hspace.hmat_size*hspace.sindex(e.second)] += matelem;
				}
			}
		}

		for (auto spin : {-1,1}) {
			for (auto m12 : mpair_p) {
				// Populate Hamiltonian that contains the parallel spin pair
				struct QN qn12[2] = {{m12.first,spin},{m12.second,spin}}; // psi_kl, lhs 
				for (auto m34 : mpair_p) {
					struct QN qn34[2] = {{m34.first,spin},{m34.second,spin}}; // psi_ij,rhs
					vector<pair<int,int>> entries = hspace.match_states(2, qn12, qn34);
					for (auto e : entries) {
						double matelem = calc_U(gaunt(hspace.l,m12.first,hspace.l,m34.first),
									 	gaunt(hspace.l,m34.second,hspace.l,m12.second),SC,2*hspace.l+1) -
									 	calc_U(gaunt(hspace.l,m12.second,hspace.l,m34.first),
									 	gaunt(hspace.l,m34.second,hspace.l,m12.first),SC,2*hspace.l+1);
						if (hspace.sindex(e.first) != hspace.sindex(e.second)) { // Off Diagonal Term
							matelem *= hspace.Psign(qn12,qn34,e.first,e.second,2,2);
						}
						hspace.fill_in_mat(mat,matelem,hspace.sindex(e.first),hspace.sindex(e.second));
						// mat[hspace.sindex(e.first)+hspace.hmat_size*hspace.sindex(e.second)] += matelem;
					}
				}
			}
		}
	}
	// Shift energy for particle hole symmetry
	double t = trace(mat,hspace.hmat_size,hspace.hmat_size);
	double eshift = hspace.pheshift(t,2);
	for (int i = 0; i < hspace.hmat_size; ++i) mat[i+hspace.hmat_size*i] -= eshift;
}

double* CFmat(int l, double tenDQ, int coord) {
	// Crystal field present when the 2 state has same spin
	try {
		if (l == 2 && coord == 6) {
			double* cfmat = new double[25]{0};
			cfmat[0] = 0.1 * tenDQ;
			cfmat[1+5*1] = -0.4 * tenDQ;
			cfmat[2+5*2] = 0.6 * tenDQ;
			cfmat[3+5*3] = -0.4 * tenDQ;
			cfmat[4+5*4] = 0.1 * tenDQ;
			cfmat[4+5*0] = 0.5 * tenDQ;
			cfmat[0+5*4] = 0.5 * tenDQ;
			return cfmat;
		} else {
			throw invalid_argument("non-d and non-octahedral crystal field not coded yet");
		}
	} catch(const exception &ex) {
			std::cout << ex.what() << "\n";
	}
	double* ph; // Place holder
	return ph;
}

void calc_CF(Hilbert& hspace, double* mat, double tenDQ, int coord) {
	//Calculate Crystal Field Matrix Element, octahedral symmetry
	if (hspace.occ_num == 0 || hspace.occ_num == hspace.orb_avail) return;
	double* cfmat = CFmat(hspace.l, tenDQ, coord);
	for (int i = 0; i < (hspace.l*2+1)*(hspace.l*2+1); ++i) {
		if (cfmat[i] != 0) {
			int ml1 = i%(hspace.l*2+1)-2;
			int ml2 = i/(hspace.l*2+1)-2;
			// Same spin creates crystal field
			for (auto spin : {-1,1}) {
				struct QN qn1[1] = {{ml1,spin}};
				struct QN qn2[1] = {{ml2,spin}};
				vector<pair<int,int>> entries = hspace.match_states(1, qn1, qn2);
				for (auto& e : entries) {
					hspace.fill_in_mat(mat,cfmat[i],hspace.sindex(e.first),hspace.sindex(e.second));
				}
			}
		}
	}
}

void calc_SO(Hilbert& hspace, double* mat, double lambda) {
	// Calculate Spin Orbit Coupling Matrix Element
	if (hspace.occ_num == 0 || hspace.occ_num == hspace.orb_avail) return;
	for (int ml = -hspace.l; ml <= hspace.l; ++ml) {
		// longitudinal
		for (auto spin: {-1,1}) {
			struct QN qn[1] = {{ml,spin}};
			vector<pair<int,int>> entries = hspace.match_states(1, qn, qn);
			for (auto& e : entries) {
				double matelem = lambda/2 * spin * ml; // No sign traversing issue due to symmetry
				hspace.fill_in_mat(mat,matelem,hspace.sindex(e.first),hspace.sindex(e.second));
			}
		}
		// transverse
		if (ml < hspace.l) {
			struct QN qn1[1] = {{ml+1,-1}};
			struct QN qn2[1] = {{ml,1}};
			vector<pair<int,int>> entries = hspace.match_states(1, qn1, qn2);
			for (auto& e : entries) {
				double matelem = lambda/2 * sqrt((hspace.l-ml)*(hspace.l+ml+1));
				hspace.fill_in_mat(mat,matelem * hspace.Psign(qn1,qn2,e.first,e.second,1,1)
					,hspace.sindex(e.first),hspace.sindex(e.second));
				hspace.fill_in_mat(mat,matelem * hspace.Psign(qn2,qn1,e.second,e.first,1,1)
					,hspace.sindex(e.second),hspace.sindex(e.first));
			}
		}
	}
}