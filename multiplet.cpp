#include <stdlib.h>
#include "multiplet.hpp"

using namespace std;

// Includes multiplet interaction

double calc_U(double* gaunt1, double* gaunt2, double* SC, int size) {
	double u = 0;
	for (int i = 0; i < size; ++i) u += gaunt1[i] * gaunt2[i] * SC[i];
	return u;
}

void calc_coulomb(Hilbert& hilbs, double* SC) {
	//Calculate Coulomb Matrix Element
	for (int i = 0; i < hilbs.atlist.size(); ++i) {
		int l = hilbs.atlist[i].l;
		int* ml_arr = new int[l*2+1]{0};
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
										gaunt(l,m34.second,l,m12.second),SC,2*l+1);
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
											gaunt(l,m34.second,l,m12.second),SC,2*l+1) -
										 	calc_U(gaunt(l,m12.second,l,m34.first),
										 	gaunt(l,m34.second,l,m12.first),SC,2*l+1);
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
}

double* CFmat(int l, double tenDQ, int coord) {
	// Crystal field present when the 2 state has same spin, hole language
	double* cfmat;
	try {
		if (l == 2 && coord == 6) {
			cfmat = new double[25]{0};
			cfmat[0] = -0.1 * tenDQ;
			cfmat[1+5*1] = 0.4 * tenDQ;
			cfmat[2+5*2] = -0.6 * tenDQ;
			cfmat[3+5*3] = 0.4 * tenDQ;
			cfmat[4+5*4] = -0.1 * tenDQ;
			cfmat[4+5*0] = -0.5 * tenDQ;
			cfmat[0+5*4] = -0.5 * tenDQ;
		} else throw invalid_argument("non-d and non-octahedral crystal field not coded yet");
	} catch(const exception &ex) {cout << ex.what() << "\n";}
	return cfmat;
}

double* CTmat(int l, double del) {
	// Need more flexible implementation in the future
	double* ctmat;
	try {
		if (l == 1) {
			ctmat = new double[9]{0};
			ctmat[0] = del;
			ctmat[1+3*1] = del;
			ctmat[2+3*2] = del;
		} else throw invalid_argument("invalid orbital for calculating charge transfer energy");
	} catch(const exception &ex) {cout << ex.what() << "\n";}
	return ctmat;
}

void calc_CF(Hilbert& hilbs, double del, double tenDQ, int coord) {
	//Calculate Crystal Field, Charge Transfer Matrix Element
	for (int i = 0; i < hilbs.atlist.size(); ++i) {
		int l = hilbs.atlist[i].l;
		double* cfmat;
		if (l == 2 && tenDQ != 0) cfmat = CFmat(l, tenDQ, coord);
		else if (hilbs.atlist[i].is_lig && del != 0) cfmat = CTmat(l, del);
		else continue;
		for (int j = 0; j < (l*2+1)*(l*2+1); ++j) {
			if (cfmat[j] == 0) continue;
			int ml1 = j%(l*2+1)-l, ml2 = j/(l*2+1)-l;
			for (auto spin : {-0.5,0.5}) {
				QN qn1(ml1,spin,i), qn2(ml2,spin,i);
				vpulli entries = hilbs.match(1,&qn1,&qn2);
				for (auto& e : entries) {
					if (e.first == e.second) hilbs.fill_hblk(cfmat[j],e.first,e.second);
					else hilbs.fill_hblk(cfmat[j]*hilbs.Fsign(&qn1,e.first,1)
						 *hilbs.Fsign(&qn2,e.second,1),e.first,e.second);
				}
			}
		}
	}
}

void calc_SO(Hilbert& hilbs, double lambda) {
	// Calculate Spin Orbit Coupling Matrix Element
	if (lambda == 0) return;
	for (int i = 0; i < hilbs.atlist.size(); ++i) {
		int l = hilbs.atlist[i].l;
		if (l != 1 || hilbs.atlist[i].is_lig) continue;// Only 2p orbitals get spin orbit coupling for now
		for (int ml = -l; ml <= l; ++ml) {
			// longitudinal phonon
			for (auto spin : {-0.5,0.5}) {
				QN qn(ml,spin,i);
				vpulli entries = hilbs.match(1,&qn,&qn);
				double matelem = -lambda * spin * ml; // No sign traversing issue due to symmetry
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

void calc_CV(Hilbert& hilbs, double* FG) {
	// Calculate core valence interactions in transition metals
	if (!hilbs.is_ex) return;
	cout << "core valence" << endl;
	for (int ci = 0; ci < hilbs.val_ati; ci++) {
		for (int vi = hilbs.val_ati; vi < hilbs.atlist.size(); vi++) {
			cout << "ci: " << ci << ", vi: " << vi << endl;
			int cl = hilbs.atlist[ci].l, vl = hilbs.atlist[vi].l;
			for (int vml = -vl; vml <= vl; ++vml) {
			for (int vmr = -vl; vmr <= vl; ++vmr) {
			for (int cml = -cl; cml <= cl; ++cml) {
			for (int cmr = -cl; cmr <= cl; ++cmr) {
				if (vml+cml != vmr+cmr) continue;
				vector<pair<double,double>> spairs = {{0.5,0.5},{0.5,-0.5},{-0.5,0.5},{-0.5,-0.5}};
				for (auto & s : spairs) {
					struct QN qnl[2] = {{vml,s.first,vi},{cml,s.second,ci}};
					struct QN qnr[2] = {{vmr,s.first,vi},{cmr,s.second,ci}};
					vpulli entries = hilbs.match(2,qnl,qnr);
					double meF = calc_U(gaunt(vl,vml,vl,vmr),gaunt(cl,cmr,cl,cml),FG,2*cl+1);
					for (auto& e : entries) 
						hilbs.fill_hblk(meF*hilbs.Fsign(qnl,e.first,2)*hilbs.Fsign(qnr,e.second,2),e.first,e.second);
					qnl[0] = {vml,s.second,vi}, qnl[1] = {cml,s.first,ci};
					qnr[0] = {vmr,s.first,vi}, qnr[1] = {cmr,s.second,ci};
					entries = hilbs.match(2,qnl,qnr);
					double meG = calc_U(gaunt(cl,cml,vl,vmr),gaunt(cl,cmr,vl,vml),FG,cl+vl+1);
					for (auto& e : entries)
						hilbs.fill_hblk(-meG*hilbs.Fsign(qnl,e.first,2)*hilbs.Fsign(qnr,e.second,2),e.first,e.second);
				}
			}}}}
			// !!!!!!Shift for particle hole symmetry
		}
	}
	return;
}

void calc_HYB(Hilbert& hilbs, double* SC) {
	// Charge Transfer
	for (int i = 0; i < hilbs.atlist.size(); ++i) {
		if (hilbs.atlist[i].is_lig) {
			cout << "do something" << endl;
		}
	}
	// Hybridization
	return;
}
