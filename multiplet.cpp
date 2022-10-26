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
	calc_coulomb(hilbs,hparam.SC); 
	if (hilbs.SO_on) calc_SO(hilbs,hparam.SO);
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
		if (l == 2 && CF != 0 && !hilbs.atlist[i].is_lig) cfmat = CFmat(l, CF);
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

void calc_SO(Hilbert& hilbs, const double lambda) {
	// Calculate Spin Orbit Coupling Matrix Element
	if (lambda == 0) return;
	for (int i = 0; i < hilbs.atlist.size(); ++i) {
		int l = hilbs.atlist[i].l;
		if (l != 1 || hilbs.atlist[i].is_lig) continue;// Only 2p core orbitals get spin orbit coupling for now
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

void calc_CV(Hilbert& hilbs, const double* FG) {
	// Calculate core valence interactions
	if (!hilbs.is_ex) return;
	// int ci = 0, vi = 1;
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
		// !!!!!!Shift for particle hole symmetry
	}
	return;
}

vecd HYBmat(Hilbert& hilbs, const HParam& hparam) {
	// Calculate Hybridization, charge transfer and crystal field energy
	// !!!! Need to factor in different sites...how to do this?
	// Can we identify what type of geomtery algorithmically
	string coord = hilbs.coord;
	// cout << "building hybdrization matrix" << endl;
	int nvo = hilbs.num_vorb/2, nco = hilbs.num_corb/2;
	vecd hybmat(nvo*nvo,0),tmatreal(nvo*nvo,0);
	double t = hparam.HYB, del = hparam.MLdelta;
	double tpd = t, tpdz = t*0.45, tpdxy = t*0.45, tpdxz = t*0.45, tpdyz = t*0.45, tpppi = 0, tppsigma = t/3, tppzpi = 0;
	// double tpd = 1, tpdz = 0.25, tpdxy = 0.225, tpdxz = 0.225, tpdyz = 0.225, tpppi = 0, tppsigma = 0.25, tppzpi = 0; // Temp place holdler
	// Temporary implementation
	// Omit tpppi since it will give imaginary number
	double kx = PI, ky = PI, kz = 0; // one site momentum = 0, p-h transformation k+Q = pi
	double t0 = tpdz/sqrt(2), t1x = tpdxz/sqrt(2), t1y = tpdyz/sqrt(2), t2p = 0.5*(tpd+tpdxy), t2m = 0.5*(tpd-tpdxy);
	double sx = 2*sin(kx/2), sy = 2*sin(ky/2), cx = 2*cos(kx/2), cy = 2*cos(ky/2), tpsigma = tppsigma; // = 2*tppsigma???
	cx = cy = 2;
	hybmat[5*nvo+0] = hybmat[0*nvo+5] = t2m*cx;//-0.5*(tpd+tpdxy)*2;
	hybmat[7*nvo+0] = hybmat[0*nvo+7] = t2p*cx;//-0.5*(tpd-tpdxy)*2;
	hybmat[8*nvo+0] = hybmat[0*nvo+8] = -t2m*cy;//-0.5*(-tpd+tpdxy)*2;
	hybmat[10*nvo+0] = hybmat[0*nvo+10] = t2p*cy;//-0.5*(-tpd-tpdxy)*2;
	hybmat[6*nvo+1] = hybmat[1*nvo+6] = -t1x*cx;//tpdxz/sqrt(2)*2;
	hybmat[9*nvo+1] = hybmat[1*nvo+9] = t1y*cy;//tpdyz/sqrt(2)*2;
	hybmat[5*nvo+2] = hybmat[2*nvo+5] = t0*cx;//-tpdz/sqrt(2)*2;
	hybmat[7*nvo+2] = hybmat[2*nvo+7] = t0*cx;//-tpdz/sqrt(2)*2;
	hybmat[8*nvo+2] = hybmat[2*nvo+8] = t0*cy;// -tpdz/sqrt(2)*2;
	hybmat[10*nvo+2] = hybmat[2*nvo+10] = -t0*cy; //-tpdz/sqrt(2)*2;
	hybmat[6*nvo+3] = hybmat[3*nvo+6] = t1x*cx; //-tpdxz/sqrt(2)*2;
	hybmat[9*nvo+3] = hybmat[3*nvo+9] = t1y*cy; //-tpdyz/sqrt(2)*2;
	hybmat[5*nvo+4] = hybmat[4*nvo+5] = t2p*cx; //-0.5*(tpd-tpdxy)*2;
	hybmat[7*nvo+4] = hybmat[4*nvo+7] = t2m*cx; //-0.5*(tpd+tpdxy)*2;
	hybmat[8*nvo+4] = hybmat[4*nvo+8] = -t2p*cy; //-0.5*(-tpd-tpdxy)*2;
	hybmat[10*nvo+4] = hybmat[4*nvo+10] = t2m*cy; //-0.5*(-tpd+tpdxy)*2;
	hybmat[10*nvo+5] = hybmat[5*nvo+10] = tpsigma*cx*cy;
	hybmat[8*nvo+7] = hybmat[7*nvo+8] = -tpsigma*cx*cy;
	hybmat[5*nvo+5] = hybmat[6*nvo+6] = hybmat[7*nvo+7] = hybmat[8*nvo+8] = hybmat[9*nvo+9] = hybmat[10*nvo+10] = del;
	
	//
	// Temporary implementation
	// tpd = 1.0, tpdz = 1, tpdxy = 0.5, tpdxz = 1, tpdyz = 1, tpppi = 0.1, tppsigma = 1, tppzpi = 10;
	// sx = 2, sy = 2, cx = 0, cy = 0, tpsigma = tppsigma; // = 2*tppsigma???
	vecc tmat(nvo*nvo,0);
	tmat[0*nvo+5] = {0,-tpd*sx};
	tmat[5*nvo+0] = -tmat[0*nvo+5];
	tmat[0*nvo+9] = {0,tpd*sy};
	tmat[9*nvo+0] =  -tmat[0*nvo+9];
	tmat[1*nvo+5] = {0,-tpdz*sx};
	tmat[5*nvo+1] = -tmat[1*nvo+5];
	tmat[1*nvo+9] = {0,-tpdz*sy};
	tmat[9*nvo+1] =  -tmat[1*nvo+9];
	tmat[2*nvo+6] = {0,tpdxy*sx};
	tmat[6*nvo+2] = -tmat[2*nvo+6];
	tmat[2*nvo+8] = {0,tpdxy*sy};
	tmat[8*nvo+2] =  -tmat[2*nvo+8];
	tmat[3*nvo+7] = {0,tpdxz*sx};
	tmat[7*nvo+3] =  -tmat[3*nvo+7];
	tmat[4*nvo+10] = {0,tpdyz*sy};
	tmat[10*nvo+4] =  -tmat[4*nvo+10];
	tmat[5*nvo+5] = tmat[6*nvo+6] = tmat[7*nvo+7] = tmat[8*nvo+8] = tmat[9*nvo+9] = tmat[10*nvo+10] = del;
	tmat[5*nvo+8] = tmat[8*nvo+5] = {-tpppi*cx*cy,0};
	tmat[5*nvo+9] = tmat[9*nvo+5] = {-tppsigma*sx*sy,0};
	tmat[6*nvo+8] = tmat[8*nvo+6] = {-tppsigma*sx*sy,0};
	tmat[6*nvo+9] = tmat[9*nvo+6] = {-tpppi*cx*cy,0};
	tmat[7*nvo+10] = tmat[10*nvo+7] = {tppzpi*cx*cy,0};

	vecc U(nvo*nvo,0);
	U[0*nvo+0] = U[0*nvo+4] = {1/sqrt(2),0};
	U[1*nvo+2] = {1,0};
	U[2*nvo+0] = {0,1/sqrt(2)};
	U[2*nvo+4] = {0,-1/sqrt(2)};
	U[3*nvo+1] = {1/sqrt(2),0};
	U[3*nvo+3] = {-1/sqrt(2),0};
	U[4*nvo+1] = U[4*nvo+3] = {0,1/sqrt(2)};
	U[5*nvo+5] = U[5*nvo+7] = U[8*nvo+8] = U[8*nvo+10] = {1/sqrt(2),0};
	U[6*nvo+5] = U[9*nvo+8] = {0,1/sqrt(2)};
	U[6*nvo+7] = U[9*nvo+10] = {0,-1/sqrt(2)};
	U[7*nvo+6] = U[10*nvo+9] = {1,0};

	// std::transform(hybmat.begin(), hybmat.end(), tmat.begin(), [](double d) {return complex<double>(d,0);});
	tmat = ed::matmult(tmat,U,nvo);
	ed::ctranspose(U,nvo,nvo);
	tmat = ed::matmult(U,tmat,nvo);

	// Eliminate i in matrix
	vecc a(nvo*nvo,0);
	for(int i = 0; i < 5; i++) a[i*nvo+i] = {1,0};
	for(int i = 5; i < 8; i++) a[i*nvo+i] = {0,1};
	for(int i = 8; i < 11; i++) a[i*nvo+i] = {1,0};
	tmat = ed::matmult(tmat,a,nvo);
	ed::ctranspose(a,nvo,nvo);
	tmat = ed::matmult(a,tmat,nvo);
	// ed::write_vec(tmat,11,11,"hyb.txt");
	// vecd tmatd(nvo*nvo);
	std::transform(tmat.begin(), tmat.end(), tmatreal.begin(), [](complex<double> c) {return c.real();});
	// Particle hole transformation???????
	return tmatreal;
}	

void calc_HYB(Hilbert& hilbs, const HParam& hparam) {
	// Charge Transfer and Hybridization
	// Temperory implementation for square planar geometry
	if (hilbs.at_per_site == 1) return;
	int nvo = hilbs.num_vorb/2, nco = hilbs.num_corb/2;
	double nd = 0;
	for (auto &at : hilbs.atlist) if(!at.is_lig && at.is_val) nd = at.num_h;
	double delta = nd * (hparam.SC[1][0] - hparam.SC[1][2]*2/63 - hparam.SC[1][4]*2/63);// This needs to be changed since 
	int nh = hilbs.is_ex ? hilbs.num_vh+1 : hilbs.num_vh;
	double U = (hparam.SC[1][0] - hparam.SC[1][2]*2/63 - hparam.SC[1][4]*2/63);
	// delta = 0.7*nh*U-0.3;
	delta = nh*U;
	delta = 19.28;

	cout << "U: " << (hparam.SC[1][0]-hparam.SC[1][2]*2/63-hparam.SC[1][4]*2/63);
	cout << ", JH: " << (hparam.SC[1][2]/14+hparam.SC[1][4]/14);
	cout << ", holes: " << nh << ", delta: " << hparam.MLdelta << ", t: " << hparam.HYB << endl;
	vecd hybmat = HYBmat(hilbs,hparam);
	// Get Hybridization Information, there should be num_orb x num_orb matrix providing hybdrization information
	// Loop through hyb matrix, TODO: Loop through each site
	for (int i = 0; i < nvo; ++i) {
		for (int j = 0; j < nvo; ++j) {
			if (abs(hybmat[i*nvo+j]) < TOL) continue;
			// Fill in spin down matrix element
			ulli lhssd = 1 << i, rhssd = 1 << j, incsd = (rhssd|lhssd);
			ulli lhssd_shift = (incsd-lhssd) << nco, rhssd_shift = (incsd-rhssd) << nco;	
			ulli lhssu = 1 << (i+nvo), rhssu = 1 << (j+nvo), incsu = (rhssu|lhssu);
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
