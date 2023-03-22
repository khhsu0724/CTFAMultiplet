#include "cluster.hpp"
using namespace std;

// Shared Cluster Implementation
// void Cluster::print_state_orbs(const vecd& occ, const vector<string>& names, int p) {
// 	for (size_t i = 0; i < vo_persite; ++i) {
// 		cout << setw(w) << (names[i]+": ") << fixed << setprecision(p);
// 		for (size_t n = i; n < occ.size(); n+=vo_persite) cout << setw(w) << occ[n];
// 		cout << endl;
// 	}
// }

void Cluster::print_state_header(const vecd& occ) {
	cout << setw(w) << "Num Holes: ";
	int nvo = vo_persite * num_sites;
	if (occ.size() % nvo != 0) throw runtime_error("Wrong occ vector size");
	for (size_t n = 0; n < occ.size(); n += nvo) {
		double num_h = 0;
		for (int i = n; i < n+nvo; ++i) num_h += occ[i];
		cout << setprecision(3) << setw(w) << num_h;
	}
	cout << endl;
	return;
};

void Cluster::print_state_orbs(const vecd& occ, const vector<string>& names, int p) {
	int num_sites_print = num_sites;
	vecd occ_print = occ;
	if (num_sites > 1 && !print_all_sites) {
		occ_print = vecd(occ.size()/num_sites,0);
		num_sites_print = 1;
		for (size_t i = 0; i < occ.size()/num_sites/vo_persite; ++i) {
		for (size_t o = 0; o < occ_print.size()/num_sites; ++o) {
		for (size_t s = 0; s < num_sites; ++s) {
			occ_print[o+i*occ_print.size()/num_sites] += 
				occ[o+vo_persite*s+i*occ.size()/num_sites];
		}}}
	}
	for (size_t s = 0; s < num_sites_print; ++s) {
		if (num_sites_print > 1) cout << setw(w) << "SITE NUMBER: " << s << endl;
		for (size_t i = s*vo_persite; i < (s+1)*vo_persite; ++i) {
			cout << setw(w) << (names[i-s*vo_persite]+": ") << fixed << setprecision(p);
			for (size_t n = i; n < occ_print.size(); n+=vo_persite*num_sites_print) {
				cout << setw(w) << occ_print[n];
			}
			cout << endl;
		}
	}
	return;
}

vecd Cluster::get_tmat_real() {
	// Swap spherical harmonics basis
	vecc tmat = get_tmat();
	vecc U = get_seph2real_mat();
    tmat = ed::matmult(tmat,U,vo_persite);
	U = ed::ctranspose(U,vo_persite,vo_persite);
	tmat = ed::matmult(U,tmat,vo_persite);
	// Swap operator basis
	vecc a = get_operator_mat();
	tmat = ed::matmult(tmat,a,vo_persite);
	a = ed::ctranspose(a,vo_persite,vo_persite);
	tmat = ed::matmult(a,tmat,vo_persite);
	if (!check_tmat_all_real(tmat))
		cerr << "WARNING: hybridization matrix contains imaginary elemenets" << endl;

	vecd tmatreal(vo_persite*vo_persite,0);
	std::transform(tmat.begin(), tmat.end(), tmatreal.begin(), 
						[](complex<double> c) {return c.real();});
	return tmatreal;
};

bool Cluster::check_tmat_all_real(const vecc& tmat) {
	for (auto& t : tmat) {
		if (abs(t.imag()) > TOL) return false;
	}
	return true;
};
// End of Shared Cluster Implementation

// Ion Implementation
Ion::Ion(std::string edge) : Cluster(edge) {
	if (edge == "K") throw std::invalid_argument("TM K edge unavailable");
	vo_persite = 5;
	co_persite = 3;
	lig_per_site = 0;
	tm_per_site = 1;
	no_HYB = true;
	return;
};

void Ion::print_eigstate(const vecd& occ, int p) {
	print_state_header(occ);
	vector<string> names = {"dx2","dz2","dxy","dxz","dyz"};
	print_state_orbs(occ,names,p);
	return;
};
// End of Ion Implementation

// Square-Planar Implementation
SquarePlanar::SquarePlanar(std::string edge) : Cluster(edge) {
	if (edge == "K") co_persite = 2;
	if (edge == "L") co_persite = 3;
	vo_persite = 11;
	lig_per_site = 2;
	tm_per_site = 1;
	no_HYB = false;
	return;
};

void SquarePlanar::set_hyb_params(const HParam& hparam) {
	tpd = hparam.tpd;
	tpp = hparam.tpp;
	del = hparam.MLdelta;
	return;
};

vecc SquarePlanar::get_seph2real_mat() {
	vecc U(vo_persite*vo_persite,0);
	U[0*vo_persite+0] = U[0*vo_persite+4] = {1/sqrt(2),0};
	U[1*vo_persite+2] = {1,0};
	U[2*vo_persite+0] = {0,1/sqrt(2)};
	U[2*vo_persite+4] = {0,-1/sqrt(2)};
	U[3*vo_persite+1] = {1/sqrt(2),0};
	U[3*vo_persite+3] = {-1/sqrt(2),0};
	U[4*vo_persite+1] = U[4*vo_persite+3] = {0,1/sqrt(2)};
	U[5*vo_persite+5] = U[5*vo_persite+7] = U[8*vo_persite+8] = U[8*vo_persite+10] = {1/sqrt(2),0};
	U[6*vo_persite+5] = U[9*vo_persite+8] = {0,1/sqrt(2)};
	U[6*vo_persite+7] = U[9*vo_persite+10] = {0,-1/sqrt(2)};
	U[7*vo_persite+6] = U[10*vo_persite+9] = {1,0};
	return U;
};

vecc SquarePlanar::get_operator_mat() {
	vecc a(vo_persite*vo_persite,0);
	for(int i = 0; i < 5; i++) a[i*vo_persite+i] = {1,0};
	for(int i = 5; i < 8; i++) a[i*vo_persite+i] = {0,1};
	for(int i = 8; i < 11; i++) a[i*vo_persite+i] = {1,0};
	return a;
};

vecc SquarePlanar::get_tmat() {
	// The ratio of tpd bonds can be tweaked individually
	double tpdz = 0.25*tpd, tpdxy = tpd*sig_pi/4, tpdxz = tpd*sig_pi, tpdyz = tpd*sig_pi;
	double tpppi = 0*tpp, tppsigma = tpp, tppzpi = 0*tpp;
	double kx = PI/2, ky = PI/2, kz = 0; // relative oxygen positionss
	double sx = 2*sin(kx), sy = 2*sin(ky), cx = 2*cos(kx), cy = 2*cos(ky);
	int nvo = vo_persite;
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
	tmat[5*nvo+5] = tmat[6*nvo+6] = tmat[7*nvo+7] = del;
	tmat[8*nvo+8] = tmat[9*nvo+9] = tmat[10*nvo+10] = del;
	tmat[5*nvo+8] = tmat[8*nvo+5] = {-tpppi*cx*cy,0};
	tmat[5*nvo+9] = tmat[9*nvo+5] = {-tppsigma*sx*sy,0};
	tmat[6*nvo+8] = tmat[8*nvo+6] = {-tppsigma*sx*sy,0};
	tmat[6*nvo+9] = tmat[9*nvo+6] = {-tpppi*cx*cy,0};
	tmat[7*nvo+10] = tmat[10*nvo+7] = {-tppzpi*cx*cy,0};
	return tmat;
};

void SquarePlanar::print_eigstate(const vecd& occ, int p) {
	print_state_header(occ);
	vector<string> names = {"dx2","dz2","dxy","dxz","dyz",
							"pxx","pxy","pxz","pyx","pyy","pyz"};
	print_state_orbs(occ,names,p);
	return;
};
// End of Square-Planar Implementation

// Octahedral Implementation
Octahedral::Octahedral(std::string edge) : Cluster(edge) {
	co_persite = 3;
	vo_persite = 14;
	lig_per_site = 3;
	tm_per_site = 1;
	no_HYB = false;
	return;
};

void Octahedral::set_hyb_params(const HParam& hparam) {
	tpd = hparam.tpd;
	tpp = hparam.tpp;
	del = hparam.MLdelta;
	return;
};

vecc Octahedral::get_seph2real_mat() {
	vecc U(vo_persite*vo_persite,0);
	U[0*vo_persite+0] = U[0*vo_persite+4] = {1/sqrt(2),0};
	U[1*vo_persite+2] = {1,0};
	U[2*vo_persite+0] = {0,1/sqrt(2)};
	U[2*vo_persite+4] = {0,-1/sqrt(2)};
	U[3*vo_persite+1] = {1/sqrt(2),0};
	U[3*vo_persite+3] = {-1/sqrt(2),0};
	U[4*vo_persite+1] = U[4*vo_persite+3] = {0,1/sqrt(2)};
	U[5*vo_persite+5] = U[5*vo_persite+7] = {1/sqrt(2),0};
	U[8*vo_persite+8] = U[8*vo_persite+10] = {1/sqrt(2),0};
	U[11*vo_persite+11] = U[11*vo_persite+13] = {1/sqrt(2),0};
	U[6*vo_persite+5] = U[9*vo_persite+8] = U[12*vo_persite+11] = {0,1/sqrt(2)};
	U[6*vo_persite+7] = U[9*vo_persite+10] = U[12*vo_persite+13] = {0,-1/sqrt(2)};
	U[7*vo_persite+6] = U[10*vo_persite+9] = U[13*vo_persite+12] = {1,0};
	return U;
};

vecc Octahedral::get_operator_mat() {
	vecc a(vo_persite*vo_persite,0);
	for(int i = 0; i < 5; i++) a[i*vo_persite+i] = {1,0};
	for(int i = 5; i < 8; i++) a[i*vo_persite+i] = {0,1};
	for(int i = 8; i < 11; i++) a[i*vo_persite+i] = {1,0};
	for(int i = 11; i < 14; i++) a[i*vo_persite+i] = {0,1};
	return a;
};

vecc Octahedral::get_tmat() {
	double tpdz = 0.25*tpd, tpdxy = tpd*tpd_sig_pi, tpdxz = tpd*tpd_sig_pi, tpdyz = tpd*tpd_sig_pi;
	// tpp_pi_1 is pxy/pyy type overlaps, tpp_pi_2 is pzy/pxy type overlaps
	double tpppi_1 = tpp_sig_pi_1*tpp, tppsigma = tpp, tpppi_2 = tpp_sig_pi_2*tpp;
	double kx = PI/2, ky = PI/2, kz = PI/2; // relative oxygen positionss
	double sx = 2*sin(kx), sy = 2*sin(ky), sz = 2*sin(kz);
	double cx = 2*cos(kx), cy = 2*cos(ky), cz = 2*sin(kz);
	double d = 1;
	std::cout << "tpd modulation: " << pow(d,-4) << ", tpp modulation: " << pow(d,-3) << std::endl;
	int nvo = vo_persite;
	vecc tmat(nvo*nvo,0);
	tmat[0*nvo+5] = {0,-tpd*sx};
	tmat[5*nvo+0] = -tmat[0*nvo+5];
	tmat[0*nvo+9] = {0,tpd*sy};
	tmat[9*nvo+0] =  -tmat[0*nvo+9];
	tmat[1*nvo+5] = {0,-tpdz*sx};
	tmat[5*nvo+1] = -tmat[1*nvo+5];
	tmat[1*nvo+9] = {0,-tpdz*sy};
	tmat[9*nvo+1] =  -tmat[1*nvo+9];
	tmat[1*nvo+13] = {0,tpd*sz*pow(d,-4)}; // Modulate the strength?	
	tmat[13*nvo+1] =  -tmat[1*nvo+13];
	tmat[2*nvo+6] = {0,tpdxy*sx};
	tmat[6*nvo+2] = -tmat[2*nvo+6];
	tmat[2*nvo+8] = {0,tpdxy*sy};
	tmat[8*nvo+2] =  -tmat[2*nvo+8];
	tmat[3*nvo+7] = {0,tpdxz*sx};
	tmat[7*nvo+3] =  -tmat[3*nvo+7];
	tmat[3*nvo+11] = {0,tpdxz*sz*pow(d,-4)};
	tmat[11*nvo+3] =  -tmat[3*nvo+11];
	tmat[4*nvo+10] = {0,tpdyz*sy};
	tmat[10*nvo+4] =  -tmat[4*nvo+10];
	tmat[4*nvo+12] = {0,tpdyz*sz*pow(d,-4)};
	tmat[12*nvo+4] =  -tmat[4*nvo+12];
	tmat[5*nvo+5] = tmat[6*nvo+6] = tmat[7*nvo+7] = del;
	tmat[5*nvo+8] = tmat[8*nvo+5] = {-tpppi_1*cx*cy,0};
	tmat[5*nvo+9] = tmat[9*nvo+5] = {-tppsigma*sx*sy,0};
	tmat[5*nvo+11] = tmat[11*nvo+5] = {-tpppi_1*cx*cz,0};
	tmat[5*nvo+13] = tmat[13*nvo+5] = {-tppsigma*sx*sz*pow(d,-3),0};
	tmat[6*nvo+8] = tmat[8*nvo+6] = {-tppsigma*sx*sy,0};
	tmat[6*nvo+9] = tmat[9*nvo+6] = {-tpppi_1*cx*cy,0};
	tmat[6*nvo+12] = tmat[12*nvo+6] = {-tpppi_2*cx*cz,0};
	tmat[7*nvo+10] = tmat[10*nvo+7] = {-tpppi_2*cx*cy,0};
	tmat[7*nvo+11] = tmat[11*nvo+7] = {-tppsigma*sx*sz*pow(d,-3),0};
	tmat[7*nvo+13] = tmat[13*nvo+7] = {-tpppi_1*cx*cz,0};
	tmat[8*nvo+8] = tmat[9*nvo+9] = tmat[10*nvo+10] = del;
	tmat[8*nvo+11] = tmat[11*nvo+8] = {-tpppi_2*cy*cz,0};
	tmat[9*nvo+12] = tmat[12*nvo+9] = {-tpppi_1*cy*cz,0};
	tmat[9*nvo+13] = tmat[13*nvo+9] = {-tppsigma*sy*sz*pow(d,-3),0};
	tmat[10*nvo+12] = tmat[12*nvo+10] = {-tppsigma*sy*sz*pow(d,-3),0};
	tmat[10*nvo+13] = tmat[13*nvo+10] = {-tpppi_1*cy*cz,0};
	tmat[11*nvo+11] = tmat[12*nvo+12] = tmat[13*nvo+13] = del;
	return tmat;
};

void Octahedral::print_eigstate(const vecd& occ, int p) {
	print_state_header(occ);
	vector<string> names = {"dx2","dz2","dxy","dxz","dyz",
							"pxx","pxy","pxz","pyx","pyy","pyz",
							"pzx","pzy","pzz"};
	print_state_orbs(occ,names,p);
	return;
};
// End of Octahedral Implementation