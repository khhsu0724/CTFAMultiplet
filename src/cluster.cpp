#include "cluster.hpp"
using namespace std;

// This file sets up the geometry of the cluster, including hybridization & orbitals
// Shared Cluster Implementation
// The matrix here uses gauge choice from PhysRevB.93.155166
void Cluster::make_atlist(std::vector<Atom>& atlist, int num_vh, 
						const std::vector<int>& sites) {
	std::vector<int> nh_per_tm = ed::distribute(num_vh,tm_per_site);
	int atind = 0;
	atlist.reserve(at_per_site()*num_sites);
	for (int x = 0; x < sites[0]; ++x) {
	for (int y = 0; y < sites[1]; ++y) {
	for (int z = 0; z < sites[2]; ++z) {
		std::vector<int> site = {x,y,z};
		for (size_t tm = 0; tm < tm_per_site; ++tm) {
			if (edge == "L") atlist.emplace(atlist.begin(),Atom("2p",atind,3,2,0,site));
			atlist.emplace_back(Atom("3d",atind,3,2,nh_per_tm[tm],site));
			atind++;
		}
		for (size_t lig = 0; lig < lig_per_site; ++lig) {
			if (edge == "K") atlist.emplace(atlist.begin(),Atom("1s",atind,2,1,0,site));
			atlist.emplace_back(Atom("2p",atind,2,1,0,site));
			atind++;
		}
	}}}
	return;
};

void Cluster::print_eigstate(const vecd& occ, bool is_print, string fname, int p) {
	multistream mout(is_print,fname);
	// Print Header
	mout << setw(w) << "Num Holes: ";
	int nvo = vo_persite * num_sites;
	if (occ.size() % nvo != 0) throw runtime_error("Wrong occ vector size");
	for (size_t n = 0; n < occ.size(); n += nvo) {
		double num_h = 0;
		for (int i = n; i < n+nvo; ++i) num_h += occ[i];
		mout << setprecision(3) << setw(w) << num_h;
	}
	mout << endl;
	// Print Occupation
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
		if (num_sites_print > 1) mout << setw(w) << "SITE NUMBER: " << s << endl;
		for (size_t i = s*vo_persite; i < (s+1)*vo_persite; ++i) {
			mout << setw(w) << (orb_names[i-s*vo_persite]+": ") << fixed << setprecision(p);
			for (size_t n = i; n < occ_print.size(); n+=vo_persite*num_sites_print) {
				mout << setw(w) << occ_print[n];
			}
			mout << endl;
		}
	}
	return;
};

vecc Cluster::get_inp_tmat() {
	// Please use first line as key for sparse or dense matrix
	// Sparse: [x,y]: (x,y)
	// Dense: (real,imag),(real,imag)
	// TODO: Read in specified datatype?
	vecc tmat;
	ifstream infile(inp_hyb_file);
	if (!infile.good()) throw invalid_argument(inp_hyb_file+" does not exist!");
	string mat_type, line;
	getline(infile,mat_type);
	transform(mat_type.begin(),mat_type.end(),mat_type.begin(),::toupper);
	if (mat_type != "DENSE" && mat_type != "SPARSE") {
		mat_type = "DENSE";
		infile.seekg(0,std::ios_base::beg);
	}
	if (mat_type == "DENSE") {
		tmat.reserve(vo_persite*vo_persite);
		dcomp matelem = 0;
		while (infile >> line) {
			ed::parse_num(line,matelem);
    		tmat.emplace_back(matelem);
		}
		if (tmat.size() != vo_persite*vo_persite) {
			cout << "Matrix should be size: " << vo_persite << "x" << vo_persite << endl;
			throw invalid_argument("Incorrect input matrix size");
		}
	} else {
		// Reading in sparse matrix
		tmat = vecc(vo_persite*vo_persite,0);
		while (getline(infile,line)) {
			if (line.empty()) continue;
			line = line.substr(0,line.find("#")); // Get rid of comments
			size_t colon =  line.find(":");
			dcomp matelem = 0;
			ed::parse_num(line.substr(colon+1),matelem);
			line = line.substr(line.find("[")+1);
			line = line.substr(0,line.find("]"));
			size_t comma =  line.find(",");
			int x = stoi(line.substr(0,comma));
			int y = stoi(line.substr(comma+1));
			if (x < 0 || x >= vo_persite || y < 0 || y >= vo_persite)
				throw invalid_argument("Invalid matrix entry from " + inp_hyb_file);
			tmat[x+vo_persite*y] += matelem;
		}
	}
	return tmat;
};

vecd Cluster::get_tmat_real() {
	// Swap spherical harmonics basis
	vecc tmat;
	if (inp_hyb_file.empty()) tmat = get_tmat();
	else tmat = get_inp_tmat();
	ed::write_vec(tmat,vo_persite,vo_persite,"./tmat.txt");
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
	orb_names = {"dx2","dz2","dxy","dxz","dyz"};
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
	orb_names = {"dx2","dz2","dxy","dxz","dyz",
				"pxx","pxy","pxz","pyx","pyy","pyz"};
	return;
};

void SquarePlanar::set_hyb_params(const HParam& hparam) {
	tpd = hparam.tpd;
	tpp = hparam.tpp;
	del = hparam.MLdelta;
	sig_pi = hparam.sig_pi;
	tpdz_ratio = hparam.tpdz_ratio;
	tppsigma_on = hparam.tppsigma_on;
	return;
};

vecc SquarePlanar::get_seph2real_mat() {
	return ed::make_blk_mat(U_d(),1,U_p(),2);
};

vecc SquarePlanar::get_operator_mat() {
	vecc a(vo_persite*vo_persite,0);
	for(int i = 0; i < 5; i++) a[i*vo_persite+i] = {1,0};
	for(int i = 5; i < 8; i++) a[i*vo_persite+i] = {0,1};
	for(int i = 8; i < 11; i++) a[i*vo_persite+i] = {1,0};
	return a;
};

vecc SquarePlanar::get_tmat() {
	// TODO: Change the hybridization matrix consistent with Oct definition
	// The ratio of tpd bonds can be tweaked individually
	double tpdz = tpd*tpdz_ratio, tpdxy = tpd*sig_pi, tpdxz = tpd*sig_pi, tpdyz = tpd*sig_pi;
	double tpppi = 0*tpp, tppsigma = tpp, tppzpi = 0*tpp; // Introduces tppzpi
	double kx = PI/2, ky = PI/2, kz = 0; // momentum at (pi,pi) for hole language
	double sx = 2*sin(kx), sy = 2*sin(ky), cx = 2*cos(kx), cy = 2*cos(ky);
	int nvo = vo_persite;
	// if (tppsigma_on) tppzpi = 20*tpp;
	cout << "tpd: " << tpd << ", tdpz: " << tpdz << ", tpdxy: " << tpdxy << ", tpdxz: " << tpdxz << endl;
	cout << "tpp: " << tpp << ", tppzpi: " << tppzpi << endl;
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
	tmat[5*nvo+9] = tmat[9*nvo+5] = {tppsigma*sx*sy,0};
	tmat[6*nvo+8] = tmat[8*nvo+6] = {tppsigma*sx*sy,0};
	tmat[6*nvo+9] = tmat[9*nvo+6] = {-tpppi*cx*cy,0};
	tmat[7*nvo+10] = tmat[10*nvo+7] = {tppzpi*cx*cy,0};
	return tmat;
};
// End of Square-Planar Implementation

// Octahedral Implementation
Octahedral::Octahedral(std::string edge) : Cluster(edge) {
	co_persite = 3;
	vo_persite = 14;
	lig_per_site = 3;
	tm_per_site = 1;
	no_HYB = false;
	orb_names = {"dx2","dz2","dxy","dxz","dyz",
				"pxx","pxy","pxz","pyx","pyy",
				"pyz","pzx","pzy","pzz"};
	return;
};

void Octahedral::set_hyb_params(const HParam& hparam) {
	tpd = hparam.tpd;
	tpp = hparam.tpp;
	del = hparam.MLdelta;
	octJT = hparam.octJT;
	tpd_sig_pi = hparam.sig_pi;
	tppsigma_on = hparam.tppsigma_on;
	return;
};

vecc Octahedral::get_seph2real_mat() {
	return ed::make_blk_mat(U_d(),1,U_p(),3);
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
	double tpdz = tpd/sqrt(3.0), tpdxy = tpd*tpd_sig_pi, tpdxz = tpd*tpd_sig_pi, tpdyz = tpd*tpd_sig_pi;
	// tpp_pi_1 is pxy/pyy type overlaps, tpp_pi_2 is pzy/pxy type overlaps
	double vppsigma = tpp, vpppi = tpp*tpp_sig_pi;
	double tppsigma_1 = 0.5*(vppsigma+vpppi), tppsigma_2 = 0.5*(vppsigma-vpppi), tpppi = vpppi;
	double kx = PI/2, ky = PI/2, kz = PI/2; // relative oxygen positionss
	double sx = 2*sin(kx), sy = 2*sin(ky), sz = 2*sin(kz);
	double cx = 2*cos(kx), cy = 2*cos(ky), cz = 2*sin(kz);
	// double d = octJT;
	double d = 1;
	std::cout << "tpd modulation: " << pow(d,-4) << ", tpp modulation: " << pow(d,-3) << std::endl;
	cout << "tpd: " << tpd << ", tdpz: " << tpdz << ", tpdxy: " << tpdxy << ", tpdxz: " << tpdxz;
	cout << ", tpdyz: " << tpdyz << endl;
	cout << "tpp: " << tpp << ", vppsigma: " << vppsigma << ", vpppi: " << vpppi << endl;
	cout << "tppsigma 1: " << tppsigma_1 << ", tppsigma 2: " << tppsigma_2 << ", tpppi: " << tpppi << endl;
	int nvo = vo_persite;
	vecc tmat(nvo*nvo,0);

	tmat[0*nvo+5] = {0,-tpd*sx};
	tmat[5*nvo+0] = -tmat[0*nvo+5];
	tmat[0*nvo+9] = {0,tpd*sy};
	tmat[9*nvo+0] =  -tmat[0*nvo+9];
	tmat[1*nvo+5] = {0,tpdz*sx};
	tmat[5*nvo+1] = -tmat[1*nvo+5];
	tmat[1*nvo+9] = {0,tpdz*sy};
	tmat[9*nvo+1] =  -tmat[1*nvo+9];
	tmat[1*nvo+13] = {0,-tpd*sz*pow(d,-4)}; // Modulate the strength?	
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
	tmat[5*nvo+8] = tmat[8*nvo+5] = {-tppsigma_1*cx*cy,0};
	tmat[5*nvo+9] = tmat[9*nvo+5] = {tppsigma_2*sx*sy,0};
	tmat[5*nvo+11] = tmat[11*nvo+5] = {-tppsigma_1*cx*cz,0};
	tmat[5*nvo+13] = tmat[13*nvo+5] = {tppsigma_2*sx*sz*pow(d,-3),0};
	tmat[6*nvo+8] = tmat[8*nvo+6] = {tppsigma_2*sx*sy,0};
	tmat[6*nvo+9] = tmat[9*nvo+6] = {-tppsigma_1*cx*cy,0};
	tmat[6*nvo+12] = tmat[12*nvo+6] = {tpppi*cx*cz,0};
	tmat[7*nvo+10] = tmat[10*nvo+7] = {tpppi*cx*cy,0};
	tmat[7*nvo+11] = tmat[11*nvo+7] = {tppsigma_2*sx*sz*pow(d,-3),0};
	tmat[7*nvo+13] = tmat[13*nvo+7] = {-tppsigma_1*cx*cz,0};
	tmat[8*nvo+8] = tmat[9*nvo+9] = tmat[10*nvo+10] = del;
	tmat[8*nvo+11] = tmat[11*nvo+8] = {tpppi*cy*cz,0};
	tmat[9*nvo+12] = tmat[12*nvo+9] = {-tppsigma_1*cy*cz,0};
	tmat[9*nvo+13] = tmat[13*nvo+9] = {tppsigma_2*sy*sz*pow(d,-3),0};
	tmat[10*nvo+12] = tmat[12*nvo+10] = {tppsigma_2*sy*sz*pow(d,-3),0};
	tmat[10*nvo+13] = tmat[13*nvo+10] = {-tppsigma_1*cy*cz,0};
	tmat[11*nvo+11] = tmat[12*nvo+12] = tmat[13*nvo+13] = del;
	return tmat;
};

// End of Octahedral Implementation

// ==========================================
// Tetrahedral Implementation
// ==========================================
Tetrahedral::Tetrahedral(std::string edge) : Cluster(edge) {
    co_persite = 3;
    vo_persite = 17; // 5 d-orbitals + (4 ligands * 3 p-orbitals)
    lig_per_site = 4;
    tm_per_site = 1;
    no_HYB = false;
	orb_names = {"dx2","dz2","dxy","dxz","dyz",
					"p0x","p0y","p0z",
					"p1x","p1y","p1z",
					"p2x","p2y","p2z",
					"p3x","p3y","p3z"};
    return;
};

void Tetrahedral::set_hyb_params(const HParam& hparam) {
    tpd = hparam.tpd;
    tpp = hparam.tpp;
    del = hparam.MLdelta;
    octJT = hparam.octJT;
    tpd_sig_pi = hparam.sig_pi;
    tppsigma_on = hparam.tppsigma_on;
    
    // NEW: Grab the 6th CF parameter for the distortion
    D_trig = hparam.CF[5]; 
    return;
};

vecc Tetrahedral::get_seph2real_mat() {
    // 1 d-site, 4 p-sites
    return ed::make_blk_mat(U_d(), 1, U_p(), 4);
};

vecc Tetrahedral::get_operator_mat() {
    vecc a(vo_persite*vo_persite,0);
    for(int i = 0; i < 5; i++) a[i*vo_persite+i] = {1,0};   // d-orbitals
    for(int i = 5; i < 8; i++) a[i*vo_persite+i] = {0,1};   // L0
    for(int i = 8; i < 11; i++) a[i*vo_persite+i] = {1,0};  // L1
    for(int i = 11; i < 14; i++) a[i*vo_persite+i] = {0,1}; // L2
    for(int i = 14; i < 17; i++) a[i*vo_persite+i] = {1,0}; // L3
    return a;
};

vecc Tetrahedral::get_tmat() {
    double V_sigma = tpd;
    double V_pi = tpd * tpd_sig_pi; 

    int nvo = vo_persite;
    vecc tmat(nvo*nvo, 0);

    // 1. Set the Ligand charge transfer energy (diagonal elements)
    for (int i = 5; i < 17; ++i) {
        tmat[i*nvo+i] = {del, 0};
    }

    // 2. Define the exact directional cosines for a perfect tetrahedron
    double rt3 = sqrt(3.0);
    double inv_rt3 = 1.0 / rt3;

    double ligands[4][3] = {
        { inv_rt3,  inv_rt3,  inv_rt3}, // L0: [ 1,  1,  1]
        {-inv_rt3, -inv_rt3,  inv_rt3}, // L1: [-1, -1,  1]
        { inv_rt3, -inv_rt3, -inv_rt3}, // L2: [ 1, -1, -1]
        {-inv_rt3,  inv_rt3, -inv_rt3}  // L3: [-1,  1, -1]
    };

	int d_x2y2 = 0;
	int d_z2   = 1;
	int d_xy   = 2;
	int d_xz   = 3;
	int d_yz   = 4;

	// NEW: Apply the trigonal distortion to the t2 manifold
	if (D_trig != 0.0) {
		// dxy - dxz coupling
		tmat[d_xy*nvo + d_xz] = {D_trig, 0};
		tmat[d_xz*nvo + d_xy] = {D_trig, 0};
		
		// dxy - dyz coupling
		tmat[d_xy*nvo + d_yz] = {D_trig, 0};
		tmat[d_yz*nvo + d_xy] = {D_trig, 0};
		
		// dxz - dyz coupling
		tmat[d_xz*nvo + d_yz] = {D_trig, 0};
		tmat[d_yz*nvo + d_xz] = {D_trig, 0};
	}

    // 3. Build the tight-binding hoppings using Slater-Koster rules
    for (int L = 0; L < 4; ++L) {
        double l = ligands[L][0];
        double m = ligands[L][1];
        double n = ligands[L][2];

        double l2 = l*l, m2 = m*m, n2 = n*n;

        // Orbital indices (Mapping to the orb_names array)
        int p_x = 5 + 3*L;
        int p_y = 6 + 3*L;
        int p_z = 7 + 3*L;

        // Calculate hoppings from p to d
        double h_px_dx2 = 0.5 * rt3 * l * (l2 - m2) * V_sigma + l * (1.0 - l2 + m2) * V_pi;
        double h_py_dx2 = 0.5 * rt3 * m * (l2 - m2) * V_sigma - m * (1.0 + l2 - m2) * V_pi;
        double h_pz_dx2 = 0.5 * rt3 * n * (l2 - m2) * V_sigma - n * (l2 - m2) * V_pi;

        double h_px_dz2 = l * (n2 - 0.5 * (l2 + m2)) * V_sigma - rt3 * l * n2 * V_pi;
        double h_py_dz2 = m * (n2 - 0.5 * (l2 + m2)) * V_sigma - rt3 * m * n2 * V_pi;
        double h_pz_dz2 = n * (n2 - 0.5 * (l2 + m2)) * V_sigma + rt3 * n * (l2 + m2) * V_pi;

        double h_px_dxy = rt3 * l2 * m * V_sigma + m * (1.0 - 2.0 * l2) * V_pi;
        double h_py_dxy = rt3 * l * m2 * V_sigma + l * (1.0 - 2.0 * m2) * V_pi;
        double h_pz_dxy = rt3 * l * m * n * V_sigma - 2.0 * l * m * n * V_pi;

        double h_px_dxz = rt3 * l2 * n * V_sigma + n * (1.0 - 2.0 * l2) * V_pi;
        double h_py_dxz = rt3 * l * m * n * V_sigma - 2.0 * l * m * n * V_pi;
        double h_pz_dxz = rt3 * l * n2 * V_sigma + l * (1.0 - 2.0 * n2) * V_pi;

        double h_px_dyz = rt3 * l * m * n * V_sigma - 2.0 * l * m * n * V_pi;
        double h_py_dyz = rt3 * m2 * n * V_sigma + n * (1.0 - 2.0 * m2) * V_pi;
        double h_pz_dyz = rt3 * m * n2 * V_sigma + m * (1.0 - 2.0 * n2) * V_pi;

        // Populate the Hermitian matrix. We use purely real components {val, 0}
        tmat[d_x2y2*nvo + p_x] = {h_px_dx2, 0}; tmat[p_x*nvo + d_x2y2] = {h_px_dx2, 0};
        tmat[d_x2y2*nvo + p_y] = {h_py_dx2, 0}; tmat[p_y*nvo + d_x2y2] = {h_py_dx2, 0};
        tmat[d_x2y2*nvo + p_z] = {h_pz_dx2, 0}; tmat[p_z*nvo + d_x2y2] = {h_pz_dx2, 0};

        tmat[d_z2*nvo + p_x] = {h_px_dz2, 0}; tmat[p_x*nvo + d_z2] = {h_px_dz2, 0};
        tmat[d_z2*nvo + p_y] = {h_py_dz2, 0}; tmat[p_y*nvo + d_z2] = {h_py_dz2, 0};
        tmat[d_z2*nvo + p_z] = {h_pz_dz2, 0}; tmat[p_z*nvo + d_z2] = {h_pz_dz2, 0};

        tmat[d_xy*nvo + p_x] = {h_px_dxy, 0}; tmat[p_x*nvo + d_xy] = {h_px_dxy, 0};
        tmat[d_xy*nvo + p_y] = {h_py_dxy, 0}; tmat[p_y*nvo + d_xy] = {h_py_dxy, 0};
        tmat[d_xy*nvo + p_z] = {h_pz_dxy, 0}; tmat[p_z*nvo + d_xy] = {h_pz_dxy, 0};

        tmat[d_xz*nvo + p_x] = {h_px_dxz, 0}; tmat[p_x*nvo + d_xz] = {h_px_dxz, 0};
        tmat[d_xz*nvo + p_y] = {h_py_dxz, 0}; tmat[p_y*nvo + d_xz] = {h_py_dxz, 0};
        tmat[d_xz*nvo + p_z] = {h_pz_dxz, 0}; tmat[p_z*nvo + d_xz] = {h_pz_dxz, 0};

        tmat[d_yz*nvo + p_x] = {h_px_dyz, 0}; tmat[p_x*nvo + d_yz] = {h_px_dyz, 0};
        tmat[d_yz*nvo + p_y] = {h_py_dyz, 0}; tmat[p_y*nvo + d_yz] = {h_py_dyz, 0};
        tmat[d_yz*nvo + p_z] = {h_pz_dyz, 0}; tmat[p_z*nvo + d_yz] = {h_pz_dyz, 0};
    }

    // We only to d-p hoppings, p-p can be added later...

    return tmat;
};

// End of Tetrahedral implementaton