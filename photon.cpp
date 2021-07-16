#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <math.h>
#include <complex>
#include <algorithm>
#include "diagonalize.h"
#include "gaunt.h"
#include "multiplet.h"
#include "hilbert.h"
#include "site.h"
#include "photon.h"

using namespace std;
typedef complex<double> dcomp;

// Contains photon based spectroscopy

bool is_pw2(int x) {return !(x == 0) && !(x & (x - 1));}

vector<double> printDistinct(double arr[], int n, bool is_print) {
    // Pick all elements one by one
    vector<double> unique_eig;
    for (int i=0; i<n; i++)
    {
        // Check if the picked element is already printed
        int j;
        for (j=0; j<i; j++)
           if (abs(arr[i] - arr[j]) < 1e-7)
               break;
        // If not printed earlier, then print it
        if (i == j) {
        	if (is_print) cout << arr[i] << " ";
        	unique_eig.push_back(arr[i]);
        }
    }
    return unique_eig;
}

double* get_eigvec(int index, int size, double* eigvec_mat) {
	double* eigvec = new double[size];
	for (int i = 0; i < size; ++i) {
		eigvec[i] = eigvec_mat[index*size+i];
	}
	return eigvec;
}

dcomp proj_pvec(int ml, vector<double>& pvec) {
	// Polarization Vector Projection to Spherical Harmonics
	try {
		if (ml == 1) return dcomp(-pvec[0]/sqrt(2),pvec[1]/sqrt(2));
		else if (ml == 0) return dcomp(pvec[2],0);
		else if (ml == -1) return dcomp(pvec[0]/sqrt(2),pvec[1]/sqrt(2));
		else throw invalid_argument("invalid ml for polarization vector projection");
		
	} catch (const exception &ex) {
		std::cout << ex.what() << "\n";
		exit(0);
	}
}

void XAS(double* SC, double* FG, double tenDQ, double lambda, vector<double>& pvec, int nedos) {
	// Scatters from p6d9 to p5d10
	double beta = 0;
	double racah_B = (SC[2]/49) - (5*SC[4]/441);
	norm_vec(pvec);

	Site cu_gs[1] = {Site("Ni")};
	Site cu_ex[1] = {Site("Ni",true,true)};

	int gs_size = cu_gs[0].hsize;
	int ex_size = cu_ex[0].hsize;

	double* gs_mat = new double[gs_size*gs_size]{0};
	double* gs_eigvec = new double[gs_size*gs_size]{0};
	double* ex_mat = new double[ex_size*ex_size]{0};
	double* ex_eigvec = new double[ex_size*ex_size]{0};
	double *gs_eig, *ex_eig;

	populate_hspace(cu_gs,1,gs_mat,SC,tenDQ,0);
	gs_eig = diagonalize(gs_mat,gs_eigvec,cu_gs[0].hsize,cu_gs[0].hsize); 
	populate_hspace(cu_ex,1,ex_mat,SC,tenDQ,lambda);
	cv_interacation(cu_ex[0],cu_ex[0],ex_mat,FG);
	ex_eig = diagonalize(ex_mat,ex_eigvec,cu_ex[0].hsize,cu_ex[0].hsize); 
	cout << cu_gs[0].hsize << " " << cu_ex[0].hsize << endl;

	dcomp* xas_mat = new dcomp[cu_gs[0].hsize*cu_ex[0].hsize];
	for (int i = 0; i < cu_gs[0].hsize*cu_ex[0].hsize; ++i) xas_mat[i] = 0;

	// Calculate XAS matrix element
	double gs_energy = gs_eig[0];
	vector<int> gs,ex; // index for ground state and excited states
	for (int init = 0; init < cu_gs[0].hsize; init++) {
		if (gs_eig[init] <= gs_energy) gs_energy = gs_eig[init];
	}
	cout << "GS energy: " << gs_energy << endl;
	// Calculate Partition function
	double partZ = 0;
	for (int init = 0; init < cu_gs[0].hsize; init++) {
		partZ += exp(-beta*gs_eig[init]);
		if (abs(gs_eig[init]-gs_energy) < 1e-7) {
			gs.push_back(init);
		}
	}

	cout << "partition function: " << partZ << endl;

	// Calculate Ground state J
	// double* gsev = new double[cu_gs[0].hsize];
	// for (int i = 0; i < cu_gs[0].hsize; ++i) gsev[i] = gs_eigvec[gs[0]*cu_gs[0].hsize+i];
	// vector<double> LS = cu_gs[0].get_hspace(3,2).momentum(cu_gs[0].simplify_state(gsev,3,2),false);
	// double Jgs = LS[0];
	// cout << "J: " << LS[0] << ", L: " << LS[1] << ", S: " << LS[2] << endl;

	// Create a map with index of Angular momentum and their index in Hilbert Space
	// for (int fin = 0; fin < cu_ex[0].hsize; fin++) {
	// 	//Calculate J
	// 	double Jex;
	// 	double* exev = new double[cu_ex[0].hsize];
	// 	for (int i = 0; i < cu_ex[0].hsize; ++i) exev[i] = ex_eigvec[fin*cu_ex[0].hsize+i];
	// 	for (auto & orb : cu_ex[0].orbs) {
	// 		LS = orb.momentum(cu_ex[0].simplify_state(exev,orb.get_n(),orb.l),false);
	// 		Jex += LS[0];
	// 		cout << "J: " << LS[0] << ", L: " << LS[1] << ", S: " << LS[2] << endl;
	// 	}
	// 	if (abs(abs(Jgs-Jex)-1) < 1e-5 || true) ex.push_back(fin);
	// }

	vector<double> E_ax,XAS;
	double emin = -15.0, emax = 15;
	for (double i = emin; i < emax; i += (emax-emin)/nedos) E_ax.push_back(i);
	XAS.resize(E_ax.size(),0.0);

	int num_transitions = 0;
	for (auto &g : gs) {
		// Evaluate Matrix Element
		for (int e = 0; e < cu_ex[0].hsize; e++) {
			if (abs(gs_eig[g]-ex_eig[e]) < emin || abs(gs_eig[g]-ex_eig[e]) > emax) continue;
			double cs = 0;
			for (int ee = 0; ee < cu_ex[0].hsize; ++ee) {
				if (abs(ex_eigvec[e*cu_ex[0].hsize+ee]) < 1e-7) continue;
				for (int gg = 0; gg < cu_gs[0].hsize; ++gg) {
					if (abs(gs_eigvec[g*cu_gs[0].hsize+gg]) < 1e-7) continue;
					int p2_mldiff = cu_gs[0].qn_get_state(2,1,gg) - cu_ex[0].qn_get_state(2,1,ee);
					int d3_mldiff = cu_ex[0].qn_get_state(3,2,ee) - cu_gs[0].qn_get_state(3,2,gg);
					struct QN qn2p[1] = {cu_ex[0].get_hspace(2,1).index2qn((int)log2(p2_mldiff))};
					struct QN qn3d[1] = {cu_gs[0].get_hspace(3,2).index2qn((int)log2(d3_mldiff))};
					if (qn2p[0].spin != qn3d[0].spin || !is_pw2(p2_mldiff) || !is_pw2(d3_mldiff) || gaunt(1,qn2p[0].ml,2,qn3d[0].ml)[1] == 0) continue;
					// cout << ", m2p: " << qn2p[0].ml << ",spin " << qn2p[0].spin << ", m3d: " << qn3d[0].ml << ",spin " << qn3d[0].spin << ", gaunt: " << gaunt(1,qn2p[0].ml,2,qn3d[0].ml)[1] << endl;
					dcomp ne = proj_pvec(qn3d[0].ml-qn2p[0].ml,pvec);
					QN* emptyqn;
					double fsign = cu_ex[0].get_hspace(2,1).Psign(qn2p,emptyqn,cu_ex[0].qn_get_state(2,1,ee),0,1,0) // 2p fermion sign
									* cu_gs[0].get_hspace(3,2).Psign(emptyqn,qn3d,0,cu_gs[0].qn_get_state(3,2,gg),0,1);
					// cout << "Fermion sign: " << fsign << endl;
					xas_mat[g+cu_gs[0].hsize*e] = gs_eigvec[g*cu_gs[0].hsize+gg]*ex_eigvec[e*cu_ex[0].hsize+ee]*gaunt(1,qn2p[0].ml,2,qn3d[0].ml)[1]*ne*fsign;
					// cout << "x: " << g << ", y: " << e << endl;
					// cout << "Ediff: " << ex_eig[e]-gs_eig[g] << ", index: " << (int)floor((ex_eig[e]-gs_eig[g]-emin)/((emax-emin)/nedos)) <<endl;
					cs += exp(-beta*gs_eig[g])*pow(abs(gs_eigvec[g*cu_gs[0].hsize+gg]*ex_eigvec[e*cu_ex[0].hsize+ee]*gaunt(1,qn2p[0].ml,2,qn3d[0].ml)[1]*ne),2); // Remember to add back in beta
				}
			}
			if (abs(cs) > 1e-4) {
				num_transitions += 1;
				cout << "cross section: " << cs << endl;
			}
			XAS[(int)floor((ex_eig[e]-gs_eig[g]-emin)/((emax-emin)/nedos))] += cs/partZ;
		}
	}

	cout << "Number of Transitions: " << num_transitions << endl;

	// Normalize XAS intensity

	for (auto &x : XAS) cout << x << " ";
	cout << endl;	
	write_mat(xas_mat,cu_gs[0].hsize,cu_ex[0].hsize,"xas.txt");



	write_mat(gs_eigvec,cu_gs[0].hsize,cu_gs[0].hsize,"gs_ev.txt");
	write_mat(ex_eigvec,cu_ex[0].hsize,cu_ex[0].hsize,"ex_ev.txt");
	cout << "ground state eigenenergies: ";
	printDistinct(gs_eig,cu_gs[0].hsize);
	cout << endl;
	cout << "excited states eigenenergies: ";
	printDistinct(ex_eig,cu_ex[0].hsize);
	cout << endl;


	return;
}
