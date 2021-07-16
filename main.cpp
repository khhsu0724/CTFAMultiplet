#include <iostream>
#include <fstream>
#include <cstdlib> 
#include <ctime> 
#include <string>
#include <vector>
#include <utility>
#include <stdlib.h>
#include <algorithm>
#include "diagonalize.h"
#include "gaunt.h"
#include "hilbert.h"
#include "site.h"
#include "multiplet.h"
#include "photon.h"

using namespace std;

// void write_mat(double* mat, size_t n, string file_dir) {
// 	ofstream matfile;
//     matfile.open (file_dir);
//     for (int i = 0; i < n; ++i) {
//     	for (int j = 0; j < n; ++j) {
//     		matfile << mat[i+n*j] << " ";
//     	}
//     	matfile << endl;
//     }
//     matfile.close();
// }

// vector<double> printDistinct(double arr[], int n, bool is_print = true) {
//     // Pick all elements one by one
//     vector<double> unique_eig;
//     for (int i=0; i<n; i++)
//     {
//         // Check if the picked element is already printed
//         int j;
//         for (j=0; j<i; j++)
//            if (abs(arr[i] - arr[j]) < 1e-10)
//                break;
//         // If not printed earlier, then print it
//         if (i == j) {
//         	if (is_print) cout << arr[i] << " ";
//         	unique_eig.push_back(arr[i]);
//         }
//     }
//     return unique_eig;
// }

void tanabe_sugano(double* SC) {
	// Terrible code now, should implement plotting & block diagonal matrices in the future
	ofstream myfile;
    myfile.open ("./tb.txt");

	Hilbert d2(3,'d',5);
	double racah_B = (SC[2]/49) - (5*SC[4]/441);
	cout << "racah_B: " << racah_B << endl;
	for (double cf = 0; cf <= 3.0; cf += 0.06) {
		double* mat = new double[d2.hmat_size*d2.hmat_size]{0};
		double* eigvec = new double[d2.hmat_size*d2.hmat_size]{0};
		calc_coulomb(d2, mat, SC);
		calc_CF(d2,mat,cf*racah_B,6);
		double* eig;
		eig = diagonalize(mat,eigvec,d2.hmat_size,d2.hmat_size);
		vector<double> unique_eig = printDistinct(eig,d2.hmat_size,false);
		cout << unique_eig.size() << ", ";
		if (unique_eig.size() != 10) {
			write_mat(mat,d2.hmat_size,d2.hmat_size,"./mat.txt");
		}
		auto min_eig = *min_element(unique_eig.begin(), unique_eig.end());
		sort(unique_eig.begin(), unique_eig.end());
		for (auto& e : unique_eig) {
			cout << e - min_eig << " ";
			myfile << e - min_eig << " ";
		}
		cout << endl;
		myfile << endl;
		delete [] mat;
		delete [] eigvec;
	}
    myfile.close();
}

void block_diag_hack(double* mat, Hilbert& hspace, int L, double S) {
	// Check block matrices, need to phase out in the future
	vector<int> bh;
	int c = 0;
	for (int i = 0; i < hspace.hmat_size; ++i) {
		string s = hspace.state2bit(hspace.hmat[i]);
	 	int ml = 0;
	 	double spin = 0;
		 for (int j = 0; j < hspace.orb_avail; ++j) {
			if (s[j] == '1') {
				ml += j % (hspace.orb_avail/2)- hspace.l;
				if (j < hspace.orb_avail/2) spin += -0.5;
				else spin += 0.5;
			}
		}
		if (ml == L && spin == S) {
			++c;
			bh.push_back(hspace.hmat[i]);
		}
	}
	double* bmat = new double[bh.size()*bh.size()];
	for (int i = 0; i < bh.size(); ++i) {
		for (int j = 0; j < bh.size(); ++j) {
			bmat[i+bh.size()*j] = mat[hspace.sindex(bh[i])+hspace.hmat_size*hspace.sindex(bh[j])];
		} 
	}

	double* eigvec = new double[bh.size()*bh.size()];
	double* eig;
	eig = diagonalize(bmat,eigvec,bh.size(),bh.size());
	for (int i = 0; i < bh.size(); ++i) {
		for (int j = 0; j < bh.size(); ++j) {
			cout << eigvec[i*bh.size()+j] << " ";
		}
		cout << endl;
	}
	for (int i = 0; i < bh.size(); ++i) {
		cout << eig[i] << " ";
	}
	cout << endl << "distinct eigenvalues: ";
	printDistinct(eig,bh.size());
	cout << endl;
}

int main(int argc, char** argv){
	// Process flags

	bool CF_on = false;
	bool SO_on = false;

	Hilbert d2(3,'d',2);
	double* mat = new double[d2.hmat_size*d2.hmat_size]{0};
	double* eigvec = new double[d2.hmat_size*d2.hmat_size]{0};
	// double SC[5] = {1.0,0,1.0*49,0,1.0*441};
	double SC[5] = {6.0,0,0.13*49,0,0.025*441}; //F0,2,4 = 6, 0.13, 0.025 eV
	double FG[4] = {6.18,6.2,6.18,4.63}; // Core-Valence Interaction Parameter
	// double SC[5] = {1.0,0,1.0*49,0,1.0*441};
	double racah_B = (SC[2]/49) - (5*SC[4]/441);
	calc_coulomb(d2, mat, SC);
	// calc_SO(d2,mat,2); 
	// calc_CF(d2,mat,1,6);

	double* eig;
	eig = diagonalize(mat,eigvec,d2.hmat_size,d2.hmat_size);
	double esum;
	for (int i = 0; i < d2.hmat_size; ++i) {
		esum += eig[i];
	}

	cout << "eigenvalues: ";
	printDistinct(eig,d2.hmat_size);
	cout << endl;
	// block_diag_hack(mat,d2,5,0.5);
	// d2.momentum_check(mat,eig,eigvec);

	// write_mat(mat,d2.hmat_size,d2.hmat_size,"p2mat.txt");
	write_mat(eigvec,d2.hmat_size,d2.hmat_size,"d4ev.txt");
	vector<double> pvec = {1,1,1};
	XAS(SC,FG,0,10.5,pvec,100);
	// block_diag_hack(mat,d2)
	// block_diag_hack(mat,d2,3,0.5);
	delete [] mat;
}