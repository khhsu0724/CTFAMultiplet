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


using namespace std;

double calc_U(double* gaunt1, double* gaunt2, double* SC, int size) {
	double u = 0;
	for (int i = 0; i < size; ++i) {
		u += gaunt1[i] * gaunt2[i] * SC[i];
	}
	return u;
}


// Implement Symmetries
void calc_coulomb(Hilbert& hspace, double* mat, double* SC) {
	//Calculate Coulomb Matrix Element
	int n = hspace.hmat_size, m = n;
	int* ml_arr = new int[hspace.l*2+1];
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
						matelem *= hspace.Psign(qn12,qn34,e.first,e.second,2);
					}
					mat[hspace.sindex(e.first)+hspace.hmat_size*hspace.sindex(e.second)] += matelem;
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
							matelem *= hspace.Psign(qn12,qn34,e.first,e.second,2);
						}
						mat[hspace.sindex(e.first)+hspace.hmat_size*hspace.sindex(e.second)] += matelem;
					}
				}
			}
		}
	}
}

double* CFmat(double l, double tenDQ) {
	// Crystal field present when the 2 state has same spin
	try {
		if (l == 2) {
			double* cfmat = new double[25];
			for (int i = 0; i < 5; ++i) cfmat[i] = 0;
			cfmat[0] = 1 * tenDQ * 0.1;
			cfmat[1+5*1] = -4 * tenDQ * 0.1;
			cfmat[2+5*2] = 6 * tenDQ * 0.1;
			cfmat[3+5*3] = -4 * tenDQ * 0.1;
			cfmat[4+5*4] = 1 * tenDQ * 0.1;
			cfmat[4+5*0] = 5 * tenDQ * 0.1;
			cfmat[0+5*4] = 5 * tenDQ * 0.1;
			return cfmat;
		} else {
			throw invalid_argument("non-d crystal field not coded yet");
		}
	} catch(const exception &ex) {
			std::cout << ex.what() << "\n";
	}
	double* ph;
	return ph;
}

void calc_CF(Hilbert& hspace, double* mat, double tenDQ) {
	//Calculate Crystal Field Matrix Element
	double* cfmat = CFmat(hspace.l, tenDQ);
	for (int i = 0; i < (hspace.l*2+1)*(hspace.l*2+1); ++i) {
		if (cfmat[i] != 0) {
			int num_states = choose(hspace.orb_avail-1,hspace.occ_num-1);
			int l1 = i%(hspace.l*2+1)-2;
			int l2 = i/(hspace.l*2+1)-2;
			if (l1 >= l2) {
				for (auto spin : {1,-1}) {
					struct QN qn1[1] = {{l1,spin}};
					struct QN qn2[1] = {{l2,spin}};
					int* state1 = hspace.generate_states(1,qn1);
					int* state2 = hspace.generate_states(1,qn2);
 					for (int s1 = 0; s1 < num_states; ++s1) {
						for (int s2 = 0; s2 < num_states; ++s2) {
							mat[hspace.sindex(state1[s1])+hspace.hmat_size*hspace.sindex(state2[s2])] 
								+= cfmat[i];
						}
					}
				}
			}
		}
	}
}

void calc_SO() {
	//Calculate Spin Orbit Coupling Matrix Element
}

void tanabe_sugano(double* SC) {
	ofstream myfile;
    myfile.open ("./tb.txt");

	Hilbert d2('d',2);
	double* mat = new double[d2.hmat_size*d2.hmat_size];
	for (double cf = 0; cf < 3.0; cf+=0.1) {
		for (int i = 0; i < d2.hmat_size*d2.hmat_size; ++i) {
			mat[i] = 0;
		}
		calc_coulomb(d2, mat, SC);
		calc_CF(d2,mat,cf);
		// double* eig;
		// eig = diagonalize(mat,d2.hmat_size,d2.hmat_size);
		// myfile << cf << " ";
		// for (int i =0; i < d2.hmat_size; ++i) {
		// 	myfile << eig[i] << " ";
		// }
		// myfile << endl;
	}
    myfile.close();
}

void printDistinct(double arr[], int n) {
    // Pick all elements one by one
    for (int i=0; i<n; i++)
    {
        // Check if the picked element is already printed
        int j;
        for (j=0; j<i; j++)
           if (abs(arr[i] - arr[j]) < 0.00001)
               break;
        // If not printed earlier, then print it
        if (i == j)
          cout << arr[i] << " ";
    }
}

void block_diag_hack(double* mat, Hilbert& hspace, int L, double S) {
	vector<int> bh;
	int c = 0;
	for (int i = 0; i < hspace.hmat_size; ++i) {
		string s = hspace.state2bit(hspace.hmat[i]);
	 	int ml = 0;
	 	double spin = 0;
		 for (int j = 0; j < hspace.orb_avail; ++j) {
			if (s[j] == '1') {
				ml += j % 5 - 2;
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
		cout << eig[i] << " ";
	}
	cout << endl << "distinct eigenvalues: ";
	printDistinct(eig,bh.size());
}

int main(int argc, char** argv){
	// Process flags

	bool CF_on = false;
	bool SO_on = false;

	Hilbert d2('d',7);
	double* mat = new double[d2.hmat_size*d2.hmat_size];
	double* eigvec = new double[d2.hmat_size*d2.hmat_size];
	for (int i = 0; i < d2.hmat_size*d2.hmat_size; ++i) {
		mat[i] = 0;
		eigvec[i] = 0;
	}
	double SC[5] = {1.0,0,1.0*49,0,1.0*441};
	calc_coulomb(d2, mat, SC);

	// d2.pretty_print(mat,{0,9},{0,9});
	double* eig;
	eig = diagonalize(mat,eigvec,d2.hmat_size,d2.hmat_size);
	double esum;
	for (int i = 0; i < d2.hmat_size; ++i) {
		esum += eig[i];
	}

	printDistinct(eig,d2.hmat_size);
	cout << endl;
	// d2.momentum_check(mat,eig,eigvec);
	// block_diag_hack(mat,d2,3,0.5);
	delete [] mat;
}