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

void write_mat(double* mat, size_t n, string file_dir) {
	ofstream matfile;
    matfile.open (file_dir);
    for (int i = 0; i < n; ++i) {
    	for (int j = 0; j < n; ++j) {
    		matfile << mat[i+n*j] << " ";
    	}
    	matfile << endl;
    }
    matfile.close();
}

vector<double> printDistinct(double arr[], int n, bool is_print = true) {
    // Pick all elements one by one
    vector<double> unique_eig;
    for (int i=0; i<n; i++)
    {
        // Check if the picked element is already printed
        int j;
        for (j=0; j<i; j++)
           if (abs(arr[i] - arr[j]) < 1e-10)
               break;
        // If not printed earlier, then print it
        if (i == j) {
        	if (is_print) cout << arr[i] << " ";
        	unique_eig.push_back(arr[i]);
        }
    }
    return unique_eig;
}

double calc_U(double* gaunt1, double* gaunt2, double* SC, int size) {
	double u = 0;
	for (int i = 0; i < size; ++i) {
		u += gaunt1[i] * gaunt2[i] * SC[i];
	}
	return u;
}

void calc_coulomb(Hilbert& hspace, double* mat, double* SC) {
	//Calculate Coulomb Matrix Element
	if (hspace.occ_num < 2 || hspace.occ_num > hspace.l - 2) return;
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
				cout << "anti, m12: " << m12.first << "," << m12.second << ", m34: " << m34.first << "," << m34.second << endl;
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
					cout << "par, m12: " << m12.first << "," << m12.second << ", m34: " << m34.first << "," << m34.second << endl;
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

double* CFmat(int l, double tenDQ, int coord) {
	// Crystal field present when the 2 state has same spin
	try {
		if (l == 2 && coord == 6) {
			double* cfmat = new double[25];
			for (int i = 0; i < 25; ++i) cfmat[i] = 0;
			cfmat[0] = 1 * tenDQ;
			cfmat[1+5*1] = -4 * tenDQ;
			cfmat[2+5*2] = 6 * tenDQ;
			cfmat[3+5*3] = -4 * tenDQ;
			cfmat[4+5*4] = 1 * tenDQ;
			cfmat[4+5*0] = 5 * tenDQ;
			cfmat[0+5*4] = 5 * tenDQ;
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
					mat[hspace.sindex(e.first)+hspace.hmat_size*hspace.sindex(e.second)] += cfmat[i];
				}
			}
		}
	}
}

void calc_SO(Hilbert& hspace, double* mat, double lambda) {
	// Calculate Spin Orbit Coupling Matrix Element
	for (int ml = -hspace.l; ml <= hspace.l; ++ml) {
		// parallel
		for (auto spin: {-1,1}) {
			struct QN qn[1] = {{ml,spin}};
			vector<pair<int,int>> entries = hspace.match_states(1, qn, qn);
			for (auto& e : entries) {
				double matelem = -lambda/2 * spin * ml; // No sign traversing issue due to symmetry
				mat[hspace.sindex(e.first)+hspace.hmat_size*hspace.sindex(e.second)] += matelem;
			}
		}
		// Antiparallel
		if (ml < hspace.l) {
			struct QN qn1[1] = {{ml+1,-1}};
			struct QN qn2[1] = {{ml,1}};
			vector<pair<int,int>> entries = hspace.match_states(1, qn1, qn2);
			for (auto& e : entries) {
				double matelem = -lambda/2 * sqrt((hspace.l-ml)*(hspace.l+ml+1));
				mat[hspace.sindex(e.first)+hspace.hmat_size*hspace.sindex(e.second)] 
					+= matelem * hspace.Psign(qn1,qn2,e.first,e.second,1);
				mat[hspace.sindex(e.second)+hspace.hmat_size*hspace.sindex(e.first)] 
					+= matelem * hspace.Psign(qn2,qn1,e.second,e.first,1);
			}
		}
	}
}

void tanabe_sugano(double* SC) {
	// Terrible code now, should implement plotting & block diagonal matrices in the future
	ofstream myfile;
    myfile.open ("./tb.txt");

	Hilbert d2('d',5);
	double racah_B = (SC[2]/49) - (5*SC[4]/441);
	cout << "racah_B: " << racah_B << endl;
	for (double cf = 0; cf <= 3.0; cf += 0.06) {
		double* mat = new double[d2.hmat_size*d2.hmat_size];
		double* eigvec = new double[d2.hmat_size*d2.hmat_size];
		for (int i = 0; i < d2.hmat_size*d2.hmat_size; ++i) {
			mat[i] = 0;
			eigvec[i] = 0;
		}
		calc_coulomb(d2, mat, SC);
		calc_CF(d2,mat,cf*racah_B,6);
		double* eig;
		eig = diagonalize(mat,eigvec,d2.hmat_size,d2.hmat_size);
		vector<double> unique_eig = printDistinct(eig,d2.hmat_size,false);
		cout << unique_eig.size() << ", ";
		if (unique_eig.size() != 10) {
			write_mat(mat,d2.hmat_size,"./mat.txt");
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

	Hilbert d2('p', 5);
	double* mat = new double[d2.hmat_size*d2.hmat_size];
	double* eigvec = new double[d2.hmat_size*d2.hmat_size];
	for (int i = 0; i < d2.hmat_size*d2.hmat_size; ++i) {
		mat[i] = 0;
		eigvec[i] = 0;
	}
	double SC[5] = {1.0,0,1.0*49,0,1*441};
	double racah_B = (SC[2]/49) - (5*SC[4]/441);
	calc_coulomb(d2, mat, SC);
	calc_SO(d2, mat, 2); 
	// calc_CF(d2,mat,0.1*racah_B,6);

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