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

double Psign(QN qn1, QN qn2, QN qn3, QN qn4) {
	// Make sure that the fermion operator is normal ordered, uses different rule to compare QN
	double p = 1;
	if (qn4 < qn3) p *= -1;
	if (qn1 < qn2) p *= -1;
	return p;
}

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
	double n = hspace.hmat_size, m = n;
	double* ml_arr = new double[hspace.l*2+1];
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
			int num_states = choose(hspace.orb_avail-2,hspace.occ_num-2);
			struct QN qn12[2] = {{m12.first,-1},{m12.second,1}}; // psi,i,j => Creation operators m4,m3
			int* state12 = hspace.generate_states(2,qn12);
			for (auto m34 : mpair_a) {
				struct QN qn34[2] = {{m34.first,-1},{m34.second,1}}; // Annihilation operators m1,m2
				int* state34 = hspace.generate_states(2,qn34);
				for (int s1 = 0; s1 < num_states; ++s1) {
					for (int s2 = 0; s2 < num_states; ++s2) {
						double matelem = calc_U(gaunt(hspace.l,m34.first,hspace.l,m12.first),
									gaunt(hspace.l,m12.second,hspace.l,m34.second),SC,2*hspace.l+1);
						if (hspace.qn2state(qn12,2) != hspace.qn2state(qn34,2)) { // Off Diagonal Term
							matelem *= Psign(qn34[0],qn34[1],qn12[1],qn12[0]);
						}
						mat[hspace.sindex(state12[s1])+hspace.hmat_size*hspace.sindex(state34[s2])] += matelem;
						// Debug

							if (hspace.sindex(state12[s1]) == 1 && hspace.sindex(state34[s2]) == 0) {
								cout << "m1(sd): " << m34.first << ", m2(su): " << m34.second << ", m3(su): " << m12.second << ", m4(sd): " << m12.first;
								cout << ", gaunt sum: " << calc_U(gaunt(hspace.l,m34.first,hspace.l,m12.first),
									gaunt(hspace.l,m12.second,hspace.l,m34.second),SC,2*hspace.l+1) 
								<< ", index: " << hspace.sindex(state12[s1])+hspace.hmat_size*hspace.sindex(state34[s2]) << endl;
							}
					}
				}
			}
		}

		for (auto spin : {-1,1}) {
			for (auto m12 : mpair_p) {
				// Populate Hamiltonian that contains the parallel spin pair
				int num_states = choose(hspace.orb_avail-2,hspace.occ_num-2);
				struct QN qn12[2] = {{m12.first,spin},{m12.second,spin}};
				int* state12 = hspace.generate_states(2,qn12);
				for (auto m34 : mpair_p) {
					struct QN qn34[2] = {{m34.first,spin},{m34.second,spin}};
					int* state34 = hspace.generate_states(2,qn34);
					for (int s1 = 0; s1 < num_states; ++s1) {
						for (int s2 = 0; s2 < num_states; ++s2) {
							double matelem = calc_U(gaunt(hspace.l,m34.first,hspace.l,m12.first),
										gaunt(hspace.l,m12.second,hspace.l,m34.second),SC,hspace.l+hspace.l+1) -
								calc_U(gaunt(hspace.l,m34.second,hspace.l,m12.first),
										gaunt(hspace.l,m12.second,hspace.l,m34.first),SC,hspace.l+hspace.l+1);
							if (hspace.qn2state(qn12,2) != hspace.qn2state(qn34,2)) { // Off Diagonal Term
								matelem *= Psign(qn34[0],qn34[1],qn12[1],qn12[0]);
							}
							mat[hspace.sindex(state12[s1])+hspace.hmat_size*hspace.sindex(state34[s2])] += matelem;

							//Debug				
							if (hspace.sindex(state12[s1]) == 1 && hspace.sindex(state34[s2]) == 0) {
								cout << "spin: " << spin;
								cout << ", m1: " << m34.first << ", m2: " << m34.second << ", m3: " << m12.second << ", m4: " << m12.first;
								cout << ", gaunt sum: " << calc_U(gaunt(hspace.l,m34.first,hspace.l,m12.first),
										gaunt(hspace.l,m12.second,hspace.l,m34.second),SC,hspace.l+hspace.l+1) -
								calc_U(gaunt(hspace.l,m34.second,hspace.l,m12.first),
										gaunt(hspace.l,m12.second,hspace.l,m34.first),SC,hspace.l+hspace.l+1) 
								<< ", index: " << hspace.sindex(state12[s1])+hspace.hmat_size*hspace.sindex(state34[s2]) << endl;
							}
						}
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
		double* eig;
		eig = diagonalize(mat,d2.hmat_size,d2.hmat_size);
		myfile << cf << " ";
		for (int i =0; i < d2.hmat_size; ++i) {
			myfile << eig[i] << " ";
		}
		myfile << endl;
	}
    myfile.close();
}

void printDistinct(double arr[], int n)
{
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

int main(int argc, char** argv){
	// Process flags

	Hilbert d2('d',3);
	double* mat = new double[d2.hmat_size*d2.hmat_size];
	for (int i = 0; i < d2.hmat_size*d2.hmat_size; ++i) {
		mat[i] = 0;
	}
	// struct QN qn_list[2] = {
	// 	{2,-1},{-1,1},
	// };
	// int* test = d2.generate_states(2,qn_list);
	// cout << choose(d2.orb_avail-2,d2.occ_num-2) << endl;
	// for (int i = 0; i < choose(d2.orb_avail-2,d2.occ_num-2); ++i) {
	// 	cout << d2.state2bit(test[i]) << endl;
	// }
	// for (int i = 0; i < d2.hmat_size; ++i) {
	// 	cout << d2.hmat[i] << endl;
	// }
	double SC[5] = {1.0,0,1.0*49,0,1*441};
	calc_coulomb(d2, mat, SC);
	d2.pretty_print(mat,{0,9},{0,9});
	double* eig;
	eig = diagonalize(mat,d2.hmat_size,d2.hmat_size);
	// calc_CF(d2,mat,0);
	printDistinct(eig,d2.hmat_size);
	cout << endl << "trace: " << trace(mat,45,45) << endl;
	// d2.pretty_print(mat,{0,8},{0,8});
	delete [] mat;
}