#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <utility>
#include <iomanip>
#include <numeric>
#include <tuple>
#include "hilbert.h"

using namespace std;

int choose(int n, int k) {
	    if (k == 0) return 1;
	    return (n*choose(n-1, k-1))/k;
}

double norm(vector<double>& vin) {
	double n = 0;
	for (auto& v : vin) {
		n += v * v;
	}
	return sqrt(n);
}

double dot(vector<double>& a, vector<double>& b) {
	try {
		if (a.size() != b.size()) invalid_argument("different vector size for dot product");
		double dp = 0;
		for (int i = 0; i < a.size(); ++i) dp += a[i]*b[i];
		return dp;
	} catch (const exception &ex) {
		cout << ex.what() << "\n";
		exit(0);
	}
}

bool operator<(const QN& qn1, const QN& qn2) {
	// State further to the "right" is smaller
	if (qn1.spin != qn2.spin) {
		return qn1.spin < qn2.spin;
	} else {
		return qn1.ml > qn2.ml;
	}
}

bool operator!=(const QN& qn1, const QN& qn2) {
	if (qn1.spin == qn2.spin && qn1.ml == qn2.ml) return false;
	return true;
}

Hilbert::Hilbert(int n, char orb, int occ_num) : n(n), orb(orb), occ_num(occ_num) {
	try {
		if(orb == 's') {
			orb_avail = 2;
			l = 0;
			if (n < 1) throw invalid_argument("invalid quantum number");
		}
		else if(orb == 'p') {
			orb_avail = 6;
			l = 1;
			if (n < 2) throw invalid_argument("invalid quantum number");
		}
		else if(orb == 'd') {
			orb_avail = 10;
			l = 2;
			if (n < 3) throw invalid_argument("invalid quantum number");
		}
		else if(orb == 'f') {
			orb_avail = 14;
			l = 3;
			if (n < 4) throw invalid_argument("invalid quantum number");
		}
		else throw invalid_argument("invalid orbital");
		if (occ_num < 0 || occ_num > orb_avail) {
			throw invalid_argument("too many electrons for the orbital");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
		exit(0);
	}
	hmat_size = choose(orb_avail, occ_num);
	QN* ph;
	hmat = generate_states(0,ph);
	t_hsize = hmat_size;
	site_ind = vector<vector<int>>(hmat_size);
	for (int i = 0; i < hmat_size; ++i) site_ind[i].push_back(i);
	return;
}

double Hilbert::Psign(QN* lhsop, QN* rhsop, int lhss, int rhss, int lopnum, int ropnum) {
	// Count how many steps the operators need to traverse
	string bra, ket;
	if (lopnum != 0) bra = state2bit(lhss);
	if (ropnum != 0) ket = state2bit(rhss);
	int p = 0;
	for (int i = lopnum - 1; i >= 0; --i) {
		//both goes in reverse order
		int lindex = orb_avail - (l-lhsop[i].ml+orb_avail*(lhsop[i].spin+1)/4) - 1;
		for (int d = 0; d <= lindex; ++d) {
			if (bra[d] != '0') p += bra[d] - '0';
		}
		bra[lindex] = (bra[lindex]+1);
	}
	for (int i = ropnum - 1; i >= 0; --i) {
		//rhs
		int rindex = orb_avail - (l-rhsop[i].ml+orb_avail*(rhsop[i].spin+1)/4) - 1;
		for (int d = 0; d <= rindex; ++d) {
			if (ket[d] != '0') p += ket[d] - '0';
		}
		ket[rindex] = (ket[rindex]+1);
	}
	return pow(-1,p);
}

// Using 2^n representation to calculate quantumn number, convert to binary
// For example: 15 = 0000001111, meaning 3d orbital, states go from -2 to 2,
// and +1/2 spin to -1/2. Angular momentum comes first

int* Hilbert::generate_states(int in_state, QN* in_state_arr) {
	// Generate array of states in integer form, orb_avail choose occ_num
	// Can generate states that already includes input electron state
	try {
		if (in_state > occ_num) {
			throw invalid_argument("input states more than occupied states");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	int* mat = new int[choose(orb_avail - in_state,occ_num - in_state)]{0};
	string bitmask(occ_num - in_state, 1);
	string newmask;
    bitmask.resize(orb_avail - in_state, 0);
    int p = 0;
   	if (in_state != 0) sort(in_state_arr,in_state_arr+in_state);
    do {
    	mat[p] = 0;
    	newmask = bitmask;
    	for (int i = 0; i < in_state; ++i) {
    		struct QN input[1] = {in_state_arr[i]};
    		newmask.insert((size_t)qn2state(input,1,false),"1");
    	}
    	for (int i = 0; i < orb_avail; ++i) {
    		if (newmask[i]) mat[p] += pow(2,i);
    	}
    	++p;
    } while (prev_permutation(bitmask.begin(), bitmask.end()));
    return mat;
}

vector<pair<int,int>> Hilbert::match_states(int snum, QN* lhs, QN* rhs) {
	// Generate array pairs of states in integer form, orb_avail choose occ_num
	// matches state that contains either rhs or lhs, but have the same rest of the states
	int unique = 0;
	if (snum != 0) {
		sort(rhs,rhs+snum);
		sort(lhs,lhs+snum);
	}
	string u(orb_avail,'0');
	vector<int> mask;
	for (int i = 0; i < snum; ++i) {
		u[(l+rhs[i].ml)+orb_avail*(1-rhs[i].spin)/4] = '1';
		u[(l+lhs[i].ml)+orb_avail*(1-lhs[i].spin)/4] = '1';
	}
	for (int i = 0; i < orb_avail; ++i) {
		if (u[i] == '0') mask.push_back(orb_avail - i - 1);
		unique += u[i] - '0';
	}
	vector<pair<int,int>> mat;

	try {
		if (snum > occ_num) {
			throw invalid_argument("input states more than occupied states");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
		exit(0);
	}

	string bitmask(occ_num - snum, 1);
    bitmask.resize(orb_avail - unique, 0);
    do {
    	int r = qn2state(rhs,snum);
    	int l = qn2state(lhs,snum);
        for (int i = 0; i < orb_avail - unique; ++i) {
            if (bitmask[i]) {
            	r += pow(2,mask[i]);
            	l += pow(2,mask[i]);
            }
        }
        mat.push_back({l,r}); 
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
	return mat;
}

int Hilbert::bit2state(string& state) {
	// convert bit format to the state of atom to integer form
	try {
		if(orb_avail == state.length()) {
			return stoi(state,nullptr,2);
		}
		else {
			throw invalid_argument("invalid state converting bit to state");
		}
	}
	catch(const exception &ex) {
			std::cout << ex.what() << "\n";
	}
	return -1;
}

string Hilbert::state2bit(int state) {
	try {
		if(state <= pow(2,orb_avail) && state > 0) {
			string bit_state (orb_avail,'0');
			for (int i = orb_avail-1; i >= 0; --i) {
				if (floor(state/pow(2,i))) {
					bit_state[orb_avail-i-1] = '1';
					state -= pow(2,i);
				}
			}
			return bit_state;
		}
		else {
			throw invalid_argument("invalid state converting state to bit");
		}
	}
	catch(const exception &ex) {
			std::cout << ex.what() << "\n";
	}
	return "";
}

int Hilbert::qn2state(QN* qn, int snum, bool is_bit) {
	// down spin and higher ml first
	try {
		int state = 0;
		for (int i = 0; i < snum; ++i) {
			if (qn[i].spin != -1 && qn[i].spin != 1) {
				throw invalid_argument("invalid spin");
			}
			if (abs(qn[i].ml) > l) {
				throw invalid_argument("invalid angular momentum");
			}
			if (is_bit) state += pow(2,(l-qn[i].ml)+orb_avail*(qn[i].spin+1)/4);
			else state += (l-qn[i].ml)+orb_avail*(qn[i].spin+1)/4;
		}
		return state;
	}
	catch(const exception &ex) {
			std::cout << ex.what() << "\n";
	}
	return -1;
}

QN Hilbert::index2qn(int index) {
	// Remember that the index goes in reverse of qn ordering
	try {
		if (index >= orb_avail) throw invalid_argument("index larger than available orbitals");
		QN qn;
		qn.ml = (orb_avail-1-index)%(orb_avail/2)-l;
		if (index < orb_avail/2) qn.spin = 1;
		else qn.spin = -1;
		return qn;
	} catch(const exception &ex) {
			std::cout << ex.what() << "\n";
			exit(0);
	}
}

int Hilbert::sindex(int state) {
	try {
		if (state > pow(2,orb_avail) || state < 1) {
			throw invalid_argument("invalid state in hilbert space");
		}
		int index = distance(hmat, find(hmat, hmat + hmat_size, state));
		if (index < 0 || index >= hmat_size) {
			throw invalid_argument("state not in hilbert space");
		}
		return index;
	}
	catch(const exception &ex) {
			std::cout << ex.what() << "\n";
	}
	return -1;
}

void Hilbert::fill_in_mat(double* mat, double matelem, int lhs, int rhs) {
	// There is interaction when 2p differs but 3d matches
	for (int i = 0; i < site_ind[0].size(); ++i) {
		mat[site_ind[lhs][i]+t_hsize*site_ind[rhs][i]] += matelem;
	}
}

void Hilbert::pretty_print(double* mat, pair<int,int> column, pair<int,int> row) {
	int w = orb_avail+2;
	cout << endl <<  "Debug Matrix" << endl << string(w,' ');
	for (int c = column.first; c <= column.second; ++c) {
		cout << setw(w) << state2bit(hmat[c]);
	}
	cout << endl << endl;
	for (int r = row.first; r <= row.second; ++r) {
		cout << setw(w) << state2bit(hmat[r]);
		for (int c = column.first; c <= column.second; ++c) {
			cout << setw(w) << setprecision(4) << mat[c+hmat_size*r];
		}
		cout << endl << endl;
	}
}

vector<double> Hilbert::momentum(double* eigvec, bool return_square) {
	// Calculate L,S for specific eigenstate
	vector<double> Lz(hmat_size,0), Lp(hmat_size,0), Lm(hmat_size,0), Sz(hmat_size,0), Sp(hmat_size,0), Sm(hmat_size,0);
	double S2,L2,J2;
	if (occ_num == 0 || occ_num == orb_avail) return vector<double>{0.0,0.0,0.0};
	for (int j = 0; j < hmat_size; ++j) {
		if (abs(eigvec[j]) > 1e-7) {
			string state = state2bit(hmat[j]);
			for (int k = 0; k < orb_avail; ++k) {
				if (state[k] == '1') {
					double ml = (k % (orb_avail/2)) - l;
					Lz[j] += ml * eigvec[j];
					// Calculate L+
					string raise = state;
					if (ml < l && state[k+1] == '0') {
						raise[k] = '0';
						raise[k+1] = '1';
						Lp[sindex(bit2state(raise))] += sqrt((l-ml)*(l+ml+1)) * eigvec[j];
					}
					// Calculate L-
					string lower = state;
					if (ml > -l && state[k-1] == '0') {
						lower[k] = '0';
						lower[k-1] = '1';
						Lm[sindex(bit2state(lower))] += sqrt((l+ml)*(l-ml+1)) * eigvec[j];
					}
					// Calculate Spin
					if (k < (orb_avail/2)) {
						// Calculate Sz
						Sz[j] += 0.5 * eigvec[j];
						// Calculate S-
						int exchange = 0;
						string slower = state;
						if (state[k+orb_avail/2] == '0') {
							for (int e = k; e < k+orb_avail/2; ++e) {
								if (state[e] == '1') ++exchange;
							}
							slower[k] = '0';
							slower[k+orb_avail/2] = '1';
							Sm[sindex(bit2state(slower))] += pow(-1,exchange) * eigvec[j];
						}
					}
					else {
						// Calculate Sz
						Sz[j] += -0.5 * eigvec[j];
						// Calculate S+
						int exchange = 0;
						string sraise = state;
						if (state[k-orb_avail/2] == '0') {
							for (int e = k; e > k-orb_avail/2; --e) {
								if (state[e] == '1') ++exchange;
							}
							sraise[k] = '0';
							sraise[k-orb_avail/2] = '1';
							Sp[sindex(bit2state(sraise))] += pow(-1,exchange) * eigvec[j];
						}
					}
				}
			}
		} 
	}
	S2 = pow(norm(Sz),2) + 0.5 * (pow(norm(Sp),2) + pow(norm(Sm),2));
	L2 = pow(norm(Lz),2) + 0.5 * (pow(norm(Lp),2) + pow(norm(Lm),2));
	J2 = L2 + S2 + 2 * dot(Lz,Sz) + dot(Lm,Sp) + dot(Lp,Sm);
	if (abs(L2) < 1e-7) L2 = 0;
	if (abs(S2) < 1e-7) S2 = 0;
	if (return_square) return vector<double>{J2,L2,S2};
	else {
		double L = (-1+sqrt(1+4*L2))/2;
		double S = (-1+sqrt(1+4*S2))/2;
		double J = (-1+sqrt(1+4*J2))/2;
		return vector<double>{J,L,S};
	}
}

void Hilbert::momentum_check(double* mat, double* eig, double* eigvec) {
	// Check Lz, L2, Sz, S2 for the matrix
	double* L2 = new double[hmat_size]{0};
	double* S2 = new double[hmat_size]{0};
	cout << setw(10) << "eigval" << setw(10) << "J2" << setw(10) << "L2" << setw(10) << "S2" << endl;
	for (int i = 0; i < hmat_size; ++i) {
		double* j_evec = new double[hmat_size]{0};
		for(int j = 0; j < hmat_size; ++j) j_evec[j] = eigvec[i*hmat_size+j];
		vector<double> LS = momentum(j_evec,true);
		cout << setw(10) << eig[i] << setw(10) << LS[0] << setw(10) << LS[1] << setw(10) << LS[2] << endl;
	}
}

double Hilbert::pheshift(double trace, int k) {
	// return value for particle-hole energy shift for Coulomb Interaction
	if (occ_num > orb_avail/2) {
		int p = choose(occ_num,k);
		return trace/p*(p-choose(orb_avail-occ_num,k))/t_hsize;
	} else return 0;
}
