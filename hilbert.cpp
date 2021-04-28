#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <iomanip>
#include "hilbert.h"

using namespace std;

int choose(int n, int k) {
	    if (k == 0) return 1;
	    return (n*choose(n-1, k-1))/k;
}

bool operator<(const QN& qn1, const QN& qn2) {
	// State further to the "right" is smaller
	if (qn1.spin != qn2.spin) {
		return qn1.spin < qn2.spin;
	} else {
		return qn1.ml > qn2.ml;
	}
}

Hilbert::Hilbert(char orb, int occ_num) : orb(orb), occ_num(occ_num) {
	try {
		if(orb == 's') {
			orb_avail = 2;
			n = 0;
		}
		else if(orb == 'p') {
			orb_avail = 6;
			n = 1;
		}
		else if(orb == 'd') {
			orb_avail = 10;
			n = 2;
		}
		else if(orb == 'f') {
			orb_avail = 14;
			n = 3;
		}
		else {
			throw invalid_argument("invalid orbital");
		}
		if (occ_num < 2 || occ_num > orb_avail-2) {
			throw invalid_argument("occupied electron number will not form a multiplet");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	l = n;
	hmat_size = choose(orb_avail, occ_num);
	QN* ph;
	hmat = generate_states(0,ph);
}

// Using 2^n representation to calculate quantumn number, convert to binary
// For example: 15 = 0000001111, meaning 3d orbital, all states in ml = 1,2
// are occupied.

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
	int* mat = new int[choose(orb_avail - in_state,occ_num - in_state)];
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

int Hilbert::bit2state(const char* state) {
	// convert bit format to the state of atom to integer form
	try {
		if(orb_avail == strlen(state)) {
			return stoi(state,nullptr,2);
		}
		else {
			throw invalid_argument("invalid state");
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
			string bit_state (10,'0');
			for (int i = orb_avail-1; i >= 0; --i) {
				if (floor(state/pow(2,i))) {
					bit_state[orb_avail-i-1] = '1';
					state -= pow(2,i);
				}
			}
			return bit_state;
		}
		else {
			throw invalid_argument("invalid state");
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
