#include "helper.hpp"

int ed::choose(int n, int k) {
	if (k > n) return 0;
	if (k == 0) return 1;
    return (n*choose(n-1, k-1))/k;
}

bool ed::is_pw2(int x) {return !(x == 0) && !(x & (x - 1));}

ulli ed::next_perm(ulli v) {
	ulli t = v | (v - 1);
	return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(v) + 1));
}

// Enumarate states that include bits for "inc" variable using recursion
void ed::enum_states(std::vector<ulli>& states, ulli n, ulli k, ulli inc, ulli s) {
	if (k > n) return;
	if (k == 0 && n < (__builtin_ctzll(inc)+1)) {
		states.emplace_back(s);
		return;
	}
	enum_states(states,n-1,k-1,inc,s+(1<<(n-1))); // bit is 1
	if (!(inc & (1<<(n-1)))) enum_states(states,n-1,k,inc,s); // bit is 0
}

// Add valence state and core state 
ulli ed::add_bits(ulli b1, ulli b2, int b1size, int b2size) {
	ulli b1d = b1 & ((1U << (b1size/2)) - 1);
	ulli b2d = b2 & ((1U << (b2size/2)) - 1);
	return (b1^b1d) << b2size | (b2^b2d) << (b1size/2) | b1d << (b2size/2) | b2d;
}

int ed::count_bits(ulli b) {
	int c;
	for (c = 0; b; c++) b &= b - 1;
	return c;
}

void ed::sph2real(double* sph, double* tet) {
/** Converts spherical harmonics to real space
	https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
	Input: number of arrays, Array of orbital (su,sd), and orbital number
	Returns: 2*(2l+1) size array for occupancy **/

	// Cast double matrix from real to complex
	int l = 0;
	if (l == 1) {
		// returns x,y,z
		dcomp trans_mat[] = {{1,0},{0,0},{1,0},
					   		 {0,1},{0,0},{0,-1},
					   		 {0,0},{sqrt(2),0},{0,0}};

	} else if (l == 2) {
		// returns xy,yz,z^2,xz,x^2-y^2
		dcomp trans_mat[] = {{0,1},{0,0},{0,0},{0,0},{0,1},
					         {0,0},{0,1},{0,0},{0,1},{0,0},
					         {0,0},{0,0},{sqrt(2),0},{0,0},{0,0},
					         {0,0},{1,0},{0,0},{-1,0},{0,0},
					         {1,0},{0,0},{0,0},{0,0},{1,0}};

	}

	return;
}
