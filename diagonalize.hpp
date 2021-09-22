#ifndef DIAGONALIZE
#define DIAGONALIZE
#include "helper.hpp"

double* dgeev(double *mat, double *vr, int n);
double* dsyev(double *mat, double *vr, int n);

class Block {
private:
	double Sz = 0, Lz = 0, K = 0;
public:
	Block(double Sz, double Lz, double K, int size): Sz(Sz), Lz(Lz), K(K), size(size) {
		// We can preallocate this with a number. i.e. Hilbert space size/number of blocks * 1.2
		ham = new double[size*size]{0};
	};
	double get_Sz() {return Sz;};
	double get_Lz() {return Lz;};
	double get_K() {return K;};	
	int size;
	double *ham, *eig, *eigvec;
	void diag_dgeev() {
		eigvec = new double[size*size]{0};
		eig = new double[size]{0};
		eig = dgeev(ham,eigvec,size);
		delete[] ham;
	}
	void diag_dsyev() {
		eig = dsyev(ham,eigvec,size);
		eigvec = ham;
	}		
};

#endif 