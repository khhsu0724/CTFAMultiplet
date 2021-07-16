#include <iostream> 
#include <fstream>
#include <cmath>
#include "diagonalize.h"


using namespace std;

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
}


double trace(double *mat, int n, int m) {
	try {
		if (n != m){
			throw invalid_argument( "Matrix is not square");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	double t = 0;
	for (int i = 0; i < n; ++i) {
		t += mat[i+n*i];
	}
	return t;
}

void norm_ev(double* ev, int n) {
	for (int i = 0; i < n; ++i) {
		double nm = 0;
		for (int j = 0; j < n; j++) {
			if (abs(ev[i*n+j]) < 1e-7) ev[i*n+j] = 0;
			nm += ev[i*n+j]*ev[i*n+j];
		}
		if (nm != 1) {
			nm = sqrt(nm);
			for (int j = 0; j < n; j++) ev[i*n+j] /= nm;
		}
	}
	return;
}

double* diagonalize(double *mat, double *vr, int n, int m) {
	try {
		if (n != m){
			throw invalid_argument( "Matrix is not square");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	// allocate data
	char Nchar='V';
	double *input = new double[n*n]{0};
	// Copy matrix for input
	for (int i = 0; i < n*m; ++i) {
		input[i] = mat[i];
	}

	double *eigReal = new double[n]{0};
	double *eigImag = new double[n]{0};
	double *vl = new double[n*n]{0};
	int lwork=6*n;
	double *work=new double[lwork]{0};
	int info;

	// calculate eigenvalues using the DGEEV subroutine
	dgeev_(&Nchar,&Nchar,&n,input,&n,eigReal,eigImag,
			vl,&n,vr,&n,work,&lwork,&info);
	norm_ev(vr,n);
	// check for errors
	try {
		if (info!=0){
			throw runtime_error( "Error: dgeev returned error code ");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}

	delete [] work;

	return eigReal;
}