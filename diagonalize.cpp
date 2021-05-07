#include <iostream>
#include <fstream>
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
	double *input = new double[n*n];
	// Copy matrix for input
	for (int i = 0; i < n*m; ++i) {
		input[i] = mat[i];
	}

	double *eigReal = new double[n];
	double *eigImag = new double[n];
	double *vl = new double[n*n];
	int lwork=6*n;
	double *work=new double[lwork];
	int info;

	// calculate eigenvalues using the DGEEV subroutine
	dgeev_(&Nchar,&Nchar,&n,input,&n,eigReal,eigImag,
			vl,&n,vr,&n,work,&lwork,&info);
	// check for errors
	try {
		if (info!=0){
			throw runtime_error( "Error: dgeev returned error code ");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	// output eigenvalues to stdout
	// cout << "--- Eigenvalues ---" << endl;
	// for (int i=0;i<n;i++){
	// 	cout << "( " << eigReal[i] << " , " << eigImag[i] << " )\n";
	// }
	// cout << endl;

	// deallocate
	// delete [] eigReal;
	// delete [] eigImag;
	delete [] work;

	return eigReal;
}