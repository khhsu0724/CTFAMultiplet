#include <iostream>
#include <fstream>
#include "diagonalize.h"


using namespace std;

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

void diagonalize_test() {
	// Generate random array
	int n,m;
	double *mat;

	n = 10;
	m = 10;

	mat = new double[n*m];

	srand(time(NULL));
	for (int i = 0; i < n*m; ++i) {
	  mat[i] = rand() % 10 + 1;
	}
	cout << "print matrix" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			cout << mat[i*n+j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	diagonalize(mat,n,m);
	delete [] mat;
}

double trace(double *mat, int n, int m) {
	try {
		if (n != m){
			throw invalid_argument( "Matrix is not square");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	int t = 0;
	for (int i = 0; i < n; ++i) {
		t += mat[i+n*i];
	}
	return t;
}

double* diagonalize(double *mat, int n, int m) {
	try {
		if (n != m){
			throw invalid_argument( "Matrix is not square");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	// allocate data
	char Nchar='N';
	double *eigReal=new double[n];
	double *eigImag=new double[n];
	double *vl,*vr;
	int one=1;
	int lwork=6*n;
	double *work=new double[lwork];
	int info;

	// calculate eigenvalues using the DGEEV subroutine
	dgeev_(&Nchar,&Nchar,&n,mat,&n,eigReal,eigImag,
			vl,&one,vr,&one,
			work,&lwork,&info);
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