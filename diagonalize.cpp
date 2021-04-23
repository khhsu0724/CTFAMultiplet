#include <iostream>
#include <fstream>
#include "diagonalize.h"


using namespace std;

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

double diagonalize(double *mat, int n, int m) {
	if (n != m) {
		cout << "Matrix is not square" << endl;
		return -1;
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
	if (info!=0){
		cout << "Error: dgeev returned error code " << info << endl;
		return -1;
	}
	// output eigenvalues to stdout
	cout << "--- Eigenvalues ---" << endl;
	for (int i=0;i<n;i++){
		cout << "( " << eigReal[i] << " , " << eigImag[i] << " )\n";
	}
	cout << endl;

	// deallocate
	delete [] mat;
	delete [] eigReal;
	delete [] eigImag;
	delete [] work;

	return 0;
}