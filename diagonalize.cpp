#include <iostream> 
#include <cmath>
#include "diagonalize.hpp"
#include "helper.hpp"

using namespace std;

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
	extern void dgeev_(char*,char*,int*,double*,int*,double*,double*, 
						double*,int*,double*,int*,double*,int*,int*);
}
extern "C" {
	extern void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
}
// arpack implementation

double* dgeev(double *mat, double *vr, int n) {
	// Slow but robust, don't use this
	char JOBVL='N',JOBVR='V';
	double *input = new double[n*n]{0};
	for (int i = 0; i < n*n; ++i) input[i] = mat[i];

	double *eigReal = new double[n]{0};
	double *eigImag = new double[n]{0};
	double *vl = new double[n*n]{0};
	int lwork = 6*n;
	double *work = new double[lwork]{0};
	int info;

	dgeev_(&JOBVL,&JOBVR,&n,input,&n,eigReal,eigImag,
			vl,&n,vr,&n,work,&lwork,&info);
	ed::norm_ev(vr,n);
	// check for errors
	try {
		if (info!=0) throw runtime_error( "Error: dgeev returned error code ");
	} catch(const exception &ex) {std::cout << ex.what() << "\n";}

	delete [] work;
	return eigReal;
}

double* dsyev(double *mat, double *vr, int n) {
	char JOBZ='V',UPLO='U';
	int info, lwork = -1;
	double *eigReal = new double[n]{0};
	double *work;
	double wkopt;

	dsyev_(&JOBZ,&UPLO,&n,mat,&n,eigReal,&wkopt,&lwork,&info); // get lwork size
	lwork = (int)wkopt;
	work = new double[lwork]{0};
	dsyev_(&JOBZ,&UPLO,&n,mat,&n,eigReal,work,&lwork,&info);
	ed::norm_ev(mat,n);

	try {
		if (info!=0) throw runtime_error( "Error: dsyev returned error code ");
	} catch(const exception &ex) {std::cout << ex.what() << "\n";}

	delete [] work;
	return eigReal;
}