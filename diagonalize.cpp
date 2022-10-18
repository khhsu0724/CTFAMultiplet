#include <iostream> 
#include "diagonalize.hpp"
using namespace std;

// Lapack functions, need dsyevx or dsyevr
// Rerference: https://epubs.stfc.ac.uk/manifestation/1885/DLTR-2007-005.pdf
// Scalapack functions
// extern void pdsyevd_( char* jobz, char* uplo, int* n, double* a, int* ia, int* ja, int* desca,
//                 double* w, double* z, int* iz, int* jz, int* descz, double* work,
//                 int* lwork, double* iwork, int* liwork, int* info );

#ifdef __INTEL_MKL__

void ed_dgees(double *_mat, double *_eigvec, double* _eigReal, size_t n) {
// 	double *eigImag = new double[n]{0};
// 	int info = 0, sdim = 0, lwork = 6*n;
// 	int *bwork;
// 	double *work = new double[lwork]{0};
// 	int (*select)(const double*,const double*);
// 	dgees_('V','N',select,&n,_mat,&n,&sdim,_eigReal,eigImag,
// 			_eigvec,&n,work,&lwork,bwork,&info);
// 	try {
// 		if (info!=0) throw runtime_error( "Error: dgees returned error code ");
// 	} catch(const exception &ex) {std::cout << ex.what() << "\n";}
// 	delete [] work;
	return;
}

void ed_dsyev(double* _mat, double* _eigval, size_t n) {
	MKL_INT N = n, lda = N, info, lwork;
	double wkopt;
	double* work;
	lwork = -1;
	dsyev("Vectors","Upper",&N,_mat,&N,_eigval,&wkopt,&lwork,&info);
	lwork = (MKL_INT)wkopt;
	work = (double*)malloc(lwork*sizeof(double));
	dsyev("Vectors","Upper",&N,_mat,&N,_eigval,work,&lwork,&info);
	try {
		if (info!=0) throw runtime_error( "Error: dsyev returned error code ");
	} catch(const exception &ex) {std::cout << ex.what() << "\n";}
	delete [] work;
	return;
}

void ed_dsyevr(double* _mat, double *_eigvec, double* _eigval, size_t n) {
	MKL_INT N = n, il, iu, info, lwork = -1, liwork = -1, iwkopt;
	double abstol = -1.0, vl, vu, wkopt;
	double* work;
	MKL_INT* iwork;
	MKL_INT* isuppz = new MKL_INT[2*N];
	dsyevr("Vectors","All","Upper",&N,_mat,&N,&vl,&vu,&il,&iu,&abstol,&N,
			_eigval,_eigvec,&N,isuppz,&wkopt,&lwork,&iwkopt,&liwork,&info);
	lwork = (MKL_INT)wkopt;
	work = (double*)malloc(lwork*sizeof(double));
	liwork = iwkopt;
	iwork = (MKL_INT*)malloc(liwork*sizeof(MKL_INT));
	N = n; // I'm not sure why I have to do this
	dsyevr("Vectors","All","Upper",&N,_mat,&N,&vl,&vu,&il,&iu,&abstol,&N,
			_eigval,_eigvec,&N,isuppz,work,&lwork,iwork,&liwork,&info);
	try {
		if (info!=0) throw runtime_error( "Error: dsyev returned error code ");
	} catch(const exception &ex) {std::cout << ex.what() << "\n";}
	delete [] work;
	delete [] iwork;
	return;
}

#else
extern "C" {
	extern void dgees_(char*,char*,int(*)(double*,double*),size_t*,double*,size_t*,
						int*,double*,double*,double*,size_t*,double*,int*,int*,int*);
}
extern "C" {
	extern void dsyev_(char*,char*,size_t*,double*,size_t*,double*,double*,int*,int*);
}

void ed_dgees(double *_mat, double *_eigvec, double* _eigReal, size_t n) {
	// Can also use for diagonlizing non-symmetric Hamiltonians
	char JOBVS='V',SORT='N';
	double *eigImag = new double[n]{0};
	int info = 0, sdim = 0, lwork = 6*n;
	int *bwork;
	double *work = new double[lwork]{0};
	int (*select)(double*,double*);
	dgees_(&JOBVS,&SORT,select,&n,_mat,&n,&sdim,_eigReal,eigImag,
			_eigvec,&n,work,&lwork,bwork,&info);
	try {
		if (info!=0) throw runtime_error( "Error: dgees returned error code ");
	} catch(const exception &ex) {std::cout << ex.what() << "\n";}
	delete [] work;
	return;
}

void ed_dsyev(double* _mat, double* _eigval, size_t n) {
	int info = 0, lwork = -1;
	char JOBZ='V',UPLO='U';
	double *work;
	double wkopt;
	dsyev_(&JOBZ,&UPLO,&n,_mat,&n,_eigval,&wkopt,&lwork,&info); // get lwork size
	lwork = (int)wkopt;
	work = new double[lwork]{0};
	dsyev_(&JOBZ,&UPLO,&n,_mat,&n,_eigval,work,&lwork,&info); // Should return orthonormal ev
	try {
		if (info!=0) throw runtime_error( "Error: dsyev returned error code ");
	} catch(const exception &ex) {std::cout << ex.what() << "\n";}
	delete [] work;
	return;
}

// TODO
void ed_dsyevr(double* _mat, double *_eigvec, double* _eigval, size_t n) {return;}

#endif

