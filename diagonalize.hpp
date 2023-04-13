#ifndef DIAGONALIZE
#define DIAGONALIZE
#include "matrix.hpp"
#include <omp.h>
#if defined __has_include && __has_include (<arpack/arpack.hpp>) 
#if __APPLE__  //Disable some warning from MacOS
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-extensions"
#include <arpack/arpack.hpp>
#pragma GCC diagnostic pop
#else 
#include <arpack/arpack.hpp>
#endif
#endif
// For sherlock...
#if defined __has_include && __has_include (<arpack-ng/arpack.hpp>) 
#include <arpack-ng/arpack.hpp>
#endif
// // For Cori...
#if defined __has_include && __has_include (<arpack.hpp>) 
#include <arpack.hpp>
#endif
// #if defined __has_include && __has_include (<arpack/arpackicb.h>) 
// #include <arpack/arpackicb.h>
// #endif

void ed_dgees(double *_mat, double *_eigvec, double* _eigReal, size_t n);
void ed_dsyev(double *_mat, double *_eigval, size_t n);
void ed_dsyevr(double *_mat, double *_eigvec, double* _eigReal, size_t n);

#ifdef __ARPACK_HPP__
#define ARPACK_TOL 1e-10
template <typename Real> 
void ed_dsarpack(Matrix<Real>* ham, Real *_eigvec, Real* _eigval, size_t n, size_t _nev) {
	std::cout << "Using Arpack" << std::endl;
	const a_uint N      = n;
	const a_uint nev    = (_nev < N) ? _nev : N;
	const a_uint ncv    = std::min(2*nev+1,N); // NCV size could be increased??
	const a_uint ldv    = N;
	const a_uint ldz    = N;
	const a_uint lworkl = ncv * (ncv + 8);
	const a_uint rvec   = 1;      // need eigenvectors

	const Real tol   = ARPACK_TOL; // small tol => more stable checks after EV computation.
	const Real sigma = 0.0;      // not referenced in this mode

	Real* resid = new Real[N];
	Real* V = new Real[ldv*ncv];
	Real* workd = new Real[3*N]{0};
	Real* workl = new Real[lworkl];
	a_int* select = new a_int[ncv];

	a_int iparam[11], ipntr[11];
	iparam[0] = 1;      // ishift
	iparam[2] = 3*N; // 3*N on input: maxit; on output: actual iteration
	iparam[3] = 1;      // NB, only 1 allowed
	iparam[6] = 1;      // mode

	a_int info = 0, ido = 0;
	do {
		arpack::saupd(ido, arpack::bmat::identity, N, arpack::which::smallest_algebraic, 
						nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, 
						workl, lworkl, info);
		ham->mvmult(workd+ipntr[0]-1,workd+ipntr[1]-1,N);
	} while (ido == 1 || ido == -1);

	// check info and number of ev found by arpack.
	if (info < 0 || iparam[4] < nev) { /*arpack may succeed to compute more EV than expected*/
		std::cout << "ERROR in saupd: iparam[4] " << iparam[4] << ", nev " << nev
		          << ", info " << info << std::endl;
		throw std::domain_error("Error inside ARPACK routines");
	}

	arpack::seupd(rvec, arpack::howmny::ritz_vectors, select, _eigval, _eigvec, 
					ldz, sigma,arpack::bmat::identity, N, arpack::which::smallest_algebraic, 
					nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, 
					workl, lworkl, info);
	if (info != 0) throw std::runtime_error("Error in seupd, info " + std::to_string(info));

	delete [] resid;
	delete [] V;
	delete [] workd;
	delete [] workl;
	delete [] select;
	return;
}
#else
template <typename Real> 
void ed_dsarpack(Matrix<Real>* ham, Real *_eigvec, Real* _eigval, size_t n, size_t& nev) {
	std::cout << "Arpack not available, switching to mkl/lapack" << std::endl;
	nev = n;
	auto _ham = ham->get_dense();
	ed_dsyevr(_ham,_eigvec,_eigval,n);
	delete [] _ham;
	_ham = nullptr;
	return;
};
#endif

class Block {
private:
	double Sz = 0, Jz = 0, K = 0;
	double *_ham, *_eig, *_eigvec; // Implicit pointer that faciliates diagonalization
public:
	size_t size, f_ind; // Accumulated first index of this Block
	size_t nev;			// Number of eigenvalues
	int diag_option;
	Matrix<double>* ham;
	uptrd eig, eigvec;
	std::vector<int> einrange;
	std::vector<ulli> rank; // rank keeps some rank information for faster hashing
	Block(double Sz, double Jz, double K, size_t size, size_t f_ind_in = 0,
		std::vector<ulli> rank_in = std::vector<ulli>(0)):
		Sz(Sz), Jz(Jz), K(K), size(size), f_ind(f_ind_in) {
		// Need to build something that prevents premature access of eig & eigvec
		if (size > std::sqrt(SIZE_MAX)) throw std::overflow_error("matrix too large");
		nev = size;
		rank = std::move(rank_in);
	};
	
	Block(Block&& blk) : Sz(blk.Sz), Jz(blk.Jz), K(blk.K), size(blk.size), f_ind(blk.f_ind), nev(blk.nev) {
		// This is a bit odd, should I copy or move?
		ham = std::move(blk.ham);
		eig = std::move(blk.eig);
		eigvec = std::move(blk.eigvec);
		_ham = std::move(blk._ham);
		_eig = std::move(blk._eig);
		_eigvec = std::move(blk._eigvec);
		rank = std::move(blk.rank);
	};

	~Block() {}
	double get_sz() {return Sz;};
	double get_jz() {return Jz;};
	double get_k() {return K;};	
	void malloc_ham(int diag_option) {
		// If matrix is too large, automatically use arpack
		if (size >= 2e5) {
			this->diag_option = 4;
			ham = new EZSparse<double>();
		} else if (size <= 1e4) {
			this->diag_option = 2;
			ham = new Dense<double>();
		} else {
			this->diag_option = diag_option;
			if (diag_option == 4) ham = new EZSparse<double>();
			else ham = new Dense<double>();
		}
		ham->malloc(size);
	};
	void diagonalize(size_t nev_in = 20) {
		// Options: 1. DGEES 2. DSYEVR (MRRR) 3. DSYEV (Shifted QR) 
		// 4. ARPACK (Lanczos)
		int option = this->diag_option;
		if (option == 1 || option == 2) {
			nev = size;
			_ham = ham->get_dense();
			_eig = new double[nev]{0};
			_eigvec = new double[size*nev]{0};
			if (option == 1) ed_dgees(_ham,_eigvec,_eig,size);
			else if (option == 2) ed_dsyevr(_ham,_eigvec,_eig,size);
			eigvec = return_uptr<double>(&_eigvec);
			eig = return_uptr<double>(&_eig);
			delete [] _ham;
			_ham = nullptr;
		} else if (option == 3) {
			nev = size;
			_ham = ham->get_dense();
			_eig = new double[size]{0};
			ed_dsyev(_ham,_eig,size);
			eigvec = return_uptr<double>(&_ham);
			eig = return_uptr<double>(&_eig);
		} else if (option == 4) {
			nev = std::min(nev_in,size/3);
			_eig = new double[nev]{0};
			_eigvec = new double[size*nev]{0};
			std::cout << "Arpack Routine" << std::endl;
			ed_dsarpack(ham,_eigvec,_eig,size,nev);
			eigvec = return_uptr<double>(&_eigvec);
			eig = return_uptr<double>(&_eig);
			ham->clear_mat();
		} else throw std::invalid_argument("invalid diagonalize method");
		return;
	}
	void init_einrange() {
		einrange = std::vector<int>(nev,0); // Hopefully there are more elegant solutions in the future
	}
};

#endif 