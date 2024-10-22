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
// For HPC self installed arpack
#if defined __has_include && __has_include (<arpack.hpp>) 
#include <arpack.hpp>
#endif

#pragma omp declare reduction(+:dcomp:omp_out=omp_out+omp_in) initializer (omp_priv=omp_orig)

void ed_dgees(double *_mat, double *_eigvec, double* _eigReal, size_t n);
void ed_dsyev(double *_mat, double *_eigval, size_t n);
void ed_dsyevr(double *_mat, double *_eigvec, double* _eigReal, size_t n);

#ifdef __ARPACK_HPP__
#define ARPACK_TOL 1e-10
template <typename Real> 
void ed_dsarpack(Matrix<Real>* ham, Real *_eigvec, Real* _eigval, size_t n, size_t _nev) {
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
	iparam[2] = 3*N; 	// 3*N on input: maxit; on output: actual iteration
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

// std::complex<double> complex_sum(std::complex<double> a, std::complex<double> b) {
//     return a + b;
// }

template <typename T>
int Lanczos(Matrix<T>* ham, const vecc& v0, vecc& alpha, vecc& betha, int niter_CFE=150) {
	// returns new niter_CFE value
	int hsize = ham->get_mat_dim();
	// std::cout << "Hdim: " << hsize  << "," << v0.size() << std::endl;
	vecc phip(hsize,0);
	vecc phipp(hsize,0);
	vecc phil(hsize,0);
	// Normalize the vector
	vecc phi = v0;
	ed::norm_vec(phi);

	for (int i = 0; i < niter_CFE; ++i) {
		ham->mvmult_cmplx(phi,phip);
		dcomp alpha_element(0);
		#pragma omp parallel for reduction (+:alpha_element)
		for (int j = 0; j < hsize; ++j) alpha_element += std::conj(phi[j])*phip[j];
		alpha[i] = alpha_element;
		if (i != 0) {
			#pragma omp parallel for 
			for (int j = 0; j < hsize; ++j) phip[j] -= betha[i] * phil[j];
		}
		#pragma omp parallel for 
		for (int j = 0; j < hsize; ++j) phipp[j] = phip[j] - alpha[i] * phi[j];
		if (i != niter_CFE-1) {
			betha[i+1] = ed::norm(phipp);
			if (std::abs(betha[i+1]) < 1E-13) {
				std::cout << "Lanzcos ended after: " << i << " steps." << std::endl;
				niter_CFE = i;
				break;
			}
			phil = phi;
			#pragma omp parallel for
			for (int j = 0; j < hsize; ++j) phi[j] = phipp[j]/betha[i+1];
		}
	}
	return niter_CFE;
};

template <typename T>
void ContFracExpan(Matrix<T>* ham, const vecc& v0, double E0, vecd& specX, vecd& specY, 
					double eps = 0.1, int niter_CFE=150) {
	// This solver specifically does not subtract elastic scattering
	// specX should specify: minE, maxE, nedos
	int nedos = specX.size();
	vecc intensity(nedos,0);
	double factor = ed::norm(v0);
	factor = factor*factor;
	// std::cout << "FACTOR: " << factor << std::endl;
	vecc alpha(niter_CFE,0);
	vecc betha(niter_CFE,0);
	niter_CFE = Lanczos(ham,v0,alpha,betha,niter_CFE);
	// for (int i = 0; i < niter_CFE; ++i) std::cout << "a: " << alpha[i] << ", b: " << betha[i] << std::endl;
	#pragma omp parallel for shared(intensity,specY,alpha,betha)
	for (int i = 0; i < nedos; ++i) {
		dcomp z = dcomp(specX[i]+E0,eps);
		intensity[i] = z - alpha[niter_CFE-1];
		for (int j = 1; j < niter_CFE; ++j) {
			intensity[i] = z - alpha[niter_CFE-j-1] - 
							(pow(betha[niter_CFE-j],2)/intensity[i]);
		}
		specY[i] += -1/PI * std::imag(factor/intensity[i]);
	}

	return;
}

template <typename T>
vecc BiCGS(Matrix<T>* ham, const vecc& b, const dcomp z, double CG_tol = 1e-8, int max_iter = 2e4) {
	// Modified from Eigen, with diagonal preconditioner 1/(ham-z)
	// See https://github.com/cryos/eigen/blob/master/Eigen/src/IterativeLinearSolvers/BiCGSTAB.h
	int hsize = ham->get_mat_dim();
	vecc r(hsize,0), r0(hsize,0), x(hsize,0), x0(hsize,0);
	// for (int j = 0; j < hsize; ++j) // Compute r0 = b - Ax, here x is started with 0
	// 	r[j] = b[j];
	r = b;
	r0 = r;
	dcomp omega = 10.0, w_old = 0;
	dcomp rho = 1.0, alpha = 1.0, rho_old = 1.0;
    vecc v(hsize,0), t(hsize,0), h(hsize,0);
    vecc s(hsize,0), p(hsize,0), shat(hsize,0), phat(hsize,0);
    double r0sqnorm = ed::norm(r0,true);
    int i = 0;
    while (ed::norm(r,true) > CG_tol and i < max_iter) {
    	// STEP 1: rho = (r0,r)
    	omega = w_old; // This is necessary, but I don't know why
    	if (i%50 == 0) {
	    	std::cout << std::setprecision(15);
	        std::cout << "iter: " << i << ", redisual: " << ed::norm(r,true);
	        if (i!=0) std::cout << ", s:" << ed::norm(s,true);
	        std::cout << std::endl;
	    }
    	rho  = 0;
    	#pragma omp parallel for reduction (+:rho)
		for (int j = 0; j < hsize; ++j) rho += std::conj(r0[j])*r[j];
		// if (std::real(rho) < 1e-12) {
		// 	std::cout << "New search direction cannot be found in BiCGStab" << std::endl;
		// 	break;
		// }
		// std::cout << "rho: " << rho << std::endl;
		if (i > 0) {
			// STEP 2: beta = (rho_i/rho_i-1)(alpha/w)
			dcomp beta = (rho/rho_old) * (alpha/omega);
			// std::cout << "beta: " << beta << std::endl;
			// STEP 3: p = r + beta(p-wv)
			#pragma omp parallel for
			for (int j = 0; j < hsize; ++j)
	            p[j] = r[j] + beta * (p[j] - omega * v[j]);
	    } else {
	    	s = vecc(hsize,0);
	    	p = r;
	    }
        // STEP 4: v = A*p (adding precondition)
        phat = ham->precond(p,z);
        ham->mvmult_cmplx(phat,v);  
        #pragma omp parallel for
        for (int j = 0; j < hsize; ++j)
        	v[j] -= phat[j]*z;
        // STEP 5: alpha = rho/(r0,v)
        dcomp rv = 0;
        #pragma omp parallel for reduction (+:rv)
		for (int j = 0; j < hsize; ++j) rv += std::conj(r0[j])*v[j];
        alpha = rho/rv;
		// std::cout << "alpha: " << alpha << std::endl;
    	// STEP 6: h = x + alpha*p
    	// #pragma omp parallel for
		// for (int j = 0; j < hsize; ++j)
        //     h[j] = x[j] + alpha * p[j];
        // STEP 7: s = r - alpha*v
    	#pragma omp parallel for
		for (int j = 0; j < hsize; ++j)
            r[j] -= alpha * v[j];
        s = r;
        if (ed::norm(s,true) <= CG_tol) {
        	// Testing
        	#pragma omp parallel for
	        for (int j = 0; j < hsize; ++j)
	        	x[j] += alpha*phat[j];
        	std::cout << "s converged" << std::endl;
        	break;
        }
        // STEP 8: t = A*s (adding preconditioning)
        shat = ham->precond(s,z);
        ham->mvmult_cmplx(shat,t);  
        #pragma omp parallel for
        for (int j = 0; j < hsize; ++j)
        	t[j] -= shat[j]*z;
        // STEP 9: w = (t,s)/(t,t)
        dcomp omega = 0;
        #pragma omp parallel for reduction (+:omega)
		for (int j = 0; j < hsize; ++j) omega += std::conj(t[j])*s[j];
		omega /= ed::norm(t,true);
		// std::cout << "omega: " << omega << std::endl;
		// STEP 10: x = h + w*s
		#pragma omp parallel for
        for (int j = 0; j < hsize; ++j)
        	x[j] += alpha * phat[j] + omega * shat[j];
        	// x[j] = h[j] + w * s[j];
        // STEP 11: r = s -w*t
        #pragma omp parallel for
        for (int j = 0; j < hsize; ++j)
        	r[j] -= omega * t[j];
        rho_old = rho;
        w_old = omega;
        ++i;
    }
    std::cout << "CG finished with " << i << " iterations" << std::endl;
    std::cout << "Ridisual: " << ed::norm(r,true) << std::endl;
    if (i >= max_iter) std::cout << "CG not converged: " << ed::norm(r,true) << std::endl;
    return x;
}

template <typename T>
class Block  {
private:
	double Sz = 0, Jz = 0, K = 0;
	double *_ham, *_eig, *_eigvec, *_ham_copy; // Implicit pointer that faciliates diagonalization
public:
	size_t size, f_ind; // Accumulated first index of this Block
	size_t nev;			// Number of eigenvalues
	int diag_option;

	Matrix<T>* ham;
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
	
	Block(Block<T>&& blk) : Sz(blk.Sz), Jz(blk.Jz), K(blk.K), size(blk.size), f_ind(blk.f_ind), nev(blk.nev) {
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
		if (size >= 1e5) {
			this->diag_option = 4;
			ham = new EZSparse<double>();
		} else if (size <= 2e3) {
			this->diag_option = 2;
			ham = new Dense<double>();
		} else {
			if (diag_option > 4 && diag_option < 1) diag_option = 4;
			this->diag_option = diag_option;
			if (diag_option == 4) ham = new EZSparse<double>();
			else ham = new Dense<double>();
		}
		ham->malloc(size);
	};
	void diagonalize(size_t nev_in = 20, bool clear_mat = true) {
		// Options: 1. DGEES 2. DSYEVR (MRRR) 3. DSYEV (Shifted QR) 
		// 4. ARPACK (Lanczos)
		int option = this->diag_option;
		if (option == 1 || option == 2) {
			nev = size;
			_ham = ham->get_dense();
			if (!clear_mat) {
				// lapack deletes hamiltonian, keeping it
				_ham_copy = new double[size*size]{0};
				#pragma omp parallel for 
				for (int i = 0; i < size*size; ++i) _ham_copy[i] = _ham[i];
				ham->reset_ham(&_ham_copy);
			}
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
			std::cout << "Arpack Eigenvalues " << nev << ", size: " << size << std::endl;
			_eig = new double[nev]{0};
			_eigvec = new double[size*nev]{0};
			ed_dsarpack(ham,_eigvec,_eig,size,nev);
			eigvec = return_uptr<double>(&_eigvec);
			eig = return_uptr<double>(&_eig);
			if (clear_mat) ham->clear_mat();
		} else throw std::invalid_argument("invalid diagonalize method");
		return;
	}
	void init_einrange() {
		einrange = std::vector<int>(nev,0); // Hopefully there are more elegant solutions in the future
	}
};

#endif 
