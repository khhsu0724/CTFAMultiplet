#ifndef DIAGONALIZE
#define DIAGONALIZE
#include "helper.hpp"
#include <memory>
#include <omp.h>
#if defined __has_include && __has_include ("mkl.h") 
#include "mkl.h" // Can we auto-detect if mkl is installed
#include "mkl_lapacke.h"
#endif


typedef std::unique_ptr<double[]> uptrd;

void ed_dgees(double *_mat, double *_eigvec, double* _eigReal, size_t n);
void ed_dsyev(double *_mat, double *_eigval, size_t n);
void ed_dsyevr(double *_mat, double *_eigvec, double* _eigReal, size_t n);

class Block {
private:
	double Sz = 0, Jz = 0, K = 0;
	double *_ham, *_eig, *_eigvec; // Implicit pointer that faciliates diagonalization
public:
	size_t size, f_ind; // Accumulated first index of this Block
	uptrd ham, eig, eigvec; // Use unique pointers that are auto-managed
	std::vector<int> einrange;
	std::vector<ulli> rank; // rank keeps some rank information for faster hashing
	Block(double Sz, double Jz, double K, size_t size, size_t f_ind_in = 0,
		std::vector<ulli> rank_in = std::vector<ulli>(0)): 
		Sz(Sz), Jz(Jz), K(K), size(size), f_ind(f_ind_in) {
		// Need to build something that prevents premature access of eig & eigvec
		if (size > std::sqrt(SIZE_MAX)) throw std::overflow_error("matrix too large");
		ham = std::make_unique<double[]>(size*size);
		rank = std::move(rank_in);
	};
	
	Block(Block&& blk) : Sz(blk.Sz), Jz(blk.Jz), K(blk.K), size(blk.size), f_ind(blk.f_ind) {
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
	void diagonalize(int option = 2) {
		_eig = new double[size]{0};
		_ham = ham.release(); // Release for diagonalize
		if (option < 3) {
			_eigvec = new double[size*size]{0};
			if (option == 1) ed_dgees(_ham,_eigvec,_eig,size);
			else if (option == 2) ed_dsyevr(_ham,_eigvec,_eig,size);
			else throw std::invalid_argument("invalid diagonalize method");
			eigvec = return_uptr(&_eigvec);
			eig = return_uptr(&_eig);
			delete [] _ham;
			_ham = nullptr;
		} else if (option == 3) {
			ed_dsyev(_ham,_eig,size);
			eigvec = return_uptr(&_ham);
			eig = return_uptr(&_eig);
		} else throw std::invalid_argument("invalid diagonalize method");
		return;
	}
	void init_einrange() {
		einrange = std::vector<int>(size,0); // Hopefully there are more elegant solutions in the future
	}
	uptrd return_uptr(double** arr) {
		uptrd tmp(*arr);
		*arr = NULL;
		return tmp;
	}
};

#endif 