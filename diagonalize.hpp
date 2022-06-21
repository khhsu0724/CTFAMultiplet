#ifndef DIAGONALIZE
#define DIAGONALIZE
#include <memory>
#include "helper.hpp"

typedef std::unique_ptr<double[]> uptrd;

void dgees(double *_mat, double *_eigvec, double* _eigReal, int n);
void dsyev(double* _mat, double* _eigval, int n);

class Block {
private:
	double Sz = 0, Lz = 0, K = 0;
	double *_ham, *_eig, *_eigvec; // Implicit pointer that faciliates diagonalization
public:
	Block(double Sz, double Lz, double K, int size): Sz(Sz), Lz(Lz), K(K), size(size) {
		// We can preallocate this with a number. i.e. Hilbert space size/number of blocks * 1.2
		// Also need to build something that prevents premature access of eig & eigvec
		ham = std::make_unique<double[]>(size*size);
	};
	Block(Block&& blk) : Sz(blk.Sz), Lz(blk.Lz), K(blk.K), size(blk.size) {
		ham = std::move(blk.ham);
		eig = std::move(blk.eig);
		eigvec = std::move(blk.eigvec);
		_ham = std::move(blk._ham);
		_eig = std::move(blk._eig);
		_eigvec = std::move(blk._eigvec);
	};

	~Block() {}
	double get_Sz() {return Sz;};
	double get_Lz() {return Lz;};
	double get_K() {return K;};	
	int size;
	uptrd ham, eig, eigvec; // Use unique pointers that are auto-managed
	std::vector<int> einrange;
	void diag_dgees() {
		_ham = ham.release();
		_eig = new double[size]{0};
		_eigvec = new double[size*size]{0};
		dgees(_ham,_eigvec,_eig,size);
		eigvec = return_uptr(&_eigvec);
		eig = return_uptr(&_eig);
		delete [] _ham;
		_ham = nullptr;
		return;
	}
	void diag_dsyev() {
		// Only eigvec, eig will point to something at the end
		_eig = new double[size]{0};
		_ham = ham.release(); // Release for diagonalize
		dsyev(_ham,_eig,size);
		eigvec = return_uptr(&_ham);
		eig = return_uptr(&_eig);
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