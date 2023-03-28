#ifndef MATRIX
#define MATRIX
#include "helper.hpp"
#if defined __has_include && __has_include ("mkl.h") 
#include "mkl.h" // Can we auto-detect if mkl is installed
#include "mkl_lapacke.h"
#include "mkl_blas.h"
#endif

typedef std::unique_ptr<double[]> uptrd;
// Different Matrix class that holds Hamiltonian

// MKL & Lapack Functions for Matrix Vector Multiplication
#ifdef __INTEL_MKL__
#define lpk_int MKL_INT
inline void _symv(char* UPLO, int* N, double* alpha, double* A, int* LDA, double* x, 
		int* incx, double* beta, double* y, int* incy) {
	dsymv(UPLO,N,alpha,A,LDA,x,incx,beta,y,incy);
	return;
}
#else
#define lpk_int int
extern "C" {
	extern void dsymv_(char*,int*,double*,double*,int*,double*,int*,
						double*,double*,int*);
}

inline void _symv(char* UPLO, int* N, double* alpha, double* A, int* LDA, double* x, 
		int* incx, double* beta, double* y, int* incy) {
	dsymv_(UPLO,N,alpha,A,LDA,x,incx,beta,y,incy);
}
#endif

template<typename T> 
std::unique_ptr<T[]> return_uptr(T** arr) {
	std::unique_ptr<T[]> tmp(*arr);
	*arr = NULL;
	return tmp;
}

// Matrix Class
template<typename T> class Matrix {
public:
	Matrix() {};
	std::string mat_type;
	int size = 0;
	virtual void fill_mat(int lind, int rind, T elem) = 0;
	virtual void malloc(int size) = 0;
	virtual T* get_dense() = 0;
	virtual void mvmult(T* vec_in, T* vec_out, int N) = 0;
};

template<typename T> 
class Dense : public Matrix<T> {
public:
	Dense() {this->mat_type = "D";};
	void fill_mat(int lind, int rind, T elem) {
		ham[lind+this->size*rind] += elem;
	};
	void malloc(int size) {
		this->size = size;
		ham = std::make_unique<T[]>(size*size);
		return;
	};
	T* get_dense() {
		T* _ham = ham.release();
		ham = nullptr;
		return _ham;
	};
	void mvmult(T* vec_in, T* vec_out, int n) {
		T* _ham = this->get_dense();
		lpk_int N = n, lda = N, inc = 1;
		T alpha = 1, beta = 0;
		char UPLO = 'U';
		_symv(&UPLO,&N,&alpha,_ham,&lda,vec_in,&inc,&beta,vec_out,&inc);
		this->ham = return_uptr<double>(&_ham);
		return;
	};
private:
	std::unique_ptr<double[]> ham;
};

template <typename T> 
class EZSparse : public Matrix<T> {
public:
	EZSparse(std::string mat_type = "S") {this->mat_type = mat_type;};
	void fill_mat(int lind, int rind, T elem) {
		// Don't save symmetric matrix element
		if (this->mat_type == "S" && lind > rind) return;
		indexi.emplace_back(lind);
		indexj.emplace_back(rind);
		val.emplace_back(elem);
		return;
	};
	void malloc(int size) {
		this->size = size;
		indexi.reserve(100*size);
		indexj.reserve(100*size);
		val.reserve(100*size);
		return;
	};
	T* get_dense() {
		T* dense = new T[this->size*this->size]{0};
		// Symmetric matrix will return the upper triangle part
		for (size_t e = 0; e < val.size(); ++e)
			dense[indexi[e]+indexj[e]*this->size] += val.at(e);
		return dense;
	};
	void mvmult(T* vec_in, T* vec_out, int N) {
		#pragma omp parallel for shared(vec_out)
		for (int i = 0; i < this->size; ++i) vec_out[i] = 0;
		if (this->mat_type == "S") {
			#pragma omp parallel for shared(vec_out)
			for (size_t e = 0; e < val.size(); ++e) {
				if (indexj[e] == indexi[e])
					#pragma omp atomic update
					vec_out[indexj[e]] += val[e] * vec_in[indexi[e]];
				else {
					#pragma omp atomic update
					vec_out[indexj[e]] += val[e] * vec_in[indexi[e]];
					#pragma omp atomic update
					vec_out[indexi[e]] += val[e] * vec_in[indexj[e]];
				}
			}
		} else {
			#pragma omp parallel for shared(vec_out)
			for (size_t e = 0; e < val.size(); ++e) {
				#pragma omp atomic update
				vec_out[indexj[e]] += val[e] * vec_in[indexi[e]];
			}
		}
		return;
	};
private:
	std::vector<size_t> indexi;
	std::vector<size_t> indexj;
	std::vector<T> val;                    
};

// A wrapper class for boost Sparse Matrix
#endif