#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <utility>
#include "hilbert.h"
#include "site.h"
#ifndef PHOTON
#define PHOTON

template <typename T> void norm_vec(std::vector<T>& invec) {
	double norm = 0;
	for (auto & i : invec) norm += pow(abs(i),2);
	norm = sqrt(norm);
	for (auto & i : invec) i = i / norm;
};

template <typename DT> void write_mat(DT* mat, size_t n, size_t m, std::string file_dir) {
	std::ofstream matfile;
    matfile.open (file_dir);
	for (int j = 0; j < m; ++j) {
    	for (int i = 0; i < n; ++i) {
    		matfile << mat[i+n*j];
    		if (i < n-1) matfile << " ";
    	}
    	matfile << endl;
    }
    matfile.close();
};

bool is_pw2(int x);
std::vector<double> printDistinct(double arr[], int n, bool is_print = true);
std::complex<double> proj_pvec(int ml, std::vector<double>& pvec);
void XAS(double* SC, double* FG, double tenDQ, double lambda, std::vector<double>& pvec, int nedos);

#endif