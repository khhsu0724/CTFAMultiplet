#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <tuple>
#include <utility>
#ifndef HELPER
#define HELPER

#define TOL 1E-7
typedef unsigned long long int ulli;

namespace ed {
	bool is_pw2(int x);
	int choose(int n, int k);
	ulli next_perm(ulli v);
	void enum_states(std::vector<ulli>& states, ulli n, ulli k, ulli inc = 0, ulli s = 0);
	ulli add_bits(ulli b1, ulli b2, int b1size, int b2size);
	int count_bits(ulli b);
	void sph2tet(double* sph, double* tet);

	template <typename T> T dot(std::vector<T> a, std::vector<T> b) {
		try {
			if (a.size() != b.size()) std::invalid_argument("different vector size for dot product");
			T dp = 0;
			for (int i = 0; i < a.size(); ++i) dp += a[i]*b[i];
			return dp;
		} catch (const std::exception &ex) {
			std::cout << ex.what() << "\n";
			exit(0);
		}
	};

	template <typename T> double norm(std::vector<T>& vin) {
		double n = 0;
		for (auto& v : vin) n += pow(abs(v),2);
		return sqrt(n);
	};

	template <typename T> void norm_vec(std::vector<T>& invec) {
		double norm = 0;
		for (auto & i : invec) norm += pow(abs(i),2);
		norm = sqrt(norm);
		for (auto & i : invec) i = i / norm;
	};

	template <typename T> void norm_ev(T* ev, int n) {
		for (int i = 0; i < n; ++i) {
			T nm = 0;
			for (int j = 0; j < n; j++) {
				// if (abs(ev[i*n+j]) < 1e-7) ev[i*n+j] = 0;
				nm += ev[i*n+j]*ev[i*n+j];
			}
			if (nm != 1) {
				nm = sqrt(nm);
				for (int j = 0; j < n; j++) ev[i*n+j] /= nm;
			}
		}
		return;
	}

	template <typename T> std::vector<T> intersection(std::vector<T>& v1, std::vector<T>& v2) {
		std::vector<T> v3;
	    std::set_intersection(v1.begin(),v1.end(),
	                          v2.begin(),v2.end(),
	                          back_inserter(v3));
	    return v3;
	};

	template <typename DT> void write_mat(DT* mat, size_t x, size_t y, std::string file_dir) {
		std::ofstream matfile;
	    matfile.open (file_dir);
		for (int j = 0; j < y; ++j) {
	    	for (int i = 0; i < x; ++i) {
	    		matfile << std::setw(8) << mat[i+x*j];
	    		if (i < x-1) matfile << " ";
	    	}
	    	matfile << "\n";
	    }
	    matfile.close();
	};

	template <typename DT> void write_vec(std::vector<DT> vec, size_t x, size_t y, std::string file_dir) {
		std::ofstream matfile;
	    matfile.open (file_dir);
		for (int j = 0; j < y; ++j) {
	    	for (int i = 0; i < x; ++i) {
	    		matfile << vec.at(i+x*j);
	    		if (i < x-1) matfile << ",";
	    	}
	    	matfile << "\n";
	    }
	    matfile.close();
	};



	template <typename T> std::vector<T> printDistinct(T arr[], int n, bool is_print = true) {
	    std::vector<T> unique_eig;
	    for (int i = 0; i < n; i++) {
	        int j;
	        for (j = 0; j < i; j++)
	           if (abs(arr[i] - arr[j]) < 1e-7)
	               break;
	        if (i == j) {
	        	if (is_print) std::cout << arr[i] << " ";
	        	unique_eig.push_back(arr[i]);
	        }
	    }
	    return unique_eig;
	};

	template <typename T> T trace(T *mat, int n) {
		T t = 0;
		for (int i = 0; i < n; ++i) t += mat[i+n*i];
		return t;
	};

	template <typename T> std::vector<T> vec_cross(std::vector<T> v1, std::vector<T> v2) {
		return std::vector<T>{v1[1]*v2[2]-v1[2]*v2[1],-v1[0]*v2[2]+v1[2]*v2[0],v1[0]*v2[1]-v1[1]*v2[0]};
	};

	template<typename T> std::vector<T> slice(std::vector<T> const &v, int m, int n) {
	    return std::vector<T> (v.begin()+m, v.begin()+n+1);;
	};

	template<typename T> std::vector<T> vec_mult(std::vector<T> v, T mult) {
		for (auto &ve : v) ve *= mult;
		return v;
	};
	// Add 2 vectors together
	template<typename T> void vec_add(std::vector<T>& v1, std::vector<T> v2) {
		if (v1.size() != v2.size()) return;
		for (size_t i = 0; i < v1.size(); ++i) v1[i] += v2[i];
	};

	template<typename T> bool veccmp(std::vector<T> v1, std::vector<T> v2) {
		if (v1.size() != v2.size()) return false;
		for (size_t i = 0; i < v1.size(); ++i) {
			if (v1[i] != v2[i]) return false;
		}
		return true;
	};

	template<typename T> int binary_search(std::vector<T> &svec, T t) {
		int ind = svec.size()/2 - 1, end = svec.size() - 1;
		while (abs(svec[ind]-t) > TOL) {
			if (ind == end) return -1;
			else if (svec[ind] < t) ind = (ind+end+1)/2;
			else if (svec[ind] > t) {
				end = ind;
				ind /= 2;
			}
		}
		return ind;
	};

	template<typename T> bool is_zero_arr(T* arr, int s) {
		for (int i = 0; i < s; ++i) if (abs(arr[i]) > TOL) return false;
		return true;
	};
}

#endif