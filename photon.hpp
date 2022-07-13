#include <vector>
#include <complex>
#include "hilbert.hpp"
#ifndef PHOTON
#define PHOTON

std::complex<double> proj_pvec(int ml, vecd& pvec);
void calc_ham(Hilbert& hilbs, std::vector<double*>& SC, double* FG, double* CF, double const& SO);
void write_XAS(vecd const& aben, vecd const& intensity, std::string file_dir, bool print = true);
void XAS(std::string input_dir, std::vector<double*>& SC, double* FG, double* CF, 
			double SO, bool HYB, vecd& pvec, int nedos, std::string edge);
void write_RIXS(vecd const& peaks, vecd const& ab, vecd const& em, std::string file_dir, 
				bool rel_em = true, bool print = true);
void RIXS(std::string input_dir, std::vector<double*>& SC, double* FG, double* CF, 
			double SO, bool HYB, vecd& pvin, vecd& pvout, int nedos, std::string edge);

#endif