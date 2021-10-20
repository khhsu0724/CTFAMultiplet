#include <vector>
#include <complex>
#include "hilbert.hpp"
#ifndef PHOTON
#define PHOTON

typedef std::vector<double> vecd;

std::complex<double> proj_pvec(int ml, vecd& pvec);
void calc_ham(Hilbert& hilbs, double* SC, double* FG, double const& tenDQ, double const& lambda);
void write_XAS(vecd const& peaks, vecd const& y, std::string file_dir, bool print = true);
void XAS(double* SC, double* FG, double tenDQ, double lambda, vecd& pvec, int nedos);

#endif