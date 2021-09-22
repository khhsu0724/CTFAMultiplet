#include <vector>
#include <complex>
#include "hilbert.hpp"
#ifndef PHOTON
#define PHOTON

std::complex<double> proj_pvec(int ml, std::vector<double>& pvec);
void calc_ham(Hilbert& hilbs, double* SC, double* FG, double const& tenDQ, double const& lambda);
void XAS(double* SC, double* FG, double tenDQ, double lambda, std::vector<double>& pvec, int nedos);

#endif