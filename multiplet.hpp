#include "hilbert.hpp"
#include "diagonalize.hpp"
#include "gaunt.hpp"
#ifndef MULTIPLET
#define MULTIPLET

double calc_U(double* gaunt1, double* gaunt2, const double* SC, int size);
void calc_coulomb(Hilbert& hilbs, const std::vector<double*>& SC);
vecd CFmat(int l, const double* CF);
void calc_CF(Hilbert& hilbs, const double* CF);
void calc_SO(Hilbert& hilbs, const double lambda);
void calc_CV(Hilbert& hilbs, const double* FG);
vecd HYBmat(Hilbert& hilbs, const HParam& hparam);
void calc_HYB(Hilbert& hilbs, const HParam& hparam);


#endif