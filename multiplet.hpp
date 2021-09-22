#include "hilbert.hpp"
#include "diagonalize.hpp"
#include "gaunt.hpp"
#ifndef MULTIPLET
#define MULTIPLET

double calc_U(double* gaunt1, double* gaunt2, double* SC, int size);
void calc_coulomb(Hilbert& hilbs, double* SC);
double* CFmat(int l, double tenDQ, int coord);
double* CTmat(int l, double del);
void calc_CF(Hilbert& hilbs, double del, double tenDQ, int coord);
void calc_SO(Hilbert& hilbs, double lambda);
void calc_CV(Hilbert& hilbs, double* FG);
void calc_HYB(Hilbert& hilbs, double* SC);


#endif