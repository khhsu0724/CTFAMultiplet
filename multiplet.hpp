#include "hilbert.hpp"
#include "diagonalize.hpp"
#include "gaunt.hpp"
#ifndef MULTIPLET
#define MULTIPLET

double calc_U(double* gaunt1, double* gaunt2, double* SC, int size);
void calc_coulomb(Hilbert& hilbs, double* SC);
double* CTmat(int l, double del);
double* CFmat(int l, double* CF);
void calc_CF(Hilbert& hilbs, double del, double* CF);
void calc_SO(Hilbert& hilbs, double lambda);
void calc_CV(Hilbert& hilbs, double* FG);
void calc_HYB(Hilbert& hilbs, double* SC);


#endif