#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include "hilbert.h"
#include "diagonalize.h"
#include "gaunt.h"
#ifndef MULTIPLET
#define MULTIPLET

double calc_U(double* gaunt1, double* gaunt2, double* SC, int size);
void calc_coulomb(Hilbert& hspace, double* mat, double* SC);
double* CFmat(int l, double tenDQ, int coord);
void calc_CF(Hilbert& hspace, double* mat, double tenDQ, int coord);
void calc_SO(Hilbert& hspace, double* mat, double lambda);

#endif