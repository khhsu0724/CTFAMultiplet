#ifndef GAUNT
#define GAUNT
	struct Gaunt_coeff;
	bool operator==(Gaunt_coeff gc1, Gaunt_coeff gc2);
	double*  gaunt(int,int,int,int);
#endif 