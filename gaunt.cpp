#include <iostream>
#include <cmath>
#include <utility>
#include <initializer_list>
#include <map>
#include <boost/functional/hash.hpp>
#include "gaunt.hpp"

using namespace std;

struct Gaunt_coeff
{
    size_t size;
    double* g;
    Gaunt_coeff(size_t size, double* g) {
    	this->size = size;
    	this->g = g;
    }
    Gaunt_coeff(const initializer_list<double> l) {
    	this->size = l.size();
    	this->g = new double[l.size()];
    	copy(l.begin(),l.end(),g);
    }
    ~Gaunt_coeff() {} //Figure out of destructor is needed
};

bool operator==(Gaunt_coeff gc1, Gaunt_coeff gc2)
{
    return gc1.size == gc2.size && memcmp(gc1.g, gc2.g, gc1.size) == 0;
}


size_t hash(Gaunt_coeff gc) 
{
    size_t h = 0;
    for (double* g = gc.g; g != gc.g + gc.size; ++g)
        boost::hash_combine(h,*g);
    return h;
}

map<pair<int,int>,Gaunt_coeff> dd_coeff =   {{{ 2, 2}, {1, 0,     -2.0/7.0, 0,       1.0/21}},
											 {{-2,-2}, {1, 0,     -2.0/7.0, 0,       1.0/21}},
											 {{ 2, 1}, {0, 0,  sqrt(6)/7.0, 0,  -sqrt(5)/21}},
											 {{-2,-1}, {0, 0,  sqrt(6)/7.0, 0,  -sqrt(5)/21}},
											 {{ 1, 2}, {0, 0, -sqrt(6)/7.0, 0,   sqrt(5)/21}}, //
											 {{-1,-2}, {0, 0, -sqrt(6)/7.0, 0,   sqrt(5)/21}}, //
											 {{ 2, 0}, {0, 0,     -2.0/7.0, 0,  sqrt(15)/21}},
											 {{-2, 0}, {0, 0,     -2.0/7.0, 0,  sqrt(15)/21}},
											 {{ 0, 2}, {0, 0,     -2.0/7.0, 0,  sqrt(15)/21}},
											 {{ 0,-2}, {0, 0,     -2.0/7.0, 0,  sqrt(15)/21}},
											 {{ 1, 1}, {1, 0,      1.0/7.0, 0,      -4.0/21}},
											 {{-1,-1}, {1, 0,      1.0/7.0, 0,      -4.0/21}},
											 {{ 1, 0}, {0, 0,      1.0/7.0, 0,  sqrt(30)/21}},
											 {{-1, 0}, {0, 0,      1.0/7.0, 0,  sqrt(30)/21}},
											 {{ 0, 1}, {0, 0,     -1.0/7.0, 0, -sqrt(30)/21}}, //
											 {{ 0,-1}, {0, 0,     -1.0/7.0, 0, -sqrt(30)/21}}, // 
											 {{ 0, 0}, {1, 0,      2.0/7.0, 0,       6.0/21}},
											 {{ 2,-2}, {0, 0,            0, 0,  sqrt(70)/21}},
											 {{-2, 2}, {0, 0,            0, 0,  sqrt(70)/21}},
											 {{ 2,-1}, {0, 0,            0, 0, -sqrt(35)/21}},  
											 {{-2, 1}, {0, 0,            0, 0, -sqrt(35)/21}}, 
											 {{-1, 2}, {0, 0,            0, 0,  sqrt(35)/21}}, // 
											 {{ 1,-2}, {0, 0,            0, 0,  sqrt(35)/21}}, //
											 {{ 1,-1}, {0, 0, -sqrt(6)/7.0, 0, -sqrt(40)/21}},
											 {{-1, 1}, {0, 0, -sqrt(6)/7.0, 0, -sqrt(40)/21}},
											};

map<pair<int,int>,Gaunt_coeff> pp_coeff =   {{{ 1, 1}, {1, 0,     -1.0/5.0}},
											 {{-1,-1}, {1, 0,     -1.0/5.0}},
											 {{ 1, 0}, {0, 0,  sqrt(3)/5.0}},
											 {{-1, 0}, {0, 0,  sqrt(3)/5.0}},
											 {{ 0, 1}, {0, 0, -sqrt(3)/5.0}}, //
											 {{ 0,-1}, {0, 0, -sqrt(3)/5.0}}, //
											 {{ 0, 0}, {1, 0,      2.0/5.0}},
											 {{ 1,-1}, {0, 0, -sqrt(6)/5.0}},
											 {{-1, 1}, {0, 0, -sqrt(6)/5.0}},
											};

map<pair<int,int>,Gaunt_coeff> ss_coeff =  	{{{ 0, 0}, {1}},
											};

map<pair<int,int>,Gaunt_coeff> sp_coeff =  	{{{ 0, 0}, {0,  1/sqrt(3)}},
											 {{ 0, 1}, {0, -1/sqrt(3)}},
											 {{ 0,-1}, {0, -1/sqrt(3)}},
											};

map<pair<int,int>,Gaunt_coeff> sd_coeff =  	{{{ 0, 0}, {0, 0,  1/sqrt(5)}},
											 {{ 0, 1}, {0, 0, -1/sqrt(5)}},
											 {{ 0,-1}, {0, 0, -1/sqrt(5)}},
											 {{ 0, 2}, {0, 0,  1/sqrt(5)}},
											 {{ 0,-2}, {0, 0,  1/sqrt(5)}},
											};

map<pair<int,int>,Gaunt_coeff> pd_coeff =  	{{{ 1, 2}, {0, -sqrt(6.0/15), 0,   sqrt(3.0/245)}},
											 {{-1,-2}, {0, -sqrt(6.0/15), 0,   sqrt(3.0/245)}},
											 {{ 1, 1}, {0,  sqrt(3.0/15), 0,  -sqrt(9.0/245)}},
											 {{-1,-1}, {0,  sqrt(3.0/15), 0,  -sqrt(9.0/245)}},
											 {{ 1, 0}, {0, -sqrt(1.0/15), 0,  sqrt(18.0/245)}},
											 {{-1, 0}, {0, -sqrt(1.0/15), 0,  sqrt(18.0/245)}},
											 {{ 0, 2}, {0, 			   0, 0,  sqrt(15.0/245)}},
											 {{ 0,-2}, {0, 			   0, 0,  sqrt(15.0/245)}},
											 {{ 0, 1}, {0, -sqrt(3.0/15), 0, -sqrt(24.0/245)}},
											 {{ 0,-1}, {0, -sqrt(3.0/15), 0, -sqrt(24.0/245)}},
											 {{ 0, 0}, {0,  sqrt(4.0/15), 0,  sqrt(27.0/245)}},
											 {{ 1,-2}, {0, 		  	   0, 0,  sqrt(45.0/245)}},
											 {{-1, 2}, {0,             0, 0,  sqrt(45.0/245)}},
											 {{ 1,-1}, {0, 			   0, 0, -sqrt(30.0/245)}},
											 {{-1, 1}, {0,             0, 0, -sqrt(30.0/245)}},
											};

double* gaunt(int l1, int ml1, int l2, int ml2) {
	// Return Gaunt coefficient in spherical harmonics
	// Input: quantumn numbers, Output: pointer to the Gaunt coefficient value array
	try {
		// catch error if ml1 > l1 or ml2 > l2
		if (ml1 > l1 || ml2 > l2) {
			throw invalid_argument("ml > l, invalid quantum number");
		}
		// d-d case
		if (l1 == 2 && l2 == 2) return dd_coeff.at({ml1,ml2}).g;
		else if (l1 == 1 && l2 == 1) return pp_coeff.at({ml1,ml2}).g;
		else if (l1 == 0 && l2 == 0) return ss_coeff.at({ml1,ml2}).g;
		else if (l1 == 0 && l2 == 1) return sp_coeff.at({ml1,ml2}).g;
		else if (l1 == 0 && l2 == 2) return sd_coeff.at({ml1,ml2}).g;
		else if (l1 == 1 && l2 == 2) return pd_coeff.at({ml1,ml2}).g;
		else if (l1 == 1 && l2 == 0) {
			double *g;
			for (int i = 0; i <= l1+l2; ++i) g[i] = sp_coeff.at({ml2,ml1}).g[i] * pow(-1,abs(ml1-ml2));
			return g;
		} else if (l1 == 2 && l2 == 0) {
			double *g;
			for (int i = 0; i <= l1+l2; ++i) g[i] = sd_coeff.at({ml2,ml1}).g[i] * pow(-1,abs(ml1-ml2));
			return g;
		} else if (l1 == 2 && l2 == 1) {
			double *g;
			for (int i = 0; i <= l1+l2; ++i) g[i] = pd_coeff.at({ml2,ml1}).g[i] * pow(-1,abs(ml1-ml2));
			return g;
		}
		else {
			throw invalid_argument("f-f gaunt coefficeints is not coded yet");
		}
	}
	catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	// Return Gaunt coefficient in Tesseral harmonics????
	return 0;
}