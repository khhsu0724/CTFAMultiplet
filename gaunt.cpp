#include <iostream>
#include <cmath>
#include <utility>
#include <initializer_list>
#include <map>
#include <boost/functional/hash.hpp>
#include "gaunt.h"

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

map<pair<int,int>,Gaunt_coeff> dd_coeff =  { {{ 2, 2}, {1, 0,     -2.0/7.0, 0,       1.0/21}},
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

double* gaunt(int l1, int ml1, int l2, int ml2) {
	// Return Gaunt coefficient in spherical harmonics
	// Input: quantumn numbers, Output: pointer to the Gaunt coefficient value array
	try {
		// catch error if ml1 > l1 or ml2 > l2
		if (ml1 > l1 || ml2 > l2) {
			throw invalid_argument("ml > l, invalid quantum number");
		}
		// d-d case
		if (l1 == 2 && l2 ==2) return dd_coeff.at({ml1,ml2}).g;
		else {
			throw invalid_argument("anything beside d-d gaunt coefficeints is not coded yet");
		}
	}
	catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	// Return Gaunt coefficient in Tesseral harmonics????
	return 0;
}