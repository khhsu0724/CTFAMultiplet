#include <iostream>
#include <cmath>
#include <utility>
#include <initializer_list>
#include <cstring>
#include <map>
#if __has_include("boost/functional/hash.hpp")
    #include <boost/functional/hash.hpp>
#elif __has_include("boost/functional/hash/hash.hpp")
    #include <boost/functional/hash/hash.hpp>
#endif
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
    	std::copy(l.begin(),l.end(),g);
    }

    // Copy constructor
    Gaunt_coeff(const Gaunt_coeff& other) {
        size = other.size;
        g = new double[size];
        std::copy(other.g, other.g+size, g);
    }

    // Copy assignment
    Gaunt_coeff& operator=(const Gaunt_coeff& other) {
        if (this != &other) {
            delete[] g;
            size = other.size;
            g = new double[size];
            std::copy(other.g, other.g + size, g);
        }
        return *this;
    }

    // Move constructor
    Gaunt_coeff(Gaunt_coeff&& other) noexcept {
        size = other.size;
        g = other.g;
        other.g = nullptr;
        other.size = 0;
    }

    // Move assignment
    Gaunt_coeff& operator=(Gaunt_coeff&& other) noexcept {
        if (this != &other) {
            delete[] g;
            size = other.size;
            g = other.g;
            other.g = nullptr;
            other.size = 0;
        }
        return *this;
    }

    // Destructor
    ~Gaunt_coeff() {
        delete[] g;
    }
};

bool operator==(Gaunt_coeff gc1, Gaunt_coeff gc2)
{
    return gc1.size == gc2.size && std::memcmp(gc1.g, gc2.g, gc1.size) == 0;
}

#ifdef BOOST_FUNCTIONAL_HASH_HASH_HPP
// Boost is likely not required 
size_t hash(Gaunt_coeff gc) 
{
    size_t h = 0;
    for (double* g = gc.g; g != gc.g + gc.size; ++g)
        boost::hash_combine(h,*g);
    return h;
}
#endif

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
											 {{ 1, 1}, {0,  sqrt(3.0/15), 0,  -3.0/sqrt(245)}},
											 {{-1,-1}, {0,  sqrt(3.0/15), 0,  -3.0/sqrt(245)}},
											 {{ 1, 0}, {0, -1.0/sqrt(15), 0,  sqrt(18.0/245)}},
											 {{-1, 0}, {0, -1.0/sqrt(15), 0,  sqrt(18.0/245)}},
											 {{ 0, 2}, {0, 			   0, 0,  sqrt(15.0/245)}},
											 {{ 0,-2}, {0, 			   0, 0,  sqrt(15.0/245)}},
											 {{ 0, 1}, {0, -sqrt(3.0/15), 0, -sqrt(24.0/245)}},
											 {{ 0,-1}, {0, -sqrt(3.0/15), 0, -sqrt(24.0/245)}},
											 {{ 0, 0}, {0,  2.0/sqrt(15), 0,  sqrt(27.0/245)}},
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
		if (abs(ml1) > l1 || abs(ml2) > l2) throw invalid_argument("ml > l, invalid quantum number");
		double *g = new double[l1+l2+1]{0};
		if (l1 == 2 && l2 == 2) return dd_coeff.at({ml1,ml2}).g;
		else if (l1 == 1 && l2 == 1) return pp_coeff.at({ml1,ml2}).g;
		else if (l1 == 0 && l2 == 0) return ss_coeff.at({ml1,ml2}).g;
		else if (l1 == 0 && l2 == 1) return sp_coeff.at({ml1,ml2}).g;
		else if (l1 == 0 && l2 == 2) return sd_coeff.at({ml1,ml2}).g;
		else if (l1 == 1 && l2 == 2) return pd_coeff.at({ml1,ml2}).g;
		else if (l1 == 1 && l2 == 0) {
			Gaunt_coeff gc = sp_coeff.at({ml2,ml1});
			for (int i = 0; i <= l1+l2; ++i) g[i] = gc.g[i] * pow(-1,ml2-ml1);
		} else if (l1 == 2 && l2 == 0) {
			Gaunt_coeff gc = sd_coeff.at({ml2,ml1});
			for (int i = 0; i <= l1+l2; ++i) g[i] = gc.g[i] * pow(-1,ml2-ml1);
		} else if (l1 == 2 && l2 == 1) {
			Gaunt_coeff gc = pd_coeff.at({ml2,ml1});
			for (int i = 0; i <= l1+l2; ++i) g[i] = gc.g[i] * pow(-1,ml2-ml1);
		}
		else throw invalid_argument("f-f gaunt coefficeints is not coded yet");
		return g;
	}
	catch(const exception &ex) {
		std::cout << ex.what() << "\n";
	}
	// Return Gaunt coefficient in Tesseral harmonics????
	return 0;
}