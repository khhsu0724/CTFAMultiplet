#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include "hilbert.h"
#ifndef SITE
#define SITE

using namespace std;

class UnitCell {
private:
public:
public:
	UnitCell() {};
	~UnitCell() {};
}

class Site {
private:
	int atnum;
	string atname;
public:
	double* pos;
	int ao_num;
	int t_hsize;
	int hsize, hindex;
	bool add_ve, rem_ce;
	vector<Hilbert> orbs;
public:
	Site() {};
	~Site() {};
	explicit Site(string atname, bool add_ve = false, bool rem_ce = false, int ao_num = 2);
	void build_index();
	void print_state(int pos);
	Hilbert& get_hspace(int n, int l);
	int qn_get_state(int n, int l, int pos);
	double* simplify_state(double* eigvec, int n, int l);
};

template <typename T> std::vector<T> intersection(std::vector<T>& v1, std::vector<T>& v2) {
	std::vector<T> v3;
    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          back_inserter(v3));
    return v3;
}

void cv_interacation(Site& core_site, Site& val_site, double* mat, double* FG);
void populate_hspace(Site* sites, int site_num, double* mat, double* SC, double tenDQ = 0, double lambda = 0);
void hybrid();

#endif