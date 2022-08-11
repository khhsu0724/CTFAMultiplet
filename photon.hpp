#include <vector>
#include <complex>
#include "hilbert.hpp"
#ifndef PHOTON
#define PHOTON

struct PM {
	bool XAS = false;
	bool RIXS = false;
	bool PE = false; // Photo Emission
	bool eloss = true;
	std::string edge;
	std::vector<double> pvin, pvout;
	PM(): XAS(false), RIXS(false), PE(false), eloss(true) {
		pvin = vecd(3,0);
		pvout = vecd(3,0);
	};
	void set_eloss(bool el) {
		if (RIXS) eloss = el;
		else if (el) std::cout << "eloss will not be set if RIXS is turned off" << std::endl;
	}
};

std::complex<double> proj_pvec(int ml, const vecd& pvec);
void calc_ham(Hilbert& hilbs, std::vector<double*>& SC, double* FG, double* CF, double const& SO);
void occupation(Hilbert& hilbs);
void write_XAS(vecd const& aben, vecd const& intensity, std::string file_dir, bool print = true);
void XAS(std::string input_dir, std::vector<double*>& SC, double* FG, double* CF, 
			double SO, bool HYB, int nedos, const PM& pm);
void write_RIXS(vecd const& peaks, vecd const& ab, vecd const& em, std::string file_dir, 
				bool eloss, bool print = true);
void RIXS(std::string input_dir, std::vector<double*>& SC, double* FG, double* CF, 
			double SO, bool HYB, int nedos, const PM& pm);

#endif