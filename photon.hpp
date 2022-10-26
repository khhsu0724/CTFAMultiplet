#include <vector>
#include <complex>
#include <chrono>
#include "hilbert.hpp"
#ifndef PHOTON
#define PHOTON

struct PM {
	bool XAS = false;
	bool RIXS = false;
	bool PE = false; // Photo Emission
	bool eloss = true;
	bool spin_flip = false;
	int nedos = 0;
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

std::string pol_str(const vecd& pvec);
std::complex<double> proj_pvec(int ml, const vecd& pvec);
void occupation(Hilbert& hilbs, const std::vector<bindex>& si);
std::vector<double> wvfnc_weight(Hilbert& hilbs, const std::vector<bindex>& si, 
								int ligNum = 3, bool print = false);
double effective_delta(Hilbert& hilbs, int ligNum = 3);
void state_composition(Hilbert& hilbs, const std::vector<bindex>& si, size_t top = 10);
void peak_occupation(Hilbert& hilbs, vecd const& peak_en, vecd const& energy, 
						vecd const& intensity, double ref_en = 0, std::string mode = "list");
void basis_overlap(Hilbert& GS, Hilbert& EX, bindex inds, std::vector<dcomp>& blap, 
					const PM& pm, bool pvout = false);
void write_XAS(vecd const& aben, vecd const& intensity, std::string file_dir, bool print = true);
void XAS(Hilbert& GS, Hilbert& EX, const PM& pm);
void write_RIXS(vecd const& peaks, vecd const& ab, vecd const& em, std::string file_dir, 
				bool eloss, bool print = true);
void RIXS(Hilbert& GS, Hilbert& EX, const PM& pm);

#endif