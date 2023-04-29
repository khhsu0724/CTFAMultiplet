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
	double em_energy = 15;
	std::string edge;
	std::vector<double> pvin, pvout, ab_range;
	PM(): XAS(false), RIXS(false), PE(false), eloss(true) {
		pvin = vecd(3,0);
		pvout = vecd(3,0);
		ab_range = vecd({-20,20});
	};
	void set_eloss(bool el) {
		if (RIXS) eloss = el;
		else if (el) std::cout << "eloss will not be set if RIXS is turned off" << std::endl;
	}
};

struct blapIndex {
	size_t g;
	size_t e;
	dcomp blap;
	blapIndex(size_t _g, size_t _e, dcomp _blap): g(_g), e(_e), blap(_blap) {};
	blapIndex(): g(0), e(0), blap(dcomp(0,0)) {};
};

std::string pol_str(const vecd& pvec);
std::complex<double> proj_pvec(int ml, const vecd& pvec);
vecd occupation(Hilbert& hilbs, const std::vector<bindex>& si, bool is_print=true, 
				std::string fname = "");
vecd wvfnc_weight(Hilbert& hilbs, const std::vector<bindex>& si, 
								int ligNum = 3, bool print = false);
double effective_delta(Hilbert& hilbs, int ligNum = 3, bool is_print = false);
void state_composition(Hilbert& hilbs, const std::vector<bindex>& si, size_t top = 10);
void basis_overlap(Hilbert& GS, Hilbert& EX, bindex inds, std::vector<blapIndex>& blap, 
					const PM& pm, bool pvout = false);
// XAS Functions
void XAS_peak_occupation(Hilbert& GS, Hilbert& EX, vecd const& peak_en, vecd const& energy, 
						vecd const& intensity, std::vector<bindex> const& gsi, double ref_en = 0, 
						std::string mode = "list", bool ref_gs = false);
void write_XAS(vecd const& aben, vecd const& intensity, std::string file_dir = "");
void XAS(Hilbert& GS, Hilbert& EX, const PM& pm);

// RIXS functions
void RIXS_peak_occupation(Hilbert& GS, Hilbert& EX, vecd const& peak_en, vecd const& ab_en, 
						vecd const& em_en, vecd const& intensity, std::vector<bindex> const& gsi,
						PM const& pm, double ref_en = 0, std::string mode = "list", bool ref_gs = false);
void write_RIXS(vecd const& peaks, vecd const& ab, vecd const& em, bool eloss, 
				std::string file_dir = "");
void RIXS(Hilbert& GS, Hilbert& EX, const PM& pm);

#endif