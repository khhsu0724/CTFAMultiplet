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
	bool cross = false;
	bool skip_ch_diag = false;
	int nedos = 1000, niterCFE = 150;
	double CG_tol = 1e-8;
	int spec_solver = 1; // 1 = exact, 2 = Classic K-H, 3 = both, 4 Lanczos/BiCGS
	int precond = 0; // 0 = no preconditioner, 1 = supply initial guess from last incident e
	double em_energy = 15;
	double gamma = 0.3;
	double eps_ab = 0.1, eps_loss = 0.1;
	double abmax = -1;
	std::string edge;
	std::vector<double> pvin, pvout, ab_range;
	std::vector<double> incident, inc_e_points;
	PM(): XAS(false), RIXS(false), PE(false), eloss(true) {
		pvin = vecd(3,0);
		pvout = vecd(3,0);
		ab_range = vecd({-20,20});
		incident = vecd({0,0,0});
	};
	void set_eloss(bool el) {
		if (RIXS) eloss = el;
		else if (el) std::cout << "eloss will not be set if RIXS is turned off" << std::endl;
		return;
	}
	void set_incident_points(bool reset = false, double gs_en = 0, double ex_en = 0) {
		if (incident.size() != 3) {
			std::cout << "Incident energy not set correctly!" << std::endl;
			exit(1);
		}
		double start = incident.at(0);
		double end = incident.at(1);
		int num = incident.at(2);
		if (num < -1) {
			std::cout << "Incident[3] needs to >= -1" << std::endl;
			exit(1);
		}
		if (num == -1) {
			if (reset) {
				std::cout << "New range: ";
				vecd incident_copy(incident);
				incident[0] = ex_en - gs_en - 1;
				incident[1] = ex_en - gs_en + incident_copy[0] - 1;
				incident[2] = incident_copy[1];
				std::cout << incident[0] << ", " << incident[1] << ", " << incident[2] << std::endl;
				set_incident_points();
			} else {
				std::cout << "Warning: auto detect range" << std::endl;
				std::cout << "Incident[1] = Energy range starting from lowest excitation " << std::endl;
				std::cout << "Incident[2] = Number of points " << std::endl;
				if (incident.at(0) <= 0 or incident.at(1) < 0) {
					std::cout << "Range/Points needs to > 0" << std::endl;
					exit(1);
				}
			}
		}
		if (num == 0) return;
	    if (num == 1) {
	        inc_e_points.push_back(start);
	        return;
	    }
	    double step = (end - start) / (1.0*(num - 1));
	    for (int i = 0; i < num; ++i) {
	        inc_e_points.push_back(start + i * step);
	    }
	    return;
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
				std::string fname = "", bool spin_res = false);
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
void write_XAS(vecd const& aben, vecd const& intensity, std::string file_dir = "", bool exact=true);
void XAS(Hilbert& GS, Hilbert& EX, const PM& pm);

// RIXS functions
void RIXS_peak_occupation(Hilbert& GS, Hilbert& EX, vecd const& peak_en, vecd const& ab_en, 
						vecd const& em_en, vecd const& intensity, std::vector<bindex> const& gsi,
						PM const& pm, double ref_en = 0, std::string mode = "list", bool ref_gs = false);
void write_RIXS(vecd const& peaks, vecd const& ab, vecd const& em, bool eloss, 
				std::string file_dir = "");
void RIXS(Hilbert& GS, Hilbert& EX, const PM& pm);
vecc gen_dipole_state(Hilbert& GS, Hilbert& EX, const PM& pm, const bindex& inds, vecc vec_in, 
						const std::vector<blapIndex>& blap, bool excite = true);

#endif