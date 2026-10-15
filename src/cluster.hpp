#include "atoms.hpp"
#ifndef CLUSTER
#define CLUSTER

// Includes information about the cluster (Geometry, Hybridization)
class Cluster {
public:
	int vo_persite = 0, co_persite = 0;
	int lig_per_site = 0, tm_per_site = 0;
	bool no_HYB;
	std::string edge, inp_hyb_file = "";
	virtual void set_hyb_params(const HParam& hparam) {return;};
	// Converts spherical harmonics to real space basis, per site
	virtual vecc get_seph2real_mat() = 0;
	// Redefine operators to bring hybridization matrix into real matrix
	virtual vecc get_operator_mat() {return vecc(0);};
	// t is the real-time propagation time (natural units, eV^-1); default 0 recovers
	// the equilibrium (no pump) hybridization matrix. Only the base class declares
	// the default argument -- every call site in this code goes through a Cluster*
	// (see Hilbert::cluster), so there is no risk of default-argument divergence.
	virtual vecc get_tmat(double t = 0.0) {return vecc(0);};
// Shared Implementation
public:
	Cluster(std::string _edge): edge(_edge) {};
	~Cluster() {};
	int at_per_site() {return lig_per_site + tm_per_site;};
	const std::vector<std::string>& get_orb_names() {return orb_names;};
	void set_print_site(bool print_site_occ) {this->print_all_sites = print_site_occ;};
	void set_num_sites(int _num_sites) {this->num_sites = _num_sites;};
	void read_inp_tmat(std::string inp_hyb_file) {this->inp_hyb_file = inp_hyb_file;};
	void make_atlist(std::vector<Atom>& atlist, int num_vh,
						const std::vector<int>& sites);
	void print_eigstate(const vecd& occ, bool is_print = true, std::string fname = "", int p = 5);
	vecc get_inp_tmat();
	vecd get_tmat_real(double t = 0.0);
protected:
	int w = 12, num_sites = 1;
	std::vector<std::string> orb_names;
	bool print_all_sites = true; // If false program will print per_site
	bool tppsigma_on = false;
	bool check_tmat_all_real(const vecc& tmatreal);
	// Pump (Peierls shift) parameters, shared by every geometry. Captured from HParam
	// by capture_pump_params(), called from each subclass's set_hyb_params(). All zero
	// by default, so get_tmat(t) is time-independent unless a pump is configured.
	double pump_freq = 0, pump_sigma = 0, pump_t0 = 0;
	double pump_A0[3] = {0,0,0};
	void capture_pump_params(const HParam& hparam) {
		pump_freq = hparam.pump_freq;
		pump_sigma = hparam.pump_sigma;
		pump_t0 = hparam.pump_t0;
		pump_A0[0] = hparam.pump_amp_x;
		pump_A0[1] = hparam.pump_amp_y;
		pump_A0[2] = hparam.pump_amp_z;
		if (!inp_hyb_file.empty() && (pump_A0[0] || pump_A0[1] || pump_A0[2]))
			throw std::invalid_argument("Pump field is not supported with an input "
				"hybridization matrix (HYBMAT) -- it has no momentum/bond-direction concept");
	};
	double pump_shift(double t, int dir) const {
		// A(t) = A0 * sin(pump_freq*(t-t0)) * exp(-(t-t0)^2/(2*pump_sigma^2)) -- matches
		// the Devereaux-group reference convention (sin, not cos) exactly. pump_freq is
		// the angular frequency directly (e.g. pi for their omega_pump=pi eV, NOT a
		// cycles-per-time value needing a further 2*pi factor).
		if (pump_sigma <= 0) return 0.0;
		double dtc = t - pump_t0;
		return pump_A0[dir] * std::exp(-dtc*dtc/(2*pump_sigma*pump_sigma)) * std::sin(pump_freq*dtc);
	};
	vecc U_d() {
		vecc U(5*5,0);
		U[0*5+0] = U[0*5+4] = {1/sqrt(2),0};
		U[1*5+2] = {1,0};
		U[2*5+0] = {0,1/sqrt(2)};
		U[2*5+4] = {0,-1/sqrt(2)};
		U[3*5+1] = {1/sqrt(2),0};
		U[3*5+3] = {-1/sqrt(2),0};
		U[4*5+1] = U[4*5+3] = {0,1/sqrt(2)};
		return U;
	};
	vecc U_p() {
		vecc U(3*3,0);
		U[0*3+0] = U[0*3+2] = {1/sqrt(2),0};
		U[1*3+0] = {0,1/sqrt(2)};
		U[1*3+2] = {0,-1/sqrt(2)};
		U[2*3+1] = {1,0};
		return U;
	};
	vecc U_s() {
		return vecc({std::complex<double>(1,0)});
	};
	// Generate site number in given order. i.e. input 0, (0,0,0), 1 (0,0,1)
	// Maybe do this in hilbert.cpp?
};

class Ion: public Cluster {
public:
	Ion(std::string edge);
	vecc get_seph2real_mat() {return U_d();};
};

class SquarePlanar : public Cluster {
public:
	double tpd = 0, tpp = 0, del = 0;
	double tpp_pi = 0; // ligand-ligand pi-bonding, see get_tmat() for the sign structure
	SquarePlanar(std::string edge);
	void set_hyb_params(const HParam& hparam);
	vecc get_seph2real_mat();
	vecc get_operator_mat();
	vecc get_tmat(double t); // default lives on the Cluster base virtual only
private:
	double sig_pi = 0.5, tpdz_ratio = 0;
};

class Octahedral : public Cluster {
public:
	double tpd = 0, tpp = 0, del = 0, octJT = 1;
	Octahedral(std::string edge);
	void set_hyb_params(const HParam& hparam);
	vecc get_seph2real_mat();
	vecc get_operator_mat();
	vecc get_tmat(double t); // default lives on the Cluster base virtual only
private:
	double tpd_sig_pi = 0.5; // This needed to be checked
	double tpp_sig_pi = 0;
};

#endif