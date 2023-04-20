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
	virtual vecc get_tmat() {return vecc(0);};
// Shared Implementation
public:
	Cluster(std::string _edge): edge(_edge) {};
	~Cluster() {};
	int at_per_site() {return lig_per_site + tm_per_site;};
	void set_print_site(bool print_site_occ) {this->print_all_sites = print_site_occ;};
	void set_num_sites(int _num_sites) {this->num_sites = _num_sites;};
	void read_inp_tmat(std::string inp_hyb_file) {this->inp_hyb_file = inp_hyb_file;};
	void make_atlist(std::vector<Atom>& atlist, int num_vh, 
						const std::vector<int>& sites);
	void print_eigstate(const vecd& occ, std::string fname = "", int p = 5);
	vecc get_inp_tmat();
	vecd get_tmat_real();
protected:
	int w = 12, num_sites = 1;
	std::vector<std::string> orb_names;
	bool print_all_sites = true; // If false program will print per_site
	bool check_tmat_all_real(const vecc& tmatreal);
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
	SquarePlanar(std::string edge);
	void set_hyb_params(const HParam& hparam);
	vecc get_seph2real_mat();
	vecc get_operator_mat();
	vecc get_tmat();
private:
	double sig_pi = 0.5;
};

class Octahedral : public Cluster {
public:
	double tpd = 0, tpp = 0, del = 0;
	Octahedral(std::string edge);
	void set_hyb_params(const HParam& hparam);
	vecc get_seph2real_mat();
	vecc get_operator_mat();
	vecc get_tmat();
private:
	double tpd_sig_pi = 0.4; // This needed to be checked
	double tpp_sig_pi_1 = 0;
	double tpp_sig_pi_2 = 0;
};

#endif