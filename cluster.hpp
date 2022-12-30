#include "atoms.hpp"
#ifndef CLUSTER
#define CLUSTER

// Includes information about the cluster (Geometry, Hybridization)
class Cluster {
public:
	int vo_persite, co_persite;
	bool no_HYB;
	std::string edge;
	virtual void set_hyb_params(const HParam& hparam) {return;}
	// Converts spherical harmonics to real space basis
	virtual vecc get_seph2real_mat() = 0;
	// Redefine operators to bring hybridization matrix into real matrix
	virtual vecc get_operator_mat() {return vecc(0);};
	virtual vecc get_tmat() {return vecc(0);};
	virtual void print_eigstate(const vecd& occ, int p = 5) = 0;

// Shared Implementation
public: 
	Cluster(std::string _edge): edge(_edge) {};
	~Cluster() {};
	void make_atlist(std::vector<Atom>& atlist, int num_vh, 
						const std::vector<int>& sites, int num_sites) {
		std::vector<int> nh_per_tm = ed::distribute(num_vh,tm_per_site);
		int atind = 0;
		// Push Core Atoms
		atlist.reserve(at_per_site()*num_sites);
		for (int x = 0; x < sites[0]; ++x) {
		for (int y = 0; y < sites[1]; ++y) {
		for (int z = 0; z < sites[2]; ++z) {
			std::vector<int> site = {x,y,z};
			for (size_t tm = 0; tm < tm_per_site; ++tm) {
				if (edge == "L") atlist.emplace(atlist.begin(),Atom("2p",atind,3,2,0,site));
				atlist.emplace_back(Atom("3d",atind,3,2,nh_per_tm[tm],site));
				atind++;
			}
			for (size_t lig = 0; lig < lig_per_site; ++lig) {
				if (edge == "K") atlist.emplace(atlist.begin(),Atom("1s",atind,2,1,0,site));
				atlist.emplace_back(Atom("2p",atind,2,1,0,site));
				atind++;
			}
		}}}
	};
	int at_per_site() {return lig_per_site + tm_per_site;};
	vecd get_tmat_real();
protected:
	int lig_per_site = 0, tm_per_site = 0, w = 12;
	void print_state_header(const vecd& occ);
	void print_state_orbs(const vecd& occ, const std::vector<std::string>& names, int p);
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
	void print_eigstate(const vecd& occ, int p = 5);
};

class SquarePlanar : public Cluster {
public:
	double tpd = 0, tpp = 0, del = 0;
	SquarePlanar(std::string edge);
	void set_hyb_params(const HParam& hparam);
	vecc get_seph2real_mat();
	vecc get_operator_mat();
	vecc get_tmat();
	void print_eigstate(const vecd& occ, int p = 5);
private:
	double sig_pi = -0.5;
};

#endif