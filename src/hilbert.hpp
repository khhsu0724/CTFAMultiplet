#include "diagonalize.hpp"
#include "cluster.hpp"
#ifndef HILBERT
#define HILBERT

// template <typename T>
class Hilbert {
public:
	// TODO: make some of the variables private
	size_t hsize = 1;
	int num_at = 0, at_per_site = 0;
	int val_ati = 0, val_ind = 0; // Index/Orbital index of first valence atom
	int num_ch = 0, num_vh = 0, num_corb = 0, num_vorb = 0; // This is number of orbital*2 (number of electron sites)
	bool SO_on = false, CV_on = false, CF_on = false, HYB_on = false;
	bool is_ex, BLOCK_DIAG = false;
	std::string coord = "none", edge, inp_hyb_file;
	std::vector<Atom> atlist; // atlist is ordered
	std::vector<Block<double>> hblks;
	std::vector<int> sites = {1,1,1};
	std::vector<int> orb_atom_ind; // return the index of atom for the orbital
	using Hashptr = bindex (Hilbert::*)(ulli s);
	using HBptr = ulli (Hilbert::*)(bindex ind);
	Hashptr hashfunc;
	HBptr hbfunc;
	Cluster* cluster = NULL;

public:
	Hilbert() {};
	Hilbert(const Hilbert &hilbs, int vh_mod = 0);
	~Hilbert() {};
	explicit Hilbert(std::string file_dir, const HParam& hparam, // HParam referenced from atoms.hpp
						std::string edge, bool is_ex = false);
	std::vector<ulli> enum_hspace(ulli inc_val = 0, ulli inc_core = 0, int vmod = 0, int cmod = 0);
	ulli qn2ulli(int snum, QN* qn, bool only_val = false, bool only_core = false);
	vpulli match(int snum, QN* lhs, QN* rhs);
	void fill_hblk(double const& matelem, ulli const& lhs, ulli const& rhs);
	void print_bits(ulli state);
	double Fsign(QN* op, ulli state, int opnum);
	double Fsign(ulli* op, ulli state, int opnum);
	int orbind(ulli s);
	int tot_site_num();
	double pheshift(double trace, int k);
	std::vector<double> get_all_eigval(bool is_err = true);
	std::vector<ulli> get_hashback_list(size_t blk_ind);

	// Functions for input file parsing and initialize
	void assign_cluster(std::string input);
	void read_from_file(std::string file_dir);

	// Nice Collection of Hash Functions
	bindex Hash(ulli s) {return (this->*hashfunc)(s);};
	ulli Hashback(bindex ind) {return (this->*hbfunc)(ind);};
	bindex norm_Hash(ulli s);
	ulli norm_Hashback(bindex ind);
	bindex sz_Hash(ulli s);
	ulli sz_Hashback(bindex ind);
	bindex jz_Hash(ulli s);
	ulli jz_Hashback(bindex ind);
	bindex ksz_Hash(ulli s);
	ulli ksz_Hashback(bindex ind);
	bindex k_Hash(ulli s);
	ulli k_Hashback(bindex ind);

	// Need to fix these
	// Maybe a copy constructor?????
	// std::vector<double> momentum(double* eigvec, bool return_square = true);
	// void momentum_check(double* mat, double* eig, double* eigvec);
};

#endif