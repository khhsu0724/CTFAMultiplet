#include "cluster.hpp"
#include <tuple>
#ifndef PHONON
#define PHONON

typedef std::vector<std::pair<int,int>> vpi;
/**
 * contains implementation of phonon mode in small clusters
 * Uses: Cluster
 * The total Hilbert space is a product of electronic + phonon part
*/

// Holds pair of bindex and prefactor
struct pn_pair {
	bindex lhs, rhs;
	double prefac;
	pn_pair() {};
	pn_pair(double prefac, bindex lhs, bindex rhs): 
			prefac(prefac), lhs(lhs), rhs(rhs) {};
	~pn_pair() = default;
};

// phonon operator
struct pn_op {
	bool is_raise; // Raising vs Lowering
	int mode;
	pn_op() {};
	pn_op(int mode, bool is_raise): mode(mode), is_raise(is_raise) {};
	~pn_op() = default;
};

// A structure that holds terms and parameters for phonon Hamiltonian
struct Phonon_Param {
	// Bare phonon parameters: wb^b
	double omega = 0;

	// Holestein parameters: gn(b^+b)
	veci h_orb_i;
	double g_h = 0;

	// Breathing parameters: gc^icj(b^+b)
	vpi b_orb_ij;
	double g_b = 0;

	// Constructors
	Phonon_Param() {};
	Phonon_Param(const Phonon_Param& pnParam):
		omega(pnParam.omega), h_orb_i(pnParam.h_orb_i), g_h(pnParam.g_h),
		b_orb_ij(pnParam.b_orb_ij), g_b(pnParam.g_b) {};
	~Phonon_Param() = default;
};

// Contains information for max number of phonons
class Phonon {
private:
public:
	int mode_id, max_phonons=0;
	Phonon_Param pnParam;
	// Collection of properties and function that returns member of this class
	Phonon(int mode_id, int max_phonons, const Phonon_Param& pnParam): 
		mode_id(mode_id), max_phonons(max_phonons), pnParam(pnParam) {
			// std::cout << "mode: " << mode_id << ", max_phonons: " << max_phonons << ", max_phonons: " << max_phonons << std::endl;
			// std::cout << "bit start: " << bit_start << ", bit end: " << bit_end << std::endl;
	};
public:
	bool is_holestein() {return pnParam.h_orb_i.size();};
	bool is_breathing() {return pnParam.b_orb_ij.size();};
	veci return_state(veci state_in, bool is_raise, bool& valid_op, double& prefac) {
		// Return the phonon state and prefactor,
		// given a series of raising/lowering operators
		if (is_raise) {
			if (state_in[mode_id] < max_phonons) {
				state_in[mode_id]++;
				prefac *= sqrt(state_in[mode_id]);
			} else valid_op = false;
		} else {
			if (state_in[mode_id] > 0) {
				prefac *= sqrt(state_in[mode_id]);
				state_in[mode_id]--;
			} else valid_op = false;
		}
		return state_in;
	};
	int return_num_phonons(const veci &state) {
		return state[mode_id];
	};
};


// Holestein term, modulate density terms
// 			 +
// H = gn  (b + b )
//       is  i   i
// class Density: public Phonons {
// public:
// 	double omega;
// }

// Breathing mode, modulates hopping terms
// 			  +    +
// H = gQP  (c c)(b + b )
//        is  j k  i   i
// class Breathing: public Phonons {
// public:
// 	double omega
	
// }


// Hosts description of phonons
// Contains a list of phonon that describes the interaction and frequency/parameter for each phonon
class PhononModes {
private:
    veci strides;
    int max_tot_phonons = 0, pn_hsize = 1;
    int w = PRINT_WIDTH;
	Cluster* cluster = NULL;
	bool check_orb_range(veci h_orb) {
        // Check for holestein mode, Is it supposed to act on core orbitals?
        for (auto & orb : h_orb) {
			if (orb < 0 || orb > (cluster->vo_persite))
				return false;
        }
		return true;
	};

	bool check_orb_range(vpi b_orb) {
        // Is it supposed to act on core orbitals?
        for (auto & orbs : b_orb) {
			if (orbs.first < 0 || orbs.first > (cluster->vo_persite))
				return false;
			if (orbs.second < 0 || orbs.second > (cluster->vo_persite))
				return false;
    	}
		return true;
	};
public:
    std::vector<std::shared_ptr<Phonon>> phonons;
	using Hashptr = bindex (PhononModes::*)(const veci &s);
	using HBptr = veci (PhononModes::*)(bindex ind);
	Hashptr hashfunc;
	HBptr hbfunc;
	PhononModes(Cluster* cluster) {
		this->cluster = cluster;
		this->cluster->print_cluster_name();
	};
    ~PhononModes() = default;
    void addPhonon(std::shared_ptr<Phonon> phonon) {
        this->phonons.push_back(phonon);
    };
public:
    void addPhonon(int max_phonons, const Phonon_Param& pnParam) {
    	try {
		    if (pnParam.h_orb_i.size() && !check_orb_range(pnParam.h_orb_i)) 
		    	throw std::out_of_range("Holestein");
		    if (pnParam.b_orb_ij.size() && !check_orb_range(pnParam.b_orb_ij)) 
		    	throw std::out_of_range("Breathing");
		} catch (const std::out_of_range& e) {
		    std::cerr << "Error: " << e.what() << " mode orbital index out of valid range."<< std::endl;
		    std::exit(EXIT_FAILURE);
		}
        this->phonons.push_back(std::make_shared<Phonon>(get_unique_modes(),max_phonons,pnParam));
        // Recount unique phonon modes and phonon per mode after this
		pn_hsize *= (max_phonons+1);
    	max_tot_phonons += max_phonons;
    	// std::cout << "Unique Modes: " << get_unique_modes() << ", pn hsize: " << pn_hsize << ", max phonon: " << max_tot_phonons << std::endl;
		// Calculate stride to help indexing
		strides.clear();
		for (int i = 0; i < phonons.size(); ++i) {
			int stride_cnt = 1;
			for (int j = phonons.size()-1; j > i; --j)
				stride_cnt *= (this->phonons.at(j)->max_phonons+1);
			strides.push_back(stride_cnt);
		}
    };
    int get_max_phonons() {return this->max_tot_phonons;};
    int get_unique_modes() {return phonons.size();};
    int get_pn_hsize() {return this->pn_hsize;};
    void check_valid_state(veci const & pnarr) {
    	try {
    		if (pnarr.size() != get_unique_modes())
    			throw std::invalid_argument("(wrong size)");
    		for (int i = 0; i < get_unique_modes(); ++i) {
    			if (pnarr[i] > this->phonons[i]->max_phonons)
    				throw std::invalid_argument("(over max phonons)");
    		}
    	} catch (const std::invalid_argument& e) {
    		std::cerr << "Error: Incorrect phonon state " << e.what() << std::endl;
    	}
    	return;
    }

    void print_phonon(bindex const & pnind) {
    	print_phonon(this->norm_Hashback(pnind));
    	return;
    }

    void print_phonon(veci const & pnarr) {
    	std::cout << "[";
    	for (auto & pn : pnarr) {
    		std::cout << pn << ",";
    	}
    	std::cout << "]";
    	return;
    }
    void check_pn_hspace() {
    	// Function that checks if hashing and hashback is correct
    	for (int p = 0; p < pn_hsize; ++p) {
    		int ind = p;
    		veci state = this->norm_Hashback(bindex(0,ind));
    		bindex index = this->norm_Hash(state);
    		std::cout << "index: " << ind << ": ";
    		this->print_phonon(state);
    		std::cout << ", hasback: " << index.second;
    		std::cout << std::endl;
    	}
    };

    std::vector<pn_pair> gen_pair_pn_index(std::vector<pn_op> phonon_ops) {
    	// Given a list operator, provide the resulting index for all states in Hilbert space
  		// Also returns prefactor coming out from the operator
    	std::vector<pn_pair> matches;
    	for (int p = 0; p < pn_hsize; ++p) {
    		veci init_state = this->norm_Hashback(bindex(0,p));
    		veci after_state = init_state;
    		// std::cout << "init: ";
    		// this->print_phonon(init_state);
    		bool valid_op = true;
			double prefac = 1;
    		for (auto & pnop : phonon_ops) {
    			if (!valid_op) break;
    			if (pnop.mode >= get_unique_modes()) {
    				std::cerr << "Error: invalid phonon operator" << std::endl;
    				std::cerr << "total modes: " << get_unique_modes();
    				std::cerr << ", operator mode: " << pnop.mode << std::endl;
    				exit(1); // Should've used try/catch
    			}
    			after_state = this->phonons[pnop.mode]->return_state(after_state,pnop.is_raise,valid_op,prefac);	
    		}
    		auto after_ind = this->norm_Hash(after_state);
    		// if (valid_op) {
    		// 	std::cout << ", prefac: " << prefac << ", after: ";
    		// 	this->print_phonon(init_state);
    		// }
			// else std::cout << ", INVALID";
    		// std::cout << std::endl;
			if (valid_op) matches.emplace_back(pn_pair(prefac,bindex(0,p),bindex(0,after_ind.second)));
    	}
    	return matches;
    };
	void print_phonons_occuptation(const vecd& occ, bool is_print, std::string fname);

    // Hash functions for phonon hilbert space
	bindex Hash(const veci &s) {return (this->*hashfunc)(s);};
	veci Hashback(bindex ind) {return (this->*hbfunc)(ind);};
	bindex norm_Hash(const veci &s);
	veci norm_Hashback(bindex ind);

};


#endif