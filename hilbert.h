#include <iostream>
#include <string>
#include <utility>
#include <tuple>
#ifndef HILBERT
#define HILBERT

int choose(int n, int k);
double norm(std::vector<double>& vin);

struct QN {
	// Quantum number for each occupations
	int ml;
	int spin;
	QN(int ml, int spin) {
		this->ml = ml;
		this->spin = spin;
	}
	QN() {}
};
bool operator<(const QN& qn1, const QN& qn2);
bool operator!=(const QN& qn1, const QN& qn2);

class Hilbert {
	// This class contains function to calculate and generate hilbert space (array of states)
	// converts bit representation of occupied state to quantum numbers
private:
	int n;
	char orb;
public:
	int get_n() const {return n;};
	char get_orb() const {return orb;};
	int l;
	int occ_num;
	int orb_avail;
	int hmat_size; // Number of combinations of this orbital
	int t_hsize; // the total hilbert space size
	int* hmat;
	int hindex; // starting index of the orbital in hilbert space
	std::vector<std::vector<int>> site_ind; // index for Hilbert Space when there are multiple sites/orbitals
public:
	Hilbert() {};
	~Hilbert() {};
	explicit Hilbert(int n, char orb, int occ_num);
	double Psign(QN* lhsop, QN* rhsop, int lhss, int rhss, int lopnum, int ropnum);
	int bit2state(std::string& state);
	std::string state2bit(int state);
	int qn2state(QN* qn, int snum, bool is_bit=true);
	int sindex(int state);
	QN index2qn(int index);
	void fill_in_mat(double* mat, double matelem, int lhs, int rhs);
	int* generate_states(int in_state, QN* in_state_arr);
	std::vector<std::pair<int,int>> match_states(int snum, QN* rhs, QN* lhs);
	void pretty_print(double* mat, std::pair<int,int> column, std::pair<int,int> row);
	std::vector<double> momentum(double* eigvec, bool return_square = true);
	void momentum_check(double* mat, double* eig, double* eigvec);
	double pheshift(double trace, int k);
};

void cv_interaction(); // Core Valence Interaction

#endif