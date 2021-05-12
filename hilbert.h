#include <iostream>
#include <string>
#include <utility>
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
	int l;
	int occ_num;
	int orb_avail;
	int hmat_size;
	int* hmat;
public:
	explicit Hilbert(char orb, int occ_num);
	double Psign(QN* lhsop, QN* rhsop, int lhss, int rhss, int opnum);
	int bit2state(std::string& state);
	std::string state2bit(int state);
	int qn2state(QN* qn, int snum, bool is_bit=true);
	int sindex(int state);
	int* generate_states(int in_state, QN* in_state_arr);
	std::vector<std::pair<int,int>> match_states(int snum, QN* rhs, QN* lhs);
	void pretty_print(double* mat, std::pair<int,int> column, std::pair<int,int> row);
	void momentum_check(double* mat, double* eig, double* eigvec);
};

#endif