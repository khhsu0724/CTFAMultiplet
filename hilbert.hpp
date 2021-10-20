#include <unordered_map>
#include "helper.hpp"
#include "diagonalize.hpp"
#ifndef HILBERT
#define HILBERT

typedef unsigned long long int ulli;
typedef unsigned long int uli;
typedef std::vector<std::pair<ulli,ulli>> vpulli;

int conv_lchar(char orb);
bool ao_order(int n1, int l1, int n2, int l2);

struct QN {
	int ml = 0, order = 0;
	double spin = 0;
	QN(int ml, double spin, int order = 0): 
		ml(ml), spin(spin), order(order) {};
	QN() {};
};
bool operator<(const QN& qn1, const QN& qn2);
bool operator!=(const QN& qn1, const QN& qn2);

struct Atom {
	int val_n,val_l,n,l,num_h;
	int sind,eind,atnum,atind; // atom variable tracks which atom it is
	int vind = -1, cind = -1;
	ulli check = 0; 	   	   // Convenient Number to check if a digit contains in this atom
	bool is_val, is_lig;
	std::string atname;
	std::vector<double> pos, kp;
	Atom() {};
	void add_orb(std::string orb, int atind) {
		this->atind = atind;
		if (orb.size() != 2) std::cerr << "invalid orbital input\n";
		n = (int) orb[0] - '0';
		l = conv_lchar(orb[1]);
		if (atname == "Zn") set_num_h(3,2,0,30);
		else if (atname == "Cu") set_num_h(3,2,1,29);
		else if (atname == "Ni") set_num_h(3,2,2,28);
		else if (atname == "Co") set_num_h(3,2,3,27);
		else if (atname == "Fe") set_num_h(3,2,4,26);
		else if (atname == "Mn") set_num_h(3,2,5,25);
		else if (atname == "Cr") set_num_h(3,2,6,24);
		else if (atname == "V")  set_num_h(3,2,7,23);
		else if (atname == "Ti") set_num_h(3,2,8,22);
		else if (atname == "Sc") set_num_h(3,2,9,21);
		else if (atname == "O")  set_num_h(2,1,0,8);
		else {
			std::cerr << "Invalid atom input from INPUT file\n";
			exit(0);
		}
	};
	void set_num_h(int val_n_in, int val_l_in, int num_h_in, int atnum_in) {
		atnum = atnum_in, val_n = val_n_in, val_l = val_l_in;
		is_val = (n == val_n && l == val_l);
		is_lig = (val_n_in != 3); 				// If n != 3 then its not ligand?
		if (is_val) num_h = num_h_in;
		else num_h = 0;
	};
	QN fast_qn(ulli s, int half_orb, int order = 0) {
		// Fast qn if 's' represents one hole
		QN qn;
		int i = (int)log2(s);
		if (i < half_orb) qn = QN(l-i+sind,-0.5,order);
		else qn = QN(l-i+half_orb+sind,0.5,order);
		return qn;
	}
	QN get_qn(ulli s, int half_orb, int order = 0) {
		QN qn(0,0,order);
		for (int i = sind; i <= eind; ++i) {
			if (s & 1<<i) {
				qn.ml += l - i + sind;
				qn.spin -= 0.5;
			}
			if (s & 1<<(i+half_orb)) {
				qn.ml += l - i + sind;
				qn.spin += 0.5;
			}
		}
		return qn;
	};
	double get_k(ulli s, int half_orb) {
		double cnt = 0;
		for (int i = sind; i <= eind; ++i) {
			if (s & 1<<i) cnt++;
			if (s & 1<<(i+half_orb)) cnt++;
		}
		return ed::norm(kp)*cnt;
	};
	int get_occ(ulli s, int half_orb) {return ed::count_bits(s & check);};
	bool contains(ulli const& s) {return s & check;};
};

class Hilbert {
public:
	int hsize = 1, num_sites = 0;
	int val_ati = 0, val_ind = 0, num_vh = 0, num_vorb = 0; // These doess not account for up/down spin
	int num_ch = 0, num_corb = 0;
	bool SO = true, CF = true, CV = true, HYB = false;
	bool is_ex;
	std::vector<Atom> atlist; // atlist is ordered 
	std::vector<double> a1, a2, a3, b1, b2, b3;
	std::vector<int> rehash, states; // Phase out
	std::vector<Block> hblks;

public:
	Hilbert() {};
	~Hilbert() {};
	explicit Hilbert(std::string file_dir, bool is_ex = false);
	std::vector<ulli> enum_hspace(int inc_val = 0, int inc_core = 0, int vmod = 0, int cmod = 0);
	ulli qn2ulli(int snum, QN* qn, bool only_val = false, bool only_core = false);
	vpulli match(int snum, QN* lhs, QN* rhs);
	int Hash(ulli s);
	ulli Hashback(int ind);
	void fill_hblk(double const& matelem, ulli const& lhs, ulli const& rhs);
	double Fsign(QN* op, ulli state, int opnum);
	std::vector<Block> make_block(std::vector<ulli>& hilb_vec);
	int orbind(ulli s);
	void read_from_file(std::string file_dir);
	double pheshift(double trace, int k);

	// Need to fix these
	// Maybe a copy constructor?????
	// std::vector<double> momentum(double* eigvec, bool return_square = true);
	// void momentum_check(double* mat, double* eig, double* eigvec);
};

#endif