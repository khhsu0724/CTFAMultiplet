#include "helper.hpp"
#ifndef ATOMS
#define ATOMS

inline int conv_lchar(char orb) {
	try {
		if (orb == 's') return 0;
		else if (orb == 'p') return 1;
		else if (orb == 'd') return 2;
		else if (orb == 'f') return 3;
		else throw std::invalid_argument("invalid orbital");
	} catch(const std::exception &ex) {
		std::cerr << ex.what() << "\n";
		exit(0);
	}
};

inline bool ao_order(int n1, int l1, int n2, int l2) {
	// Check if n1,l1 is a higher/equal orbital than n2,l2
	if (n1+l1 > n2+l2) return true;
	else if (n1+l1 < n2+l2) return false;
	else return (n1 >= n2);
};

struct QN {
	int ml = 0, order = 0;
	double spin = 0;
	QN(int ml, double spin, int order = 0): 
		ml(ml), spin(spin), order(order) {};
	QN() {};
};

inline bool operator<(const QN& qn1, const QN& qn2) {
	// State further to the "right" is smaller
	if (qn1.order != qn2.order) return qn1.order < qn2.order;
	if (qn1.spin != qn2.spin) return qn1.spin < qn2.spin;
	else return qn1.ml > qn2.ml;
};

inline bool operator!=(const QN& qn1, const QN& qn2) {
	if (qn1.spin == qn2.spin && qn1.ml == qn2.ml && qn1.order == qn2.order) return false;
	return true;
};

struct HParam {
	// Hilbert space parameters
	double nedos = 0; // PHASE OUT
	double HFscale = 1.0;
	double MLdelta = 0, octJT = 1, sig_pi = 0.5;
	double tpd = 0, tpp = 0, tpdz_ratio = 0.25;
	bool tppsigma_on = false;
	double SO[3]{0}, CF[5]{0};
	double SC2[5]{0}, SC1[3]{0}, FG[4]{0}, SC2EX[5]{0};
	int gs_diag_option = 2, ex_diag_option = 2;
	bool block_diag = true, HYB = true, effective_delta = true;
	bool print_site_occ = false;
	int ex_nev = 0, gs_nev = 0;
	std::vector<double*> SC;
	HParam() {
		SO[2] = -1; // Set this to -1 for checking
		SC.emplace_back(SC1);
		SC.emplace_back(SC2);
	};
	// Copy constructor??
};

struct Atom {
	// TODO: make some of the variables private
	int val_n,val_l,n,l,num_h,atnum;
	int sind,eind,atind; // atom variable tracks which atom it is
	int vind = -1, cind = -1;
	ulli check = 0; 	   	   // Convenient Number to check if a digit contains in this atom
	bool is_val, is_lig;
	std::string atname;
	std::vector<int> site;
	Atom() {};
	Atom(std::string orb, int atind, int val_n, int val_l, int num_h, std::vector<int> site):
		atind(atind), site(site), atname("null") {
		if (orb.size() != 2) std::cerr << "invalid orbital input\n";
		n = (int) orb[0] - '0';
		l = conv_lchar(orb[1]);
		set_num_h(val_n,val_l,num_h,0);
	};
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
		int i = __builtin_ctzll(s); // Should be the same as log2
		if (i < half_orb) qn = QN(l-i+sind,-0.5,order);
		else qn = QN(l-i+half_orb+sind,0.5,order);
		return qn;
	};
	QN get_qn(ulli s, int half_orb, int order = 0) {
		QN qn(0,0,order);
		for (int i = sind; i <= eind; ++i) {
			if (s & BIG1<<i) {
				qn.ml += l - i + sind;
				qn.spin -= 0.5;
			}
			if (s & BIG1<<(i+half_orb)) {
				qn.ml += l - i + sind;
				qn.spin += 0.5;
			}
		}
		return qn;
	};
	int get_occ(ulli s, int half_orb) {return ed::count_bits(s & check);};
	bool contains(ulli const& s) {return s & check;};
};

#endif