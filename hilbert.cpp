#include <fstream>
#include <regex>
#include <algorithm>
#include <cmath>
#include <bitset>
#include "hilbert.hpp"

using namespace std;

int conv_lchar(char orb) {
	try {
		if (orb == 's') return 0;
		else if (orb == 'p') return 1;
		else if (orb == 'd') return 2;
		else if (orb == 'f') return 3;
		else throw invalid_argument("invalid orbital");
	} catch(const exception &ex) {
		cerr << ex.what() << "\n";
		exit(0);
	}
}

bool ao_order(int n1, int l1, int n2, int l2) {
	// Check if n1,l1 is a higher/equal orbital than n2,l2
	if (n1+l1 > n2+l2) return true;
	else if (n1+l1 < n2+l2) return false;
	else return (n1 >= n2);
}

bool operator<(const QN& qn1, const QN& qn2) {
	// State further to the "right" is smaller
	if (qn1.order != qn2.order) return qn1.order < qn2.order;
	if (qn1.spin != qn2.spin) return qn1.spin < qn2.spin;
	else return qn1.ml > qn2.ml;
}

bool operator!=(const QN& qn1, const QN& qn2) {
	if (qn1.spin == qn2.spin && qn1.ml == qn2.ml && qn1.order == qn2.order) return false;
	return true;
}

Hilbert::Hilbert(string file_dir, bool is_ex): is_ex(is_ex) {
	read_from_file(file_dir);
	if (is_ex) {
		SO = true;
		num_vh--;
		num_ch++;
	} else SO = false;

	// Fix this so the code doesnt have to rely on a vector allocated
	vector<ulli> hspace = enum_hspace();
	// for (auto & h : hspace) {
	// 	cout << h << ", " << bitset<16>(h) << ", Hash: " << Hash(h) << ", Hashback: " << Hashback(Hash(h)) << endl;
	// 	if (h != Hashback(Hash(h))) cout << "error hashing" << endl;
	// }
	// for (auto & at : atlist) cout << at.atname << ", " << at.is_lig << endl;
	hsize = ed::choose(num_vorb,num_vh)*ed::choose(num_corb,num_ch);
	make_block(hspace);
}

vector<ulli> Hilbert::enum_hspace(int inc_val, int inc_core, int vmod, int cmod) {
	vector<ulli> val,core,hspace;
	ed::enum_states(val,num_vorb,num_vh+vmod,inc_val);
	ed::enum_states(core,num_corb,num_ch+cmod,inc_core);
	for (auto &c : core) for (auto &v : val) hspace.push_back(ed::add_bits(v,c,num_vorb,num_corb));
	return hspace;
}


// Using 2^n representation to calculate quantumn number, convert to binary
// and +1/2 spin to -1/2. Angular momentum comes first, hole language
inline ulli Hilbert::qn2ulli(int snum, QN* qn, bool only_val, bool only_core) {
	ulli s = 0;
	int num_orb = (only_val) ? num_vorb : (only_core) ? num_corb : num_corb + num_vorb;
	int shift = (only_val) ? val_ind : 0;
	for (size_t i = 0; i < snum; ++i) {
		Atom *qnat = &atlist[qn[i].order];
		if (qnat->l - qn[i].ml < 0) cerr << "quantum number larger than site's ml\n";
		int ind = qnat->sind - shift + (qnat->l - qn[i].ml) + (qn[i].spin+0.5)*num_orb/2;
		s += 1 << ind;
	}
	return s;
}

vpulli Hilbert::match(int snum, QN* lhs, QN* rhs) {
	// Match input states within the Hilbert Space, atomic number needs to match
	ulli rhs_val = 0, lhs_val = 0, rhs_core = 0, lhs_core = 0;
	for (size_t i = 0; i < snum; ++i) {
		if (atlist[rhs[i].order].is_val) rhs_val += qn2ulli(1,&rhs[i],true,false);
		else rhs_core += qn2ulli(1,&rhs[i],false,true);
		if (atlist[lhs[i].order].is_val) lhs_val += qn2ulli(1,&lhs[i],true,false);
		else lhs_core += qn2ulli(1,&lhs[i],false,true);
	}
	ulli inc_val = rhs_val|lhs_val, inc_core = rhs_core|lhs_core;
	vector<ulli> ms = enum_hspace(inc_val,inc_core,ed::count_bits(inc_val)-ed::count_bits(rhs_val)
									,ed::count_bits(inc_core)-ed::count_bits(rhs_core));
	vpulli mp;
	ulli rhs_all = qn2ulli(snum,rhs), lhs_all = qn2ulli(snum,lhs);
	ulli inc_all = rhs_all | lhs_all;
	for (auto &s : ms) mp.emplace_back(pair<ulli,ulli>(s-inc_all+lhs_all,s-inc_all+rhs_all));
	return mp;
}

int Hilbert::Hash(ulli s) { // Add symmetries
	// Hash function that convert a state in bits to index
	int cind = 0, vind = 0;
	int cnt_c = 0, cnt_v = 0;
	int chash = num_corb/2, vhash = num_vorb/2;
	for (int i = 0; i < num_corb/2; ++i)
		if (s & (1 << i)) cind += ed::choose(i,++cnt_c);
	for (int i = 0; i < num_corb/2; ++i)
		if (s & (1 << (i+num_corb/2+num_vorb/2))) cind += ed::choose(i+chash,++cnt_c);
	for (int i = 0; i < num_vorb/2; ++i)
		if (s & (1 << (i+num_corb/2))) vind += ed::choose(i,++cnt_v);
	for (int i = 0; i < num_vorb/2; ++i)
		if (s & (1 << (i+num_corb+num_vorb/2))) vind += ed::choose(i+vhash,++cnt_v);
	// if ((cnt_cu+cnt_cd) != num_ch || (cnt_vu+cnt_vd) != num_vh) return -1;
	return vind + cind * ed::choose(num_vorb,num_vh);
}

ulli Hilbert::Hashback(int ind) { // Add symmetries
	// Hash function that convert index to state in bitset
	int cind = ind / ed::choose(num_vorb,num_vh), ch = num_ch;
	int vind = ind % ed::choose(num_vorb,num_vh), vh = num_vh;
	ulli c = 0, v = 0;
	for (int i = num_vorb-1; i >= 0; --i) {
		if (vind >= ed::choose(i,vh)) {
			v |= (1 << i);
			vind -= ed::choose(i,vh--);
		}
	}
	for (int i = num_corb-1; i >= 0; --i) {
		if (cind >= ed::choose(i,ch)) {
			c |= (1 << i);
			cind -= ed::choose(i,ch--);
		}
	}
	return ed::add_bits(v,c,num_vorb,num_corb);
}

double Hilbert::Fsign(QN* op, ulli state, int opnum) {
	// Calculate Fermion Sign Problem, return +1 or -1
	int p = 0;
	for (size_t i = 0; i < opnum; ++i) {
		ulli o = qn2ulli(1,op+i);
		if (!(state & o)) return 0; //Annihilate on vacuum
		state |= o;
		p += ed::count_bits(state/o);
	}
	return pow(-1,p);
}

vector<Block> Hilbert::make_block(vector<ulli>& hilb_vec) {
	// Allocate memory for the block matrices
	vector<Block> vb;
	double max_sz = (num_vh+num_ch <= num_vorb+num_corb) ? 
					(num_vh+num_ch)/2 : (num_vorb+num_corb-num_vh-num_ch)/2;
	// if (SO && !is_ex) {
	// 	// Group by K
	// } else if (CF || HYB) {
	// 	// Group by Sz and K

	// } else {
	// 	// Group by Sz K and L
	// }


	// !!!!Come up with a faster way to reserve memory for blocks!!!!
	hblks.push_back(Block(0,0,0,hsize));
	return vb;
}

void Hilbert::fill_hblk(double const& matelem, ulli const& lhs, ulli const& rhs) {
	// Find block index for lhs and rhs
	int lind = 0, rind = 0;
	// cout << "size: " << hblks[lind].size << ", Hash: " << Hash(lhs) << ", " << Hash(rhs);
	// cout << ", bits: " << bitset<16>(lhs) << ", " << bitset<16>(rhs) << endl;
	if (lind == rind) hblks[lind].ham[Hash(lhs)+hblks[lind].size*Hash(rhs)] += matelem;
	return;
}

int Hilbert::orbind(ulli s) {
	int i = 0;
	while (!(atlist[i].check & s)) ++i;
	return i;
}

void Hilbert::read_from_file(string file_dir) {
	string line;
	ifstream input(file_dir);
	vector<double> cell_param;
	try {
		if (input.is_open()) {
			bool read_cell = false, read_atoms = false;
			int atind = 0;
			while (getline (input,line)) {
				if (line[0] == '#') continue;
				if (line[0] == '/') {
					read_cell = false;
					read_atoms = false;
					continue;
				}
				if (read_cell) {
					string cp;
					for (int s = 0; s < line.size(); ++s) {
						if (line[s] != ' ') {
							cp.push_back(line[s]);
							if (s == line.size()-1) cell_param.emplace_back(stod(cp));
						}
						else {
							cell_param.emplace_back(stod(cp));
							cp = "";
						}
					}
				}
				if (read_atoms) {
					string at;
					Atom atin;
					int entry = 0;
					bool skipline = false;
					atind++;
					for (int s = 0; s < line.size(); ++s) {
						if (line[s] == '#') {
							atind--;
							skipline = true;
							break;
						}
						if (line[s] != ' ' && line[s] != '	') {
							at.push_back(line[s]);
							if (s == line.size()-1) {
								atin.add_orb(at,atind);
								atlist.emplace_back(atin);
								if (!atin.is_val) rotate(atlist.begin(),atlist.end()-1,atlist.end());
							}
						}
						else {
							if (at != "") {
								entry++;
								if (regex_match(at,regex(".*[a-zA-Z]+.*"))) {
									if (entry == 1) atin.atname = at;
									else {
										atin.add_orb(at,atind);
										atlist.emplace_back(atin);
										if (!atin.is_val) rotate(atlist.begin(),atlist.end()-1,atlist.end());
									}
								} else atin.pos.emplace_back(stod(at));
							}
							at = "";
						}
					}
					if (atin.pos.size() != 3 && !skipline) throw invalid_argument("Invalid atom position");
				}
				if (line == "&CELL") read_cell = true;
				else if (line == "&ATOM_POSITIONS") read_atoms = true;
			}
			input.close();
			int ind = 0;
			for (size_t i = 0; i < atlist.size(); ++i) {
				atlist[i].sind = ind;
				atlist[i].eind = ind + (2*atlist[i].l+1)-1;
				ind = atlist[i].eind + 1;
				if (atlist[i].is_val) {
					num_vh += atlist[i].num_h;
					num_vorb += 2*(2*atlist[i].l+1);
				} else {
					num_ch += atlist[i].num_h;
					num_corb += 2*(2*atlist[i].l+1);
					val_ind = atlist[i].eind + 1;
					val_ati = i + 1;
				}
				size_t j = 0;
				for (j = 0; j < i; ++j) if(ed::veccmp(atlist[i].pos,atlist[j].pos)) break;
				if (i == j) num_sites++;
			}
			if (cell_param.size() != 9)  throw invalid_argument("Invalid unit cell parameter");
			a1 = ed::slice(cell_param,0,2);
			a2 = ed::slice(cell_param,3,5);
			a3 = ed::slice(cell_param,6,8);
			double V = 	ed::dot(a1,ed::vec_cross(a2,a3));
			b1 = ed::vec_mult(ed::vec_cross(a2,a3),2*M_PI/V);
			b2 = ed::vec_mult(ed::vec_cross(a3,a1),2*M_PI/V);
			b3 = ed::vec_mult(ed::vec_cross(a1,a2),2*M_PI/V);
			for (auto &at : atlist) {
				for (size_t i = 0; i < 3; ++i) at.kp.emplace_back(at.pos[i]*(b1[i]+b2[i]+b3[i]));
				for (size_t i = at.sind; i <= at.eind; i++) at.check += (1 << i | 1 << (i+num_corb/2+num_vorb/2));
			}
			ind = 0;
			while(!atlist[ind].is_val) {
				int i = ind + 1;
				while (atlist[i].atind != atlist[ind].atind) ++i;
				atlist[ind].vind = i;
				atlist[i].cind = ind++;
			}
		}
		else throw invalid_argument("Cannot open INPUT file");
	} catch (const exception &ex) {
		cerr << ex.what() << "\n";
		exit(0);
	}
	return;
}

// double Hilbert::pheshift(double trace, int k) {
// 	// return value for particle-hole energy shift for Coulomb Interaction
// 	if (occ_num > orb_avail/2) {
// 		int p = ed::choose(occ_num,k);
// 		return trace/p*(p-ed::choose(orb_avail-occ_num,k))/t_hsize;
// 	} else return 0;
// }

// NEED TO FIX THESE

// vector<double> Hilbert::momentum(double* eigvec, bool return_square) {
// 	// Calculate L,S for specific eigenstate
// 	vector<double> Lz(hmat_size,0), Lp(hmat_size,0), Lm(hmat_size,0), Sz(hmat_size,0), Sp(hmat_size,0), Sm(hmat_size,0);
// 	double S2,L2,J2;
// 	if (occ_num == 0 || occ_num == orb_avail) return vector<double>{0.0,0.0,0.0};
// 	for (int j = 0; j < hmat_size; ++j) {
// 		if (abs(eigvec[j]) > 1e-7) {
// 			string state = state2bit(hmat[j]);
// 			for (int k = 0; k < orb_avail; ++k) {
// 				if (state[k] == '1') {
// 					double ml = (k % (orb_avail/2)) - l;
// 					Lz[j] += ml * eigvec[j];
// 					// Calculate L+
// 					string raise = state;
// 					if (ml < l && state[k+1] == '0') {
// 						raise[k] = '0';
// 						raise[k+1] = '1';
// 						Lp[sindex(bit2state(raise))] += sqrt((l-ml)*(l+ml+1)) * eigvec[j];
// 					}
// 					// Calculate L-
// 					string lower = state;
// 					if (ml > -l && state[k-1] == '0') {
// 						lower[k] = '0';
// 						lower[k-1] = '1';
// 						Lm[sindex(bit2state(lower))] += sqrt((l+ml)*(l-ml+1)) * eigvec[j];
// 					}
// 					// Calculate Spin
// 					if (k < (orb_avail/2)) {
// 						// Calculate Sz
// 						Sz[j] += 0.5 * eigvec[j];
// 						// Calculate S-
// 						int exchange = 0;
// 						string slower = state;
// 						if (state[k+orb_avail/2] == '0') {
// 							for (int e = k; e < k+orb_avail/2; ++e) {
// 								if (state[e] == '1') ++exchange;
// 							}
// 							slower[k] = '0';
// 							slower[k+orb_avail/2] = '1';
// 							Sm[sindex(bit2state(slower))] += pow(-1,exchange) * eigvec[j];
// 						}
// 					}
// 					else {
// 						// Calculate Sz
// 						Sz[j] += -0.5 * eigvec[j];
// 						// Calculate S+
// 						int exchange = 0;
// 						string sraise = state;
// 						if (state[k-orb_avail/2] == '0') {
// 							for (int e = k; e > k-orb_avail/2; --e) {
// 								if (state[e] == '1') ++exchange;
// 							}
// 							sraise[k] = '0';
// 							sraise[k-orb_avail/2] = '1';
// 							Sp[sindex(bit2state(sraise))] += pow(-1,exchange) * eigvec[j];
// 						}
// 					}
// 				}
// 			}
// 		} 
// 	}
// 	S2 = pow(ed::norm(Sz),2) + 0.5 * (pow(ed::norm(Sp),2) + pow(ed::norm(Sm),2));
// 	L2 = pow(ed::norm(Lz),2) + 0.5 * (pow(ed::norm(Lp),2) + pow(ed::norm(Lm),2));
// 	J2 = L2 + S2 + 2 * ed::dot(Lz,Sz) + ed::dot(Lm,Sp) + ed::dot(Lp,Sm);
// 	if (abs(L2) < 1e-7) L2 = 0;
// 	if (abs(S2) < 1e-7) S2 = 0;
// 	if (return_square) return vector<double>{J2,L2,S2};
// 	else {
// 		double L = (-1+sqrt(1+4*L2))/2;
// 		double S = (-1+sqrt(1+4*S2))/2;
// 		double J = (-1+sqrt(1+4*J2))/2;
// 		return vector<double>{J,L,S};
// 	}
// }

// void Hilbert::momentum_check(double* mat, double* eig, double* eigvec) {
// 	// Check Lz, L2, Sz, S2 for the matrix
// 	double* L2 = new double[hmat_size]{0};
// 	double* S2 = new double[hmat_size]{0};
// 	cout << setw(10) << "eigval" << setw(10) << "J2" << setw(10) << "L2" << setw(10) << "S2" << endl;
// 	for (int i = 0; i < hmat_size; ++i) {
// 		double* j_evec = new double[hmat_size]{0};
// 		for(int j = 0; j < hmat_size; ++j) j_evec[j] = eigvec[i*hmat_size+j];
// 		vector<double> LS = momentum(j_evec,true);
// 		cout << setw(10) << eig[i] << setw(10) << LS[0] << setw(10) << LS[1] << setw(10) << LS[2] << endl;
// 	}
// }

