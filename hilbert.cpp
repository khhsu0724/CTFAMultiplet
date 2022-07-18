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

Hilbert::Hilbert(string file_dir, vector<double*>& SC, double* FG, double* CF, double const& SO
				, bool HYB, string edge, bool is_ex): is_ex(is_ex), HYB_on(HYB), edge(edge) {
	read_from_file(file_dir);
	if (is_ex) {
		num_vh--;
		num_ch++;
	}
	Assign_Hash(FG,CF,SO);
	// Fix this so the code doesnt have to rely on a vector allocated
	vector<ulli> hspace = enum_hspace();
	// DEBUG
	// for (auto & h : hspace) {
	// 	cout << h << ", " << bitset<26>(h) << ", Hash: " << Hash(h) << ", Hashback: " << Hashback(Hash(h)) << endl;
	// 	if (h != Hashback(Hash(h))) cout << "error hashing" << endl;
	// }
	// for (auto & at : atlist) cout << at.atname << ", " << at.is_lig << endl;
	// DEBUG
	hsize = ed::choose(num_vorb,num_vh)*ed::choose(num_corb,num_ch);
	make_block(hspace);
}

vector<ulli> Hilbert::enum_hspace(ulli inc_val, ulli inc_core, int vmod, int cmod) {
	vector<ulli> val,core,hspace;
	ed::enum_states(val,num_vorb,num_vh+vmod,inc_val);
	ed::enum_states(core,num_corb,num_ch+cmod,inc_core);
	for (auto &c : core) for (auto &v : val) hspace.emplace_back(ed::add_bits(v,c,num_vorb,num_corb));
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
	// if snum <= 1, then the we don't need to compute this for loop?
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

void Hilbert::fill_hblk(double const& matelem, ulli const& lhs, ulli const& rhs) {
	// Find block index for lhs and rhs
	int lind = 0, rind = 0;

	// DEBUG
	// if (num_ch == 0) {
	// 	double lspin = 0, rspin = 0;
	// 	for (int j = num_corb/2; j < (num_vorb+num_corb)/2; ++j) {
	// 		if (lhs & 1<<j) lspin -= 0.5;
	// 		if (rhs & 1<<j) rspin -= 0.5;
	// 	}
	// 	lspin = 2 * lspin + 0.5 * num_vh;
	// 	rspin = 2 * rspin + 0.5 * num_vh;
	// 	if (lspin != rspin) cout << "lhs: " << bitset<16>(lhs) << ", rhs: " << bitset<16>(rhs) << endl;
	// }
	// DEBUG

	if (lind == rind) hblks[lind].ham[Hash(lhs)+hblks[lind].size*Hash(rhs)] += matelem;
	return;
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

int Hilbert::orbind(ulli s) {
	int i = 0;
	while (!(atlist[i].check & s)) ++i;
	return i;
}

int Hilbert::tot_site_num() {
	// Returns number of total sites
	int tsn = 1;
	for (const auto& s : sites) tsn *= s;
	return tsn; 
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
	hblks.emplace_back(std::move(Block(0,0,0,hsize)));
	return vb;
}

bool Hilbert::build_coordination(bool nh_read, int& tm_per_site, int& lig_per_site) {
	// This function builds information on atom relative 
	// positions and what sites they are on
	if (coord == "none") {
		tm_per_site = 1;
		if (edge == "K") throw invalid_argument("TM K edge unavailable");
	}
	else if (coord == "sqpl") {
		tm_per_site = 1;
		lig_per_site = 2;
	}
	this->at_per_site = tm_per_site + lig_per_site;
	if (nh_read) make_atlist(edge,tm_per_site,lig_per_site);
	return true;
}

void Hilbert::make_atlist(string edge, const int& tm_per_site, const int& lig_per_site) {
	int num_sites = tot_site_num(), atind = 0, siteind = 0;
	vector<int> dist = ed::distribute(num_vh,num_sites);
	for (int x = 0; x < this->sites[0]; ++x) {
	for (int y = 0; y < this->sites[1]; ++y) {
	for (int z = 0; z < this->sites[2]; ++z) {
		vector<int> site = {x,y,z};
		for (size_t tm = 0; tm < tm_per_site; ++tm) {
			if (edge == "L") atlist.emplace(atlist.begin(),Atom("2p",atind,3,2,0,site));
			atlist.emplace_back(Atom("3d",atind,3,2,dist[siteind],site));
			atind++;
		}
		for (size_t lig = 0; lig < lig_per_site; ++lig) {
			if (edge == "K") atlist.emplace(atlist.begin(),Atom("1s",atind,2,1,0,site));
			atlist.emplace_back(Atom("2p",atind,2,1,0,site));
			atind++;
		}
		siteind++;
	}}}
	this->num_at = atind;
	return;
}

void Hilbert::read_from_file(string file_dir) {
	string line;
	ifstream input_file(file_dir);
	try {
		if (input_file.is_open()) {
			bool read_cell = false, read_atoms = false, coord_built = false, nh_read = false;
			int atind = 0, cpsize = 0, tm_per_site = 0, lig_per_site = 0;
			while (getline(input_file,line)) {
				if (line[0] == '#') continue;
				if (line[0] == '/') {
					if (read_cell) coord_built = build_coordination(nh_read,tm_per_site,lig_per_site);
					read_cell = false;
					read_atoms = false;
					continue;
				}
				if (read_cell) {
					string p;
					bool skip = false;
					for (int s = 0; s < line.size() && !skip; ++s) {
						if (line[s] != ' ' && line[s] != '	' && line[s] != '=') p.push_back(line[s]);	
						else {
							transform(p.begin(),p.end(),p.begin(),::toupper);
							if (p == "COORDINATION") {
								string input = line;
								size_t ep, sp = input.find('"');
								if (sp != string::npos) {
								 	ep = input.find('"',++sp);
								 	if (ep != string::npos) input = input.substr(sp,ep-sp);
								 	else throw invalid_argument("No end quote");
								} else throw invalid_argument("Quotation needed for argument");
								transform(input.begin(),input.end(),input.begin(),::toupper);
								if (input == "SQUARE PLANAR" || input == "SQPL") coord = "sqpl";
								else if (input == "" || input == "NONE") coord = "none";
								else throw invalid_argument("unavailable coordination");
 								skip = true;
							}
							else if (p == "SITES") {
								vector<int> site_vec;
								size_t eqpos = line.find('=');
								string input = line.substr(eqpos+1,line.find('#')-eqpos-1);
								size_t sp = input.find(' ');
								while (sp!= string::npos) {
							        if (sp != 0) site_vec.push_back(stoi(input.substr(0,sp)));
							        input = input.substr(++sp);
							        sp = input.find(' ');
							    }
							    if (input.find_first_not_of(' ') != string::npos) site_vec.push_back(stoi(input));
							    if (site_vec.size() > 3) throw invalid_argument("incorrect site number");
							    for (size_t i = 0; i < site_vec.size(); ++i) sites[i] = site_vec[i];
								skip = true;
							}                    
							else if (p == "HOLES") {
								nh_read = true;
								size_t eqpos = line.find('=');
								string input = line.substr(eqpos+1,line.find('#')-eqpos-1);
								eqpos = input.find_first_not_of(' ');
								input = input.substr(eqpos,line.size()-eqpos);
								size_t sp = input.find(' ');
								if (sp == string::npos) num_vh = stoi(input);
								else {
								    num_vh = stoi(input.substr(0,sp));
							        if(input.substr(sp+1).find_first_not_of(' ') != string::npos) {
							        	throw invalid_argument("more than 1 argument");
							        }
								}
								skip = true;
							}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
						}
					}
					// call build coord and make sure sites & coordination are both read
				}
				if (read_atoms) {
					//gets how many line, edge determine what orbitals are included
					if (nh_read)  continue;
					if (!coord_built) throw invalid_argument("&Cell missing or needs to precede &Atoms");
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
								atin.add_orb(at,atind-1);
								atlist.emplace_back(atin);
								if (!atin.is_val) rotate(atlist.begin(),atlist.end()-1,atlist.end());
							}
						}
						else {
							if (at != "") {
								entry++;
								if (regex_match(at,regex(".*[a-zA-Z]+.*"))) {
									if (entry == 1) {
										num_at++;
										atin.atname = at;
									} else {
										atin.add_orb(at,atind-1);
										atlist.emplace_back(atin);
										if (!atin.is_val) rotate(atlist.begin(),atlist.end()-1,atlist.end());
									}
								} else atin.kp.emplace_back(stod(at));
							}
							at = "";
						}
					}
					if (atin.kp.size() != 3 && !skipline) throw invalid_argument("Invalid atom position");
				}
				if (line == "&CELL") {
					read_cell = true;
					read_atoms = false;
				}
				else if (line == "&ATOMS") {
					read_atoms = true;
					read_cell = false;
				}
			}
			input_file.close();
			int ind = 0;
			for (size_t i = 0; i < atlist.size(); ++i) {
				atlist[i].sind = ind;
				atlist[i].eind = ind + (2*atlist[i].l+1)-1;
				ind = atlist[i].eind + 1;
				if (atlist[i].is_val) {
					if (!nh_read) num_vh += atlist[i].num_h;
					num_vorb += 2*(2*atlist[i].l+1);
				} else {
					if (!nh_read) num_ch += atlist[i].num_h;
					num_corb += 2*(2*atlist[i].l+1);
					val_ind = atlist[i].eind + 1;
					val_ati = i + 1;
				}
				// size_t j = 0;
				// for (j = 0; j < i; ++j) if(ed::veccmp(atlist[i].kp,atlist[j].kp)) break;
			}
			if (num_vh > num_vorb) throw invalid_argument("too many holes from input");
			for (auto &at : atlist) {
				for (size_t i = at.sind; i <= at.eind; i++) at.check += (1 << i | 1 << (i+num_corb/2+num_vorb/2));
			}
			ind = 0;
			while(!atlist[ind].is_val) {
				int i = ind + 1;
				while (atlist[i].atind != atlist[ind].atind) ++i;
				atlist[ind].vind = i;
				atlist[i].cind = ind++;
			}
			//DEBUG
			// for (auto & at : atlist) {
			// 	cout << "atom ind: " << at.atind << ", val_n: " << at.val_n << ", val_l: " << at.val_l <<  ", n: " << at.n << ", l: " << at.l << ", num hole: " << at.num_h << endl;
			// 	cout << "atom sind: " << at.sind << ", eind: " << at.eind << ", vind: " << at.vind << ", cind: " << at.cind;
			// 	if (at.is_lig) cout << ", is ligand";
			// 	if (at.is_val) cout << ", is valence";
			// 	cout << ", sites: ";
			// 	for (auto s:at.site) cout << s << ",";
			// 	cout << "check: " << at.check << endl << endl;
			// }
			// DEBUG
		} else throw invalid_argument("Cannot open INPUT file");
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

// Here is a nice collection of hash functions
void Hilbert::Assign_Hash(double* FG, double* CF, double const& SO) {
	// Check the incoming parameter to see which hash function to use
	if (SO != 0) SO_on = true;
	if (!ed::is_zero_arr(CF,5)) CF_on = true;
	if (is_ex && !ed::is_zero_arr(FG,4)) CV_on = true;

	hashfunc = &Hilbert::norm_Hash;
	hbfunc = &Hilbert::norm_Hashback;
	return;
}

size_t Hilbert::norm_Hash(ulli s) {
	// Hash function that convert a state in bits to index
	size_t cind = 0, vind = 0;
	size_t cnt_c = 0, cnt_v = 0;
	size_t hc = num_corb/2, hv = num_vorb/2;
	for (size_t i = 0; i < hc; ++i)
		if (s & (1 << i)) cind += ed::choose(i,++cnt_c);
	for (size_t i = 0; i < hc; ++i)
		if (s & (1 << (i+hc+hv))) cind += ed::choose(i+hc,++cnt_c);
	for (size_t i = 0; i < hv; ++i)
		if (s & (1 << (i+hc))) vind += ed::choose(i,++cnt_v);
	for (size_t i = 0; i < hv; ++i)
		if (s & (1 << (i+num_corb+hv))) vind += ed::choose(i+hv,++cnt_v);
	return vind + cind * ed::choose(num_vorb,num_vh);
}

ulli Hilbert::norm_Hashback(size_t ind) {
	// Hash function that convert index to state in bitset
	size_t cind = ind / ed::choose(num_vorb,num_vh), ch = num_ch;
	size_t vind = ind % ed::choose(num_vorb,num_vh), vh = num_vh;
	ulli c = 0, v = 0;
	for (size_t i = num_vorb-1; i --> 0;) {
		if (vind >= ed::choose(i,vh)) {
			v |= (1 << i);
			vind -= ed::choose(i,vh--);
		}
	}
	for (size_t i = num_vorb-1; i --> 0;) {
		if (cind >= ed::choose(i,ch)) {
			c |= (1 << i);
			cind -= ed::choose(i,ch--);
		}
	}
	return ed::add_bits(v,c,num_vorb,num_corb);
}

