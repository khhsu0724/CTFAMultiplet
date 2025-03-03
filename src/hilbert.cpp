#include <fstream>
#include <regex>
#include <cmath>
#include <bitset>
#include "hilbert.hpp"

using namespace std;
// This file contains all operations related to Hilbert space

Hilbert::Hilbert(string file_dir, const HParam& hparam, string edge, bool is_ex): 
					is_ex(is_ex), edge(edge) {
	read_from_file(file_dir);
	num_vh -= is_ex;
	num_ch += is_ex;
	HYB_on = (hparam.tpd || hparam.tpp || hparam.MLdelta) && hparam.HYB;
	SO_on = !ed::is_zero_arr(hparam.SO,2);
	CF_on = !ed::is_zero_arr(hparam.CF,5);
	CV_on = (is_ex && !ed::is_zero_arr(hparam.FG,4));
	cluster->set_print_site(hparam.print_site_occ);
	assign_phonon(file_dir);
	// Calculate Hilbert space size
	if (hparam.block_diag && !hparam.SO[1]) BLOCK_DIAG = !(edge == "L" && is_ex) || !hparam.SO[0];
	int hv = num_vorb/2, hc = num_corb/2;
	if (this->have_phonon()) phonon_hsize = phonon_modes->get_pn_hsize();
	this->hsize = ed::choose(num_vorb,num_vh) * ed::choose(num_corb,num_ch) * phonon_hsize;
	cout << "phonon hilbert space size: " << phonon_hsize << ", hsize: " << this->hsize << endl;
	if (num_vh > num_vorb) throw invalid_argument("too many holes from input");
	size_t bfind = 0;
	if (BLOCK_DIAG) {
		// Sz ordered blocks, when spin orbit coupling is not on
		// cout << "BLOCK DIAG" << endl;
		SO_on = false;
		hashfunc = &Hilbert::sz_Hash;
		hbfunc = &Hilbert::sz_Hashback;
		int max_2csz = hc - abs(hc-num_ch);
		int max_2sz = hv - abs(hv-num_vh) + max_2csz;
		int core_comb = ed::choose(num_corb,num_ch);
		hblks.reserve(max_2sz+1);
		for (int sz = -max_2sz; sz <= max_2sz; sz+=2) {
			vector<ulli> rank(core_comb,-1);
			ulli c = (BIG1 << max_2csz) - 1; // This swaps between core/hole
			size_t blksize = 0;
			for (int i = 0; i < core_comb; ++i) {
				int vsd = (-sz+max_2sz)/2 - ed::count_bits(c%(BIG1<<hc));
				int vsu = (sz+max_2sz)/2 - ed::count_bits(c>>hc);
				c = ed::next_perm(c);
				if (vsd < 0 || vsu < 0) continue;
				rank[i] = blksize/phonon_hsize; // Rank is fermionic hashing
				blksize += ed::choose(hv,vsd) * ed::choose(hv,vsu) * phonon_hsize;
			}
			hblks.emplace_back(double(sz)/2,0,0,blksize,bfind,rank);
			bfind += blksize;
		}
	} else if (false) {
		// Jz ordered blocks, should not be used
		hashfunc = &Hilbert::jz_Hash;
		hbfunc = &Hilbert::jz_Hashback;
		int min_2jz = 0, max_2jz = 0;
		vector<ulli> hspace = enum_hspace();
		vector<int> jz2_arr(hspace.size(),0);
		for (auto & h : hspace) {
			int jz2 = 0;
			jz2 +=  ed::count_bits(h/(BIG1<<(hc+hv+1))) - ed::count_bits(h%(BIG1<<(hc+hv+1)));
			for (auto & at : atlist) {
				for (int i = at.sind; i <= at.eind; ++i) {
					if (h & (BIG1<<i)) jz2 += 2*(i-at.sind-at.l);
					if (h & (BIG1<<(i+hv+hc))) jz2 += 2*(i-at.sind-at.l);
				}
			}
			jz2_arr[&h-&hspace[0]] = jz2;
			if (jz2 < min_2jz) min_2jz = jz2;
			if (jz2 > max_2jz) max_2jz = jz2;
		}
		// Calculate block size
		vector<int> blk_size((max_2jz-min_2jz)/2+1,0);
		for (int i = 0; i < hspace.size(); ++i) blk_size[(jz2_arr[i]-min_2jz)/2] += phonon_hsize;
		// Make Blocks
		hblks.reserve(blk_size.size());
		for (int j = min_2jz; j <= max_2jz; j += 2) {
			hblks.emplace_back(0,double(j)/2,0,blk_size[(j-min_2jz)/2]);
			bfind += blk_size[(j-min_2jz)/2];
		}
		for (int b = 0; b < hblks.size(); ++b) {
			vector<ulli> rank;
			for (int i = 0; i < hspace.size(); ++i) {
				if ((jz2_arr[i]-min_2jz)/2 == b) rank.push_back(hspace[i]);
			}
			if (is_ex) std::sort(rank.begin(), rank.end(), std::greater<ulli>());
			hblks[b].rank = std::move(rank);
		}		
	}
	else {
		cout << "NORM DIAG" << endl;
		hashfunc = &Hilbert::norm_Hash;
		hbfunc = &Hilbert::norm_Hashback;
		hblks.emplace_back(0,0,0,this->hsize);
	}
	for (auto & hblk : hblks) hblk.set_phonon_hsize(phonon_hsize);
	// DEBUG: repeating n times per phonon mode
	// vector<ulli> hspace = enum_hspace();
	// if (is_ex) cout << "excited state" << endl;
	// else cout << "groud state" << endl;
	// for (auto & h : hspace) {
	// 	auto i = Hash(h);
	// 	Hashback(i);
	// 	cout << h << ", " << bitset<34>(h) << ", Hash: " << i.first << ", " << i.second << ", Hashback: " << Hashback(Hash(h)) << endl;
	// 	if (h != Hashback(Hash(h))) cout << "error hashing" << endl;
	// }
	// cout << "ATOMS: " << endl;
	// for (auto & at : atlist) cout << at.atname << ", " << at.is_lig << endl;
	// exit(0);
	// DEBUG
	return;
}

Hilbert::Hilbert(const Hilbert &hilbs, int vh_mod, int pnhsize) {
	// Copy constructor, vh_mod can be used for dummy hilbert spaces
	num_vh = hilbs.num_vh + vh_mod;
	num_ch = hilbs.num_ch;
	num_vorb = hilbs.num_vorb;
	num_corb = hilbs.num_corb;
	SO_on = hilbs.SO_on;
	CF_on = hilbs.CF_on;
	CV_on = hilbs.CV_on;
	HYB_on = hilbs.HYB_on;
	BLOCK_DIAG = hilbs.BLOCK_DIAG;
	num_at = hilbs.num_at;
	at_per_site = hilbs.at_per_site;
	val_ati = hilbs.val_ati;
	val_ind = hilbs.val_ind;
	coord = hilbs.coord;
	edge = hilbs.edge;
	atlist = hilbs.atlist;
	sites = hilbs.sites;
	// I think I should copy the phonons here
	if (pnhsize == 0) phonon_hsize = hilbs.phonon_hsize;
	else phonon_hsize = pnhsize;
	if (!vh_mod) cout << "USING COPY CONSTRUCTOR IS NOT ADVISED" << endl;
	hsize = ed::choose(num_vorb,num_vh) * ed::choose(num_corb,num_ch) * phonon_hsize;

	assign_cluster(coord);
	// This is not great, we should be able to specify Hash function?
	hashfunc = &Hilbert::norm_Hash;
	hbfunc = &Hilbert::norm_Hashback;
	// TODO: complete copy constructor for no vh_mod
	return;
}

// Enumerate all possible bit state in the fermionic Hilbert space, 
// with the inclusion any existing core/valence state 
// (vmod changes many valence states are included)
// Important: inc_val should be input as a binary with !!2*nvo!! digits
//.           same with inc_core. The code will split these up and 
//			  combine these to a binary that can is in the Hilbert space
//			  See interaction.cpp->holestein for a simple example
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
		s += BIG1 << ind;
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
	// Find block index for lhs and rhs (lhs and rhs are the fermionic part)
	// repeat this for the bosonic part
	bindex lind = Hash(lhs);
	bindex rind = Hash(rhs);
	if (lind.first == rind.first) {
		int fermion_hsize = hblks[lind.first].fermion_hsize;
		for (int i = 0; i < hblks[lind.first].phonon_hsize; ++i) {
			hblks[lind.first].ham->fill_mat(lind.second+i*fermion_hsize,
											rind.second+i*fermion_hsize,matelem);
		}
	} else throw out_of_range("invalid block matrix element entry");
	return;
}

void Hilbert::fill_hblk_pn(double const& matelem, bindex const& lind_pn, bindex const& rind_pn,
							ulli const& lhs_fm, ulli const& rhs_fm) {
	// Find block index for lhs and rhs for phonon part
	// lhs_pn and rhs_pn is for phonon, lhs_fm and rhs_fm is for fermion
	// !!!!! PHONON BLOCKS ARE NOT IMPLEMENTED !!!!!
	int ipn = lind_pn.second;
	int jpn = rind_pn.second;
	if (lhs_fm != 0 && rhs_fm != 0) {
		bindex lind_fm = Hash(lhs_fm);
		bindex rind_fm = Hash(rhs_fm);
		// print_bits(lhs_fm);
		// print_bits(rhs_fm);
		int ifm = lind_fm.second;
		int jfm = rind_fm.second;
		// cout << "SIZE: " << hblks[lind_fm.first].size << "," << hblks[lind_fm.first].fermion_hsize << endl;
		if (lind_fm.first == rind_fm.first) {
			// cout << "pnind: " << ipn << "," << jpn << ", fmind: " << ifm << "," << jfm << endl;
			int fermion_hsize = hblks[lind_fm.first].fermion_hsize;
			hblks[lind_fm.first].ham->fill_mat(ifm+ipn*fermion_hsize,
												jfm+jpn*fermion_hsize,matelem);
		} else throw out_of_range("invalid block matrix element entry");

	} else {
		// No fermionic term, go through all fermionic basis state
		if (lhs_fm != 0 || rhs_fm != 0)
			throw out_of_range("invalid fermionic basis state for matrix element entry");
		for (auto& blk : hblks) {
			int fermion_hsize = blk.fermion_hsize;
			for (int ifm = 0; ifm < fermion_hsize; ++ifm) {
				blk.ham->fill_mat(ifm+ipn*fermion_hsize,ifm+jpn*fermion_hsize,
									matelem);
			}
		}
	}
	return;
}

void Hilbert::print_bits(ulli state) {
	// const int bitnum = (num_corb+num_vorb)*2;
	cout << bitset<22>(state) << endl; // This is a bit crude
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

// TODO: avoid duplication function as above
double Hilbert::Fsign(ulli* op, ulli state, int opnum) {
	int p = 0;
	for (size_t i = 0; i < opnum; ++i) {
		ulli o = *(op+i); 
		if (!(state & o)) return 0;
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

vector<double> Hilbert::get_all_eigval(bool is_err) {
	vector<double> all_eig(hsize,0);
	int last_ind = 0;
	for (auto & blk : hblks) {
		if (blk.eig == nullptr) {
			if (is_err) throw out_of_range("invalid access of eigval");
		} else {
			for (size_t i = 0; i < blk.nev; ++i) {
				all_eig[last_ind] = blk.eig[i];
				last_ind++;
			}
		}
	}
	all_eig.resize(last_ind);
	return all_eig;
}

vector<ulli> Hilbert::get_hashback_list(size_t blk_ind) {
	// Get list of bit represented state in specified block
	int blksize = hblks[blk_ind].size;
	vector<ulli> hblist(blksize);
	#pragma omp parallel for shared(hblist)
	for (size_t i = 0; i < blksize; ++i) {
		hblist[i] = this->Hashback(bindex(blk_ind,i));
	}
	return hblist;
}


// Functions for file parsing
void Hilbert::assign_cluster(string input) {
	transform(input.begin(),input.end(),input.begin(),::toupper);
	if (input == "" || input == "NONE" || input == "ION") {
		coord = "ion";
		cluster = new Ion(edge);
	} else if (input == "SQUARE PLANAR" || input == "SQPL") {
		coord = "sqpl"; 
		cluster = new SquarePlanar(edge);
	} else if (input == "OCTAHEDRAL" || input == "OCT") {
		coord = "oct"; 
		cluster = new Octahedral(edge);
	} else throw invalid_argument("unavailable coordination");
	return;
}

void Hilbert::assign_phonon(string file_dir) {
	// Ideally this should be read from a separate &PHONON tag
	if (cluster == NULL or cluster->lig_per_site == 0) {
		cout << "Cannot have phonon modes" << endl;
		return;
	}
	phonon_modes = new PhononModes(this->cluster);
	phonon_modes->read_phonon_from_input(file_dir);
	phonon_modes->hashfunc = &PhononModes::norm_Hash;
	phonon_modes->hbfunc = &PhononModes::norm_Hashback;
	return;
}

void Hilbert::read_from_file(string file_dir) {
	string line;
	ifstream input_file(file_dir);
	try {
		if (input_file.is_open()) {
			bool read_cell = false;
			int nh_read = 0;
			int atind = 0, cpsize = 0, tm_per_site = 0, lig_per_site = 0;
			while (getline(input_file,line)) {
				if (line[0] == '#') continue;
				if (line[0] == '/') {
					read_cell = false;
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
								assign_cluster(input);
 								skip = true;
							}
							if (p == "HYBMAT") {
								string input = line;
								size_t ep, sp = input.find('"');
								if (sp != string::npos) {
								 	ep = input.find('"',++sp);
								 	if (ep != string::npos) input = input.substr(sp,ep-sp);
								 	else throw invalid_argument("No end quote");
								} else throw invalid_argument("Quotation needed for argument");
								// Read in file
								if (input.empty()) cout << "No input matrix found, default to &CONTROL" << endl;
								cout << "Found input hybridization matrix, tpd/tpp/MLCT ignored" << endl;
								inp_hyb_file = input;
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
							else if (p == "HOLES" || p == "TOT_HOLES" || p == "HOLES_PER_SITE") {
								if (nh_read) cout << "WARNING: Number of holes already read, overwritting" << endl;
								if (p == "HOLES_PER_SITE") nh_read = 2;
								else nh_read = 1;
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
				}
				if (line == "&CELL") read_cell = true;
			}
			input_file.close();
			if (!nh_read) throw invalid_argument("Number of holes needed for input");
			else if (nh_read == 2) num_vh *= tot_site_num();
			if (cluster == NULL) assign_cluster("");
			cluster->read_inp_tmat(inp_hyb_file);
			this->at_per_site = cluster->at_per_site();
			this->num_at = cluster->at_per_site() * tot_site_num();
			cluster->set_num_sites(tot_site_num());
			cluster->make_atlist(atlist,num_vh,sites);
			int ind = 0;
			for (size_t i = 0; i < atlist.size(); ++i) {
				atlist[i].sind = ind;
				atlist[i].eind = ind + (2*atlist[i].l+1)-1;
				ind = atlist[i].eind + 1;
				if (atlist[i].is_val) num_vorb += 2*(2*atlist[i].l+1);
				else {
					num_corb += 2*(2*atlist[i].l+1);
					val_ind = atlist[i].eind + 1;
					val_ati = i + 1;
				}
			}
			// TEMPORARY SAFE GAURD 
			if (num_corb+num_vorb > 64) 
				throw invalid_argument("Cannot have more than 32 orbitals for now");
			// TEMPORARY SAFE GAURD 
			for (auto &at : atlist) {
				for (size_t i = at.sind; i <= at.eind; i++) {
					at.check += (BIG1 << i | BIG1 << (i+num_corb/2+num_vorb/2));
				}
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
			// 	cout << " check: " << bitset<34>(at.check) << endl << endl;
			// }
			// cout << "num_ch: " << num_ch << ", num_vh: " << num_vh << ", num_corb: " << num_corb << ", num_vorb: " << num_vorb << endl;
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
// TODO: A function to label quantum number for states

// Here is a nice collection of hash functions

bindex Hilbert::norm_Hash(ulli s) {
	// Hash function that convert a state in bits to index
	size_t cind = 0, vind = 0;
	size_t cnt_c = 0, cnt_v = 0;
	size_t hc = num_corb/2, hv = num_vorb/2;
	for (size_t i = 0; i < hc; ++i)
		if (s & (BIG1 << i)) cind += ed::choose(i,++cnt_c);
	for (size_t i = 0; i < hc; ++i)
		if (s & (BIG1 << (i+hc+hv))) cind += ed::choose(i+hc,++cnt_c);
	for (size_t i = 0; i < hv; ++i)
		if (s & (BIG1 << (i+hc))) vind += ed::choose(i,++cnt_v);
	for (size_t i = 0; i < hv; ++i)
		if (s & (BIG1 << (i+num_corb+hv))) vind += ed::choose(i+hv,++cnt_v);
	return bindex(0, vind+cind*ed::choose(num_vorb,num_vh));
}

ulli Hilbert::norm_Hashback(bindex ind) {
	// Hash function that convert index to state in bitset
	ind.second = ind.second % hblks[ind.first].fermion_hsize;
	size_t edchoose = ed::choose(num_vorb,num_vh);
	size_t cind = ind.second / edchoose, ch = num_ch;
	size_t vind = ind.second % edchoose, vh = num_vh;
	ulli c = 0, v = 0;
	for (size_t i = num_vorb; i --> 0;) {
		if (vind >= ed::choose(i,vh)) {
			v |= (BIG1 << i);
			vind -= ed::choose(i,vh--);
		}
	}
	for (size_t i = num_corb; i --> 0;) {
		if (cind >= ed::choose(i,ch)) {
			c |= (BIG1 << i);
			cind -= ed::choose(i,ch--);
		}
	}
	return ed::add_bits(v,c,num_vorb,num_corb);
}

bindex Hilbert::sz_Hash(ulli s) {
	// Hash function that uses sz as main quantum number
	size_t cind = 0, vsdind = 0, vsuind = 0;
	size_t cnt_c = 0, cnt_v = 0;
	size_t hc = num_corb/2, hv = num_vorb/2;
	int nsd = ed::count_bits(s % (BIG1 << (hc+hv)));
	int nsu = num_vh + num_ch - nsd;
	int max_2sz = int(2*hblks.back().get_sz());
	size_t blk_ind = (max_2sz-nsd+nsu)/2;
	for (size_t i = 0; i < hc; ++i)
		if (s & (BIG1 << i)) cind += ed::choose(i,++cnt_c);
	for (size_t i = 0; i < hc; ++i)
		if (s & (BIG1 << (i+hc+hv))) cind += ed::choose(i+hc,++cnt_c);
	for (size_t i = 0; i < hv; ++i)
		if (s & (BIG1 << (i+hc))) vsdind += ed::choose(i,++cnt_v);
	size_t vsdchoose = ed::choose(hv,cnt_v);
	cnt_v = 0;
	for (size_t i = 0; i < hv; ++i)
		if (s & (BIG1 << (i+num_corb+hv))) vsuind += ed::choose(i,++cnt_v);
	return bindex(blk_ind,hblks[blk_ind].rank[cind]+vsdind+vsuind*vsdchoose);
}

ulli Hilbert::sz_Hashback(bindex ind) {
	// Hashback function that uses sz as main QN
	ind.second = ind.second % hblks[ind.first].fermion_hsize;
	int max_2sz = int(2*hblks.back().get_sz());
	auto& blk = hblks[ind.first];
	auto r = std::find_if(blk.rank.rbegin(), blk.rank.rend(),
				[&](size_t e){return (e >= 0) && (ind.second >= e);});
	if (r == blk.rank.rend()) r = blk.rank.rend() - 1;
	size_t vind = ind.second - *r, hv = num_vorb/2;
	size_t cind = blk.rank.rend() - r - 1, ch = num_ch;
	ulli c = 0, vsd = 0, vsu = 0;
	for (size_t i = num_corb; i --> 0;) {
		if (cind >= ed::choose(i,ch)) {
			c |= (BIG1 << i);
			cind -= ed::choose(i,ch--);
		}
	}	
	size_t nvhsd = (num_vh+num_ch+max_2sz)/2-ind.first-ed::count_bits(c%(BIG1<<(num_corb/2)));
	size_t nvhsu = num_vh - nvhsd, vsdchoose = ed::choose(hv,nvhsd);
	size_t vsdind = vind % vsdchoose, vsuind = vind / vsdchoose;
	for (size_t i = hv; i --> 0;) {
		if (vsdind >= ed::choose(i,nvhsd)) {
			vsd |= (BIG1 << i);
			vsdind -= ed::choose(i,nvhsd--);
		}
		if (vsuind >= ed::choose(i,nvhsu)) {
			vsu |= (BIG1 << i);
			vsuind -= ed::choose(i,nvhsu--);
		}
	}
	return ed::add_bits(vsd|(vsu<<hv),c,num_vorb,num_corb);
}

bindex Hilbert::jz_Hash(ulli s) {
	// Hash function that uses jz as main QN
	int hv = num_vorb/2, hc = num_corb/2;
	int min_2jz = int(2*hblks[0].get_jz()), jz2 = 0;
	jz2 +=  ed::count_bits(s/(BIG1<<(hc+hv+1))) - ed::count_bits(s%(BIG1<<(hc+hv+1)));
	for (auto & at : atlist) {
		for (int i = at.sind; i <= at.eind; ++i) {
			if (s & (BIG1<<i)) jz2 += 2*(i-at.sind-at.l);
			if (s & (BIG1<<(i+hv+hc))) jz2 += 2*(i-at.sind-at.l);
		}
	}
	size_t bind = (jz2-min_2jz)/2;
	return bindex(bind,ed::binary_search(hblks[bind].rank,s,
						std::less<ulli>(),std::equal_to<ulli>()));
}

ulli Hilbert::jz_Hashback(bindex ind) {
	// Hashback function that uses jz as main QN
	ind.second = ind.second % hblks[ind.first].fermion_hsize;
	return hblks[ind.first].rank[ind.second];

}
