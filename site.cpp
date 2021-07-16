#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <utility>
#include "site.h"
#include "hilbert.h"
#include "multiplet.h"

using namespace std;
typedef vector<pair<int,int>> vpi;

Site::Site(string atname, bool add_ve, bool rem_ce, int ao_num) : 
	atname(atname), ao_num(ao_num), hsize(1), add_ve(add_ve), rem_ce(rem_ce) {
	try {
		if (atname == "Cu") {
			atnum = 29;
			orbs.resize(ao_num);
			orbs[0] = Hilbert(2,'p',6 - int(rem_ce));
			orbs[1] = Hilbert(3,'d',9 + int(add_ve));
		}
		else if (atname == "Ni") {
			atnum = 28;
			orbs.resize(ao_num);
			orbs[0] = Hilbert(2,'p',6 - int(rem_ce));
			orbs[1] = Hilbert(3,'d',8 + int(add_ve));
		}
		else if (atname == "Sc") {
			atnum = 21;
			orbs.resize(ao_num);
			orbs[0] = Hilbert(2,'p',2 - int(rem_ce));
			orbs[1] = Hilbert(3,'d',0 + int(add_ve));
		}
		else {
			throw invalid_argument("atom not coded yet");
		}
	} catch(const exception &ex) {
		std::cout << ex.what() << "\n";
		exit(0);
	}
	for (auto& o : orbs) hsize *= o.hmat_size;
	t_hsize = hsize;
	for (auto& o : orbs) o.t_hsize = hsize;
	build_index();
}

void Site::build_index() {
	// Reset site index
	for (auto &orb: orbs) orb.site_ind = vector<vector<int>>(orb.hmat_size);
	for (int j = 0; j < hsize; ++j) {
		// here each for loop is a state in the Hilbert Space
		int t = j;
		for (auto &orb: orbs) {
			int index = t % orb.hmat_size;
			t /= orb.hmat_size;
			orb.site_ind[index].emplace_back(j);
		}
	}
	return;
}

void Site::print_state(int pos) {
	try {
		if (pos >= hsize) throw invalid_argument("index larger than Hilbert Space");
		for (auto & orb: orbs) {
			int index = pos % orb.hmat_size;
			pos /= orb.hmat_size;
			cout << orb.get_n() << orb.get_orb() << ": " << orb.state2bit(orb.hmat[index]) << " ";
		}
		// cout << endl;o
	} catch (const exception &ex) {
		std::cout << ex.what() << "\n";
		exit(0);
	}
	return;
}

Hilbert& Site::get_hspace(int n, int l) {
	try {
		for (auto & orb: orbs) {
			if (orb.get_n() == n && orb.l == l) return orb;
		}
		throw range_error("specified orbital not included in the site");
	} catch (const exception &ex) {
		std::cout << ex.what() << "\n";
		exit(0);
	}
}

int Site::qn_get_state(int n, int l, int pos) {
	try {
		if (pos >= hsize) throw invalid_argument("index larger than Hilbert Space");
		for (auto & orb: orbs) {
			int index = pos % orb.hmat_size;
			if (orb.get_n() == n && orb.l == l) return orb.hmat[index];
			pos /= orb.hmat_size;
		}
		throw range_error("specified orbital not included in the site");
	} catch (const exception &ex) {
		std::cout << ex.what() << "\n";
		exit(0);
	}
	return -1;
}

double* Site::simplify_state(double* eigvec, int n, int l) {
	// Transforming eigenvector to orbital space to calculate momentum
	// Input size: site.hsize, return size: Hilbert:hspace
	try {
		for (auto & orb: orbs) {
			if (orb.get_n() == n && orb.l == l) {
				double* simp_vec = new double[orb.hmat_size]{0};
				for (int i = 0; i < orb.hmat_size; ++i) {
					for (auto &ind : orb.site_ind[i]) {
						simp_vec[i] += pow(eigvec[ind],2);
					}
					simp_vec[i] = sqrt(simp_vec[i]);
				}
				return simp_vec;
			}
		}
		throw range_error("specified orbital not included in the site");
	} catch (const exception &ex) {
		std::cout << ex.what() << "\n";
		exit(0);
	}
}

void cv_interacation(Site& core_site, Site& val_site, double* mat, double* FG) {
	Hilbert val = val_site.get_hspace(3,2), core = core_site.get_hspace(2,1);
	// electron language, loop through all electron combination
	int nvh = val.orb_avail - val.occ_num, nch = core.orb_avail - core.occ_num;
	if (nvh == 0 || nch == 0) return; 

	vpi mpair;
	for (int vml = -val.l; vml <= val.l; ++vml) {
		for (int vmr = -val.l; vmr <= val.l; ++vmr) {
			// valence electron pairs
			for (int cml = -core.l; cml <= core.l; ++cml) {
				for (int cmr = -core.l; cmr <= core.l; ++cmr) {
					if (vml+cml != vmr+cmr) continue;
					vpi spairs = {{1,1},{1,-1},{-1,1},{-1,-1}};
					// cout << vml << ", " << vmr << ", " << cml << ", " << cmr << endl;
					for (auto & s : spairs) {
						struct QN valqnl[1] = {{vml,s.first}};
						struct QN valqnr[1] = {{vmr,s.first}};
						struct QN coreqnl[1] = {{cml,s.second}};
						struct QN coreqnr[1] = {{cmr,s.second}};
						vpi val_entries = val.match_states(1, valqnl, valqnr);
						vpi core_entries = core.match_states(1, coreqnl, coreqnr);
						double meF = calc_U(gaunt(val.l,vml,val.l,vmr),
									gaunt(core.l,cmr,core.l,cml),FG,2*core.l+1);
						for (auto & ve : val_entries) {
							for (auto & ce : core_entries) {
								// cout << val.state2bit(ve.first) << ", " << val.state2bit(ve.second) << ", " << core.state2bit(ce.first) << ", " << core.state2bit(ce.second) << endl;
								double matelem = meF * val.Psign(valqnl,valqnr,ve.first,ve.second,1,1)
													 * core.Psign(coreqnl,coreqnr,ce.first,ce.second,1,1);

									// cout << "meF: " << meF * val.Psign(valqnl,valqnr,ve.first,ve.second,1,1)
									// 	 		* core.Psign(coreqnl,coreqnr,ce.first,ce.second,1,1) << endl;
								for (auto& indl : intersection(val.site_ind[val.sindex(ve.first)],core.site_ind[core.sindex(ce.first)])) {
									for (auto& indr : intersection(val.site_ind[val.sindex(ve.second)],core.site_ind[core.sindex(ce.second)])) {
										// cout << "index: " << indl << ", " << indr << endl;
										mat[indl+val.t_hsize*indr] += matelem;
									}
								}
							}
						}
						valqnl[0] = {vml,s.second};
						valqnr[0] = {vmr,s.first};
						coreqnl[0] = {cml,s.first};
						coreqnr[0] = {cmr,s.second};
						val_entries = val.match_states(1, valqnl, valqnr);
						core_entries = core.match_states(1, coreqnl, coreqnr);
						double meG = calc_U(gaunt(core.l,cml,val.l,vmr),
									gaunt(core.l,cmr,val.l,vml),FG,core.l+val.l+1);
						for (auto & ve : val_entries) {
							for (auto & ce : core_entries) {
								// cout << val.state2bit(ve.first) << ", " << val.state2bit(ve.second) << ", " << core.state2bit(ce.first) << ", " << core.state2bit(ce.second) << endl;
								double matelem = -meG * val.Psign(valqnl,valqnr,ve.first,ve.second,1,1)
													  * core.Psign(coreqnl,coreqnr,ce.first,ce.second,1,1); // Traverse one extra time

								// cout << ", meG: " << -meG * val.Psign(valqnl,valqnr,ve.first,ve.second,1,1)
								// 					  * core.Psign(coreqnl,coreqnr,ce.first,ce.second,1,1) << endl;
								for (auto& indl : intersection(val.site_ind[val.sindex(ve.first)],core.site_ind[core.sindex(ce.first)])) {
									for (auto& indr : intersection(val.site_ind[val.sindex(ve.second)],core.site_ind[core.sindex(ce.second)])) {
										// cout << "index: " << indl << ", " << indr << endl;
										mat[indl+val.t_hsize*indr] += matelem;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	// Shift energy for particle hole symmetry
	double tr = trace(mat,core_site.t_hsize,core_site.t_hsize);
	double eshift;
	if ((val.occ_num > val.orb_avail/2) && (core.occ_num > core.orb_avail/2)) {
		double p = val.occ_num * core.occ_num;
		eshift = tr*(p-(val.orb_avail-val.occ_num)*(core.orb_avail-core.occ_num))/p/val_site.t_hsize;
	}
	else if (val.occ_num > val.orb_avail/2) eshift = val.pheshift(tr,1);
	else if (core.occ_num > core.orb_avail/2) eshift = core.pheshift(tr,1);
	else eshift = 0;
	for (int i = 0; i < core_site.t_hsize; ++i) mat[i+core_site.t_hsize*i] -= eshift;
}


void populate_hspace(Site* sites, int site_num, double* mat, double* SC, double tenDQ, double lambda) {
	for (int i = 0; i < site_num; ++i) {
		for (auto & orb : sites[i].orbs) {
			calc_coulomb(orb,mat,SC);
			if (lambda != 0 && orb.l == 1) calc_SO(orb,mat,lambda); // Enforce Spin orbit coupling only on p
			if (tenDQ != 0 && orb.l == 2) calc_CF(orb,mat,tenDQ,6); // Enforce crystal field only on d
		}
	}
	return;
}