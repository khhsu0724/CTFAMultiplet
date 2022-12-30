#include <complex>
#include <bitset>
#include "multiplet.hpp"
#include "photon.hpp"

#define PRINT_TOL 1E-4

using namespace std;
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
		std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
		initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp declare reduction(vec_dcomp_plus : std::vector<dcomp> : \
		std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<dcomp>())) \
		initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

// Contains photon based spectroscopy
string pol_str(const vecd& pvec) {
	if (pvec[0]) return "X";
	if (pvec[1]) return "Y";
	if (pvec[2]) return "Z";
	return "";
}

dcomp proj_pvec(int ml, const vecd& pvec) {
	// Polarization Vector Projection to Spherical Harmonics
	try {
		if (ml == 1) return dcomp(-pvec[0]/sqrt(2),pvec[1]/sqrt(2));
		else if (ml == 0) return dcomp(pvec[2],0);
		else if (ml == -1) return dcomp(pvec[0]/sqrt(2),pvec[1]/sqrt(2));
		else throw invalid_argument("invalid ml for polarization vector projection");
		
	} catch (const exception &ex) {
		std::cout << ex.what() << "\n";
		exit(0);
	}
}

vecd occupation(Hilbert& hilbs, const vector<bindex>& si, bool is_print) {
	// Calculation occupation of orbitals, only valid when matrix is diagonalized
	for (auto& blk : hilbs.hblks) if (blk.eigvec == nullptr) throw runtime_error("matrix not diagonalized for occupation");
	int nvo = hilbs.cluster->vo_persite, nco = hilbs.cluster->co_persite;
	vecc U = hilbs.cluster->get_seph2real_mat();
	vector<double> occ_totc(nvo,0);
	Hilbert minus_1vh(hilbs,-1);
	for (auto &s  : si) {
		Block& blk = hilbs.hblks[s.first];
		// Pick an operator
		for (size_t c = 0; c < nvo; ++c) {
			vector<complex<double>> wvfnc(minus_1vh.hsize,0);
			for (size_t ci = 0; ci < nvo; ++ci) {
				if (abs(U[c*nvo+ci]) < TOL) continue;
				for (size_t j = 0; j < blk.size; ++j) {
					double sj = blk.eigvec[s.second*blk.size+j];
					if (abs(sj) < TOL) continue;
					// Operate on spin up and spin down
					for (int sud = nco; sud <= 2*nco+nvo; sud += nco+nvo) {
						ulli op_state = hilbs.Hashback(bindex(s.first,j)), op;
						if (1 << (ci+sud) & op_state) op = 1<<(ci+sud);
						else continue;
						bindex ind = minus_1vh.Hash(op_state-op);
						int fsgn = hilbs.Fsign(&op,op_state,1);
						wvfnc[ind.second] += fsgn*sj*U[c*nvo+ci];
					}
				}
			}
			for (auto &w : wvfnc) occ_totc[c] += abs(pow(w,2))/si.size();
		}
	}
	double tot_h = 0;
	for (auto &h : occ_totc) tot_h += h;
	if (is_print) hilbs.cluster->print_eigstate(occ_totc);
	return occ_totc;
}

vector<double> wvfnc_weight(Hilbert& hilbs, const vector<bindex>& si, int ligNum, bool print) {
	// Returns fraction of degenerate states with their ligand occupation
	vector<double> wvfnc(hilbs.hsize), dLweight(ligNum,0);
	for (auto& s : si) {
		auto &blk = hilbs.hblks[s.first];
		for (size_t i = 0; i < blk.size; ++i) {
			// Calculate Wavefunction
			if (abs(blk.eigvec[s.second*blk.size+i]) < TOL) continue;
			wvfnc[blk.f_ind+i] += pow(blk.eigvec[s.second*blk.size+i],2);
			ulli state = hilbs.Hashback(bindex(s.first,i));
			int Lcnt = 0;
			// This needs to account geometry
			for (size_t j = 5; j < 11; ++j) {
				if ((1<<(j+hilbs.num_corb/2)) & state) Lcnt++;
				if ((1<<(j+hilbs.num_corb+hilbs.num_vorb/2)) & state) Lcnt++;
			}
			if (Lcnt < ligNum) dLweight[Lcnt] += pow(blk.eigvec[s.second*blk.size+i],2)/si.size();
		}
	}
	if (print) {
		cout << "Ground State composition";
		for (size_t i = 0; i < ligNum; ++i) {
			cout << ", d" << 10-hilbs.num_vh+i;
			if (i == 1) cout << "L: ";
			else if (i > 1) cout << "L" << i << ": ";
			else cout << ": ";
			cout  << fixed << setprecision(5) << dLweight[i];
		}
		cout << endl;
	}
	return dLweight;
}

double effective_delta(Hilbert& hilbs, int ligNum, bool is_print) {
	// Check effective delta, needs t = 0
	double d, dL;
	vector<bool> firstLh(ligNum,true);
	vector<double> distinct = ed::printDistinct(hilbs.get_all_eigval(false),0.0,hilbs.hsize,is_print);
	for (size_t e = 0; e < distinct.size(); ++e) {
		vector<bindex> gei;
		for (size_t i = 0; i < hilbs.hblks.size(); ++i) {
			if (hilbs.hblks[i].eig == nullptr) continue;
			for (size_t j = 0; j < hilbs.hblks[i].size; ++j) {
				if (abs(hilbs.hblks[i].eig[j]-distinct[e]) < TOL) {
					gei.push_back({i,j});
				}
			}
		}
		vector<double> dLweight = wvfnc_weight(hilbs,gei,ligNum,false);
		for (size_t i = 0; i < ligNum; ++i) {
			if (abs(dLweight[i]-1) < TOL && firstLh[i]) {
				if (i == 0) d = distinct[e];
				if (i == 1) dL = distinct[e];
				firstLh[i] = false;
				if (is_print) {
					cout << "Energy: " << distinct[e] << ", degeneracy: " << gei.size() << endl;
					wvfnc_weight(hilbs,gei,ligNum,true);
				}
			}
		}
	}
	if (is_print) cout << "Effective delta: " << dL - d << endl;
	return dL - d;
}

void state_composition(Hilbert& hilbs, const vector<bindex>& si, size_t top) {
	// prints out many-body wavefunction composition of vector of eigen-states
	vector<double> wvfnc(hilbs.hsize);
	for (auto& s : si) {
		auto &blk = hilbs.hblks[s.first];
		for (size_t i = 0; i < blk.size; ++i) {
			if (abs(blk.eigvec[s.second*blk.size+i]) < TOL) continue;
			wvfnc[blk.f_ind+i] += pow(blk.eigvec[s.second*blk.size+i],2);
		}
	}
	double wv_sum = 0;
	// Find index of state in wavefunction with top weight
	vector<size_t> indices = ed::top_n(wvfnc,top);
	for (auto& blk : hilbs.hblks) {
		for (size_t i = 0; i < blk.size; ++i) {
			wv_sum += wvfnc[blk.f_ind+i]/si.size();
			if (abs(wvfnc[blk.f_ind+i]) < TOL) continue;
			for (size_t ind = 0; ind < top; ++ind) {
				if ((blk.f_ind+i) == indices[ind]) {
					cout << "weight: " << wvfnc[blk.f_ind+i]*100.00/si.size() << "%, state: ";
					cout << bitset<26>(hilbs.Hashback(bindex(&blk-&hilbs.hblks[0],i))) << endl;
					break;
				}
			}
		}
	}
	return;
}

void basis_overlap(Hilbert& GS, Hilbert& EX, bindex inds, vector<blapIndex>& blap, 
					const PM& pm, bool pvout) {
	// Calculate basis state overlap base on the indices of the blocks
	size_t gbi = inds.first, exi = inds.second;
	blap.clear();
	blap.reserve(EX.hblks[exi].size*2); // rough guess
	vecd pvec = pvout ? pm.pvout : pm.pvin;
	int cl = GS.atlist[0].l, vl = GS.atlist[GS.atlist[0].vind].l;
	int half_orb = (EX.num_vorb+EX.num_corb)/2;
	#pragma omp parallel for shared(blap) collapse(2)
	for (size_t g = 0; g < GS.hblks[gbi].size; g++) {
		for (size_t e = 0; e < EX.hblks[exi].size; e++) {
			ulli gs = GS.Hashback(bindex(gbi,g)), exs = EX.Hashback(bindex(exi,e));
			ulli ch = exs - (gs & exs), vh = gs - (gs & exs);
			int coi = EX.orbind(ch), voi = GS.atlist[coi].vind;
			if (!ed::is_pw2(vh) || !GS.atlist[voi].contains(vh)) continue;
			QN chqn = EX.atlist[coi].fast_qn(ch,half_orb,coi);
			QN vhqn = GS.atlist[voi].fast_qn(vh,half_orb,voi);
			if (chqn.spin != vhqn.spin || abs(vhqn.ml-chqn.ml) > 1) continue;
			dcomp blap_val = gaunt(cl,chqn.ml,vl,vhqn.ml)[1] 
				* GS.Fsign(&vh,gs,1) * EX.Fsign(&ch,exs,1) * proj_pvec(vhqn.ml-chqn.ml,pvec);
			if (blap_val != dcomp(0.0,0.0)) {
				#pragma omp critical
				{
					blap.emplace_back(g,e,blap_val);
				}
			}
		}
	}
	return;
}

void XAS_peak_occupation(Hilbert& GS, Hilbert& EX, vecd const& peak_en, vecd const& energy, 
						vecd const& intensity, vector<bindex> const& gsi, double ref_en, 
						string mode, bool ref_gs) {
	// Calculate peak occupation using energy
	// If mode = top, it will use the peak_en as the top n number of peaks
	// TODO: Range of peak energy
	vecd peak_en_copy;
	vector<size_t> indices;
	if (mode == "list") {
		peak_en_copy = peak_en;
		// Find indices for intensity??
	}
	else if (mode == "top") {
		size_t top = (size_t)peak_en[0];
		vecd intensity_copy = intensity;
		indices = ed::top_n(intensity_copy,top);
		peak_en_copy.resize(top);
		for (size_t i = 0; i < top; ++i) peak_en_copy[i] = energy[indices[i]];
	}
	vecd occ_gs = ref_gs ? occupation(GS,gsi,false) : vecd(GS.num_vorb/2,0);
	for (size_t p = 0; p < peak_en_copy.size(); p++) {
		vector<bindex> peaks;
		for (auto &blk : EX.hblks) {	
			for (size_t i = 0; i < blk.size; ++i) {
				if (abs(blk.eig[i]-ref_en-peak_en_copy[p]) < PRINT_TOL) 
					peaks.push_back(bindex(&blk-&EX.hblks[0],i));
			}
		}
		if (mode == "top" && abs(intensity[indices[p]]) < PRINT_TOL) break;
		cout << "PEAK " << p+1;
		if (mode == "top") cout << ", intensity: " << intensity[indices[p]] << endl;
		else cout << endl;
		cout << "absorption energy: " << peak_en_copy[p] << ", degeneracy: " << peaks.size() << endl;
		vecd occ = occupation(EX,peaks,false);
		if (ref_gs) {
			cout << "----------------------------Rel to GS--" << endl;
			vecd occ_rel = occ;
			for (size_t i = 0; i < occ_rel.size(); ++i) occ_rel[i] -= occ_gs[i];
			occ.insert(occ.end(),occ_rel.begin(),occ_rel.end());
		} else cout << "---------------------------------------" << endl;
		EX.cluster->print_eigstate(occ);
		cout << "---------------------------------------" << endl << endl;
	}
	return;
}

void write_XAS(vecd const& aben, vecd const& intensity, string file_dir, bool print) {
	std::ofstream xasfile;
	xasfile.open(file_dir);
	xasfile << setw(15) << "peaks (eV)" << setw(15) << "intensity" << endl;
	if (print) cout << fixed << setprecision(5) << setw(15) << "peaks" << setw(15) << "intensity" << endl;
	for (int i = 0; i < intensity.size(); ++i) {
		if (intensity[i] < TOL) continue;
		xasfile << setw(15) << aben[i] << setw(15) << intensity[i] << endl;
		if (print) {
			if (intensity[i] < PRINT_TOL) continue;
			cout << fixed << setprecision(5) << setw(15) << aben[i] << setw(15) << intensity[i] << endl;
		}
	}
	xasfile.close();
}

void XAS(Hilbert& GS, Hilbert& EX, const PM& pm) {
	// Note: different diagonalize routine might yield different results, 
	// if the width of delta function is not small enough.
	double beta = 0, nedos = pm.nedos;//, racah_B = (SC[2]/49) - (5*SC[4]/441);

	if (GS.hblks[0].eig == nullptr) throw runtime_error("Hamiltonian not diagonalized");
	if (EX.hblks[0].eig == nullptr) throw runtime_error("Hamiltonian not diagonalized");
	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<bindex> gsi; // index for ground state
	// Calculate Partition function
	double Z = 0, emin = pm.ab_range[0], emax = pm.ab_range[1], SDegen = 0;
	for (size_t i = 0; i < GS.hblks.size(); ++i) {
	for (size_t j = 0; j < GS.hblks[i].size; ++j) {
		if (abs(GS.hblks[i].eig[j]-gs_en) < TOL) {
			Z += exp(-beta*GS.hblks[i].eig[j]);
			if (abs(GS.hblks[i].get_sz()) >= SDegen) SDegen = abs(GS.hblks[i].get_sz());
			gsi.push_back({i,j});
		}
	}}

	cout << "Calculating cross section..." << endl;
	auto start = chrono::high_resolution_clock::now();
	vector<blapIndex> blap;

	vecd xas_aben(nedos,0), xas_int(nedos,0);
	for (auto &exblk : EX.hblks) {
		size_t last_gs_block = gsi[0].first;
		basis_overlap(GS,EX,bindex(last_gs_block,&exblk-&EX.hblks[0]),blap,pm);
		for (auto &g  : gsi) {
			Block& gsblk = GS.hblks[g.first];
			ed::print_progress((&exblk-&EX.hblks[0])*gsi.size()+(&g-&gsi[0])+1,EX.hblks.size()*gsi.size());
			if (!GS.SO_on && !EX.SO_on && gsblk.get_sz() != exblk.get_sz()) continue; // Spin order blocks
			if (g.first != last_gs_block) {
				basis_overlap(GS,EX,bindex(g.first,&exblk-&EX.hblks[0]),blap,pm);
				last_gs_block = g.first;
			}
			#pragma omp parallel for reduction (vec_double_plus:xas_int) schedule(dynamic)
			for (size_t ei = 0; ei < exblk.size; ++ei) {
				if (exblk.eig[ei]-gs_en < emin || exblk.eig[ei]-gs_en > emax) continue;
				dcomp cs = 0;
				for (auto & b : blap) {
					size_t gsind = g.second*gsblk.size+b.g;
					size_t exind = ei*exblk.size+b.e;
					cs +=  gsblk.eigvec[gsind] * exblk.eigvec[exind] * b.blap;
				}
				if (abs(cs) > TOL) { // Is this arbritrary?????	
					size_t peak_pos = round((exblk.eig[ei]-gs_en-emin)/((emax-emin)/nedos));
					xas_aben[peak_pos] = exblk.eig[ei]-gs_en;
					// Peaks position might be different if delta function is too wide, we can fix this by taking the smaller value
					xas_int[peak_pos] += exp(-beta*gs_en)*pow(abs(cs),2);
				}
			}
		}
	}
	// Calculate peak intensity, TODO: output peaks?
	XAS_peak_occupation(GS,EX,vecd({10}),xas_aben,xas_int,gsi,gs_en,"top",true);

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";
	cout << "Writing results..." << endl << endl;
	write_XAS(xas_aben,xas_int,"XAS_"+pm.edge+"edge_"+pol_str(pm.pvin)+".txt");
	return;
}


void RIXS_peak_occupation(Hilbert& GS, Hilbert& EX, vecd const& peak_en, vecd const& ab_en, 
						vecd const& em_en, vecd const& intensity, vector<bindex> const& gsi,
						PM const& pm, double ref_en, string mode, bool ref_gs) {
	// Returns top n peaks
	vecd peak_ab_en, peak_em_en;
	vector<size_t> indices;
	if (mode == "list") {
		// input should be in order of: ab/em/ab/em...
		for (size_t i = 0; i < peak_en.size(); ++i) {
			peak_ab_en.push_back(peak_en[i]);
			peak_em_en.push_back(peak_en[i+1]);
			// Find indices of intensity??
		}
	} else if (mode == "top") {
		size_t top = (size_t)peak_en[0];
		peak_ab_en.resize(top);
		peak_em_en.resize(top);
		vecd intensity_copy = intensity;
		indices = ed::top_n(intensity_copy,top);
		for (size_t i = 0; i < top; ++i) {
			size_t abind = indices[i] % pm.nedos;
			size_t emind = indices[i] / pm.nedos;
			peak_ab_en[i] = ab_en[abind];
			peak_em_en[i] = em_en[emind];
		}
	}
	vecd occ_gs = ref_gs ? occupation(GS,gsi,false) : vecd(GS.num_vorb/2,0);
	for (size_t p = 0; p < peak_ab_en.size(); p++) {
		vector<bindex> ab_peaks, em_peaks;
		// Get absorption peaks
		for (auto &blk : EX.hblks) {	
			for (size_t i = 0; i < blk.size; ++i) {
				if (abs(blk.eig[i]-ref_en-peak_ab_en[p]) < PRINT_TOL) 
					ab_peaks.push_back(bindex(&blk-&EX.hblks[0],i));
			}
		}
		// Get emission peaks
		double loss_en = pm.eloss ? peak_em_en[p] : peak_ab_en[p]-peak_em_en[p];
		for (auto &blk : GS.hblks) {	
			for (size_t i = 0; i < blk.size; ++i) {
				if (abs(blk.eig[i]-ref_en-loss_en) < PRINT_TOL) 
					em_peaks.push_back(bindex(&blk-&GS.hblks[0],i));
			}
		}
		if (mode == "top" && abs(intensity[indices[p]]) < PRINT_TOL) break;
		cout << "PEAK " << p+1;
		if (mode == "top") cout << ", intensity: " << intensity[indices[p]] << endl;
		else cout << endl;
		cout << "absorption energy: " << peak_ab_en[p] << ", degeneracy: " << ab_peaks.size() << endl;
		cout << "emission energy: " << peak_em_en[p] << ", degeneracy: " << em_peaks.size() << endl;
		vecd occ_ex = occupation(EX,ab_peaks,false);
		vecd occ_fs = occupation(GS,em_peaks,false);
		occ_ex.insert(occ_ex.end(),occ_fs.begin(),occ_fs.end());
		cout << "--------------------Abs----------Em----" << endl;
		EX.cluster->print_eigstate(occ_ex);
		if (ref_gs) {
			cout << "--------Relative to Ground state-------" << endl;
			for (size_t i = 0; i < occ_ex.size()/2; ++i) {
				occ_ex[i] -= occ_gs[i];
				occ_ex[i+11] -= occ_gs[i];
			}
			EX.cluster->print_eigstate(occ_ex);
		}
		cout << "---------------------------------------" << endl << endl;
	}
	return;
}

void write_RIXS(vecd const& peaks, vecd const& ab, vecd const& em, string file_dir, bool eloss, bool print) {
	std::ofstream rixsfile;
	rixsfile.open(file_dir);
	if (eloss) rixsfile << setw(18) << "absorption (eV)" << setw(18) << "energy loss (eV)" << setw(18) << "intensity" << endl;
	else rixsfile << setw(18) << "absorption (eV)" << setw(18) << "emission (eV)" << setw(18) << "intensity" << endl;
	if (print) {
		if (eloss) cout << setw(18) << "absorption (eV)" << setw(18) << "energy loss (eV)" << setw(18) << "intensity" << endl;
		else cout << setw(18) << "absorption (eV)" << setw(18) << "emission (eV)" << setw(18) << "intensity" << endl;
	}
	for (int y = 0; y < em.size(); ++y) {
		for (int x = 0; x < ab.size(); ++x) {
			if (peaks[y*ab.size()+x] < TOL) continue;
			rixsfile << setw(18) << ab[x];
			rixsfile << setw(18) << em[y];
			rixsfile << setw(18) << peaks[y*ab.size()+x] << endl;
			if (print) {
				if (peaks[y*ab.size()+x] < PRINT_TOL) continue;
				cout << setw(18) << ab[x];
				cout << setw(18) << em[y];
				cout << setw(18) << peaks[y*ab.size()+x] << endl;
			}
	 	}
	}
	rixsfile.close();
	return;
}

void RIXS(Hilbert& GS, Hilbert& EX, const PM& pm) {
	double beta = 0, hbar = 6.58e-16, nedos = pm.nedos;
	// double Racah_B = (SC[2]/49) - (5*SC[4]/441);
	dcomp igamma(0,1*0.1);

	if (GS.hblks[0].eig == nullptr) throw runtime_error("Hamiltonian not diagonalized");
	if (EX.hblks[0].eig == nullptr) throw runtime_error("Hamiltonian not diagonalized");

	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<bindex> gsi,exi; // index for ground state and excited states
	double Z = 0, ab_emin = pm.ab_range[0], ab_emax = pm.ab_range[1], elm_min, elm_max;
	if (pm.eloss) { // Energy Loss
		elm_max = pm.em_energy;
		elm_min = 0;
	} else { // Emission Energy
		elm_max = ab_emax + pm.em_energy;
		elm_min = ab_emin - pm.em_energy;
	}
	for (size_t i = 0; i < GS.hblks.size(); ++i) {
	for (size_t j = 0; j < GS.hblks[i].size; ++j) {
		if (abs(GS.hblks[i].eig[j]-gs_en) < TOL) {
			gsi.push_back({i,j});
			Z += exp(-beta*GS.hblks[i].eig[j]);
		}
	}} 

	// Figure out energy levels for excited states within the spectra range
	vecd exen;
	for (auto &exblk : EX.hblks) {
		exblk.init_einrange();
		for (int ei = 0; ei < exblk.size; ++ei) {
			if (exblk.eig[ei] < gs_en + ab_emin || exblk.eig[ei] > gs_en + ab_emax) {
				exblk.einrange[ei] = -1;
				continue;
			}
			bool is_dup = false;
			for (auto &ex : exen) {
				if (abs(exblk.eig[ei] - ex) < TOL) {
					is_dup = true;
					break;
				}
			}
			if (!is_dup) exen.push_back(exblk.eig[ei]);
		}
	}

	sort(exen.begin(),exen.end());
	if (exen.size() > 1) { // Skip this if there are less than 1 excited state
		for (auto &exblk : EX.hblks) {
		#pragma omp parallel for shared(exblk)
		for (size_t ei = 0; ei < exblk.size; ++ei) {
			if (exblk.einrange[ei] == -1) continue;
			exblk.einrange[ei] = ed::binary_search(exen,exblk.eig[ei],std::greater<double>(),
					[&](double a, double b){return abs(a-b) < TOL;});
		}}
	}

	cout << "Calculating cross section..." << endl;
	auto start = chrono::high_resolution_clock::now();

	// Calculate and store <v|D|i>
	vecd rixs_em(nedos,0),rixs_ab(nedos,0),rixs_loss(nedos,0),rixs_peaks(nedos*nedos,0);
	for (int i = 0; i < nedos; ++i) rixs_ab[i] = ab_emin + i*(ab_emax-ab_emin)/nedos;
	vector<dcomp> rixskern(gsi.size()*EX.hsize,0);
	vector<blapIndex> blap;
	vector<bindex> peak_occ_ab, peak_occ_em; // Calculate peak occupation
	cout << "Calculating <v|D|i>" << endl;
	for (auto &exblk : EX.hblks) {
		size_t last_gs_block = gsi[0].first;
		basis_overlap(GS,EX,bindex(last_gs_block,&exblk-&EX.hblks[0]),blap,pm);
		for (auto &g : gsi) {
			Block& gsblk = GS.hblks[g.first];
			ed::print_progress((&exblk-&EX.hblks[0])*gsi.size()+(&g-&gsi[0])+1,EX.hblks.size()*gsi.size());
			if (!GS.SO_on && !EX.SO_on && gsblk.get_sz() != exblk.get_sz()) continue;
			if (g.first != last_gs_block) {
				basis_overlap(GS,EX,bindex(g.first,&exblk-&EX.hblks[0]),blap,pm);
				last_gs_block = g.first;
			}	
			#pragma omp parallel for reduction (vec_dcomp_plus:rixskern) schedule(dynamic)
			for (size_t ei = 0; ei < exblk.size; ++ei) {
				if (exblk.eig[ei]-gs_en < ab_emin || exblk.eig[ei]-gs_en > ab_emax) continue;
				dcomp csvi = 0;
				for (auto & b : blap) {
					size_t gsind = g.second*gsblk.size+b.g;
					size_t exind = ei*exblk.size+b.e;
					csvi +=  gsblk.eigvec[gsind] * exblk.eigvec[exind] * b.blap;
				}
				if (abs(csvi) > TOL) rixskern[(&g-&gsi[0])*EX.hsize+ei+exblk.f_ind] += csvi;
			}
		}
	}

	// Loop through final states, calculate <f|D|v><v|D|i>, add to spectra
	cout << endl << "Calculating <f|D|v><v|D|i>" << endl;
	for (auto &fsblk : GS.hblks) {
	for (auto &exblk : EX.hblks) {	
		if (!GS.SO_on && !EX.SO_on && fsblk.get_sz() != exblk.get_sz()) continue; // Spin order blocks
		cout << "FS blk: " << &fsblk-&GS.hblks[0] << ", EX blk: " << &exblk-&EX.hblks[0] << endl;
		basis_overlap(GS,EX,bindex(&fsblk-&GS.hblks[0],&exblk-&EX.hblks[0]),blap,pm,true);
		for (size_t fi = 0; fi < fsblk.size; ++fi) {
			ed::print_progress((double)fi+1,(double)fsblk.size);
			if (pm.eloss) {if (fsblk.eig[fi]-gs_en > elm_max) continue;}
			vector<dcomp> csum(gsi.size()*exen.size(),0);
			bool is_gs = abs(fsblk.eig[fi]-gs_en) < TOL;
			size_t gsind;
			if (is_gs) gsind = distance(gsi.begin(), find(gsi.begin(), gsi.end(), bindex(&fsblk-&GS.hblks[0],fi)));
			#pragma omp parallel for reduction (vec_dcomp_plus:csum) schedule(dynamic)
			for (size_t ei = 0; ei < exblk.size; ++ei) {
				if (exblk.einrange[ei] == -1) continue; // Not in absorption range
				if (!pm.eloss) {if (exblk.eig[ei]-fsblk.eig[fi] > elm_max || exblk.eig[ei]-fsblk.eig[fi] < elm_min) continue;}
				// Breaks if no ground state connects with intermediate state
				bool skip_ei = true;
				for (size_t gi = 0; gi < gsi.size(); ++gi) {
					if (abs(rixskern[gi*EX.hsize+ei+exblk.f_ind]) > TOL) {
						skip_ei = false;
						break;
					}
				}
				if (skip_ei) continue;
				// Calculate Cross Section
				dcomp csvf = 0;
				if (is_gs && pm.pvin == pm.pvout) csvf = rixskern[gsind*EX.hsize+ei+exblk.f_ind];
				else {
					for (size_t b = 0; b < blap.size(); ++b) {
						size_t fsind = fi*fsblk.size+blap[b].g;
						size_t exind = ei*exblk.size+blap[b].e;
						csvf +=  fsblk.eigvec[fsind] * exblk.eigvec[exind] * blap[b].blap;
					}
				}
				if (abs(csvf) < TOL) continue;
				for (size_t gi = 0; gi < gsi.size(); ++gi) {
					if (abs(rixskern[gi*EX.hsize+ei+exblk.f_ind]) < TOL) continue;
					csum[gi*exen.size()+exblk.einrange[ei]] += conj(csvf) * rixskern[gi*EX.hsize+ei+exblk.f_ind];
				}
			}
			/* Classical Implementation of Kramers-Heisenberg Equation
	 		// Sweep through the absorption frequency
	 		// for (int g = 0; g < gsi.size(); g++) {
	 		// 	for (int ai = 0; ai < nedos; ++ai) {
	 		// 		dcomp cs(0,0);
				// 	for (int e = 0; e < exen.size(); e++) {
				// 		if (abs(csum[g*exen.size()+e]) < TOL) continue;
		 	// 			cs += csum[g*exen.size()+e]/(rixs_ab[ai]-exen[e]+gs_en-igamma);
				// 	}
				// 	int emind = round((rixs_ab[ai]-fsblk.eig[fi]+gs_en-em_emin)/(em_emax-em_emin)*nedos);
				// 	rixs_em[emind] = rixs_ab[ai]-fsblk.eig[fi]+gs_en;
				// 	int lossind = round((fsblk.eig[fi]-gs_en-el_min)/(el_max-el_min)*nedos);
				// 	rixs_loss[lossind] = fsblk.eig[fi]-gs_en;
				// 	if (fsblk.eig[fi]-gs_en > el_max) continue;
				// 	rixs_peaks[lossind+ai*nedos] += exp(-beta*gs_en)*pow(abs(cs),2);
				// }
	 		// }
	 		*/
			for (size_t g = 0; g < gsi.size(); g++) {
				for (size_t e = 0; e < exen.size(); e++) {
					dcomp cs(0,0);
					if (abs(csum[g*exen.size()+e]) < TOL) continue;
					int abind = round((exen[e]-gs_en-ab_emin)/(ab_emax-ab_emin)*nedos);
					int elmind;
					rixs_ab[abind] = exen[e]-gs_en;
					if (pm.eloss) {
						elmind = round((fsblk.eig[fi]-gs_en)/elm_max*nedos);
						rixs_em[elmind] = fsblk.eig[fi]-gs_en;
					} else {
						elmind = round((exen[e]-fsblk.eig[fi]-elm_min)/(elm_max-elm_min)*nedos);
						rixs_em[elmind] = exen[e]-fsblk.eig[fi];
					}
					rixs_peaks[elmind*nedos+abind] += exp(-beta*gs_en)*pow(abs(csum[g*exen.size()+e]),2);
				}
			}
		}
	}}
	RIXS_peak_occupation(GS,EX,vecd({6}),rixs_ab,rixs_em,rixs_peaks,gsi,pm,gs_en,"top",true);

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";
	cout << "Writing results..." << endl << endl;
	write_RIXS(rixs_peaks,rixs_ab,rixs_em,"RIXS_"+pm.edge+"edge_"+pol_str(pm.pvin)+"_"+pol_str(pm.pvout)+".txt",pm.eloss,true);
	return;
}