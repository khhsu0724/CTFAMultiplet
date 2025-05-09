#include <complex>
#include <bitset>
#include "multiplet.hpp"
#include "photon.hpp"

// This file handles calculating cross section/Occupation/CI configuration

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

vecd occupation(Hilbert& hilbs, const vector<bindex>& si, bool is_print, string fname, bool spin_res) {
	// Calculation occupation of orbitals, only valid when matrix is diagonalized
	// There is a bug when calculating octahedral cluster????
	for (auto& blk : hilbs.hblks) if (blk.eigvec == nullptr) throw runtime_error("matrix not diagonalized for occupation");
	int nvo = hilbs.cluster->vo_persite, nco = hilbs.cluster->co_persite;
	vecc U = hilbs.cluster->get_seph2real_mat();
	vector<double> occ_totsu(nvo,0), occ_totsd(nvo,0);
	Hilbert minus_1vh(hilbs,-1);
	int gs_count = 0;
	for (auto &s  : si) {
		auto& blk = hilbs.hblks[s.first];
		if (spin_res && s.first >= (hilbs.hblks.size()+1)/2) continue;
		gs_count++;
		// Pick an operator
		for (size_t c = 0; c < nvo; ++c) {
			vector<dcomp> wvfncsu(minus_1vh.hsize,0);
			vector<dcomp> wvfncsd(minus_1vh.hsize,0);
			for (size_t ci = 0; ci < nvo; ++ci) {
				if (abs(U[c*nvo+ci]) < TOL) continue;
				for (size_t j = 0; j < blk.size; ++j) {
					double sj = blk.eigvec[s.second*blk.size+j];
					if (abs(sj) < TOL) continue;
					// Operate on spin up and spin down
					ulli op_state = hilbs.Hashback(bindex(s.first,j));
					ulli opsu = BIG1 << (ci+nco), opsd = BIG1 << (ci+2*nco+nvo);
					if (opsu & op_state) {
						bindex ind = minus_1vh.Hash(op_state-opsu);
						int fsgn = hilbs.Fsign(&opsu,op_state,1);
						wvfncsu[ind.second] += fsgn*sj*U[c*nvo+ci];
					}
					if (opsd & op_state) {
						bindex ind = minus_1vh.Hash(op_state-opsd);
						int fsgn = hilbs.Fsign(&opsd,op_state,1);
						wvfncsd[ind.second] += fsgn*sj*U[c*nvo+ci];
					}
				}
			}
			for (auto &w : wvfncsu) occ_totsu[c] += abs(pow(w,2));
			for (auto &w : wvfncsd) occ_totsd[c] += abs(pow(w,2));
		}
	}
	vector<double> occ_tot(nvo,0);
	for (int i = 0; i < nvo; ++i) {
		occ_totsu[i] = occ_totsu[i]/gs_count;
		occ_totsd[i] = occ_totsd[i]/gs_count;
	}
	for (int i = 0; i < nvo; ++i) occ_tot[i] = occ_totsd[i] + occ_totsu[i];
	if (spin_res) {
		// cout << "Spin up: " << endl;
		hilbs.cluster->print_eigstate(occ_totsu,is_print,fname);
		// cout << "Spin down: " << endl;
		hilbs.cluster->print_eigstate(occ_totsd,is_print,fname);
	} else {
		hilbs.cluster->print_eigstate(occ_tot,is_print,fname);
	}
	return occ_tot;
}

vector<double> wvfnc_weight(Hilbert& hilbs, const vector<bindex>& si, int ligNum, bool print) {
	// Returns fraction of degenerate states with their ligand occupation
	ligNum = (ligNum < hilbs.num_vh) ? ligNum : hilbs.num_vh;
	vector<double> wvfnc(hilbs.hsize), dLweight(ligNum+1,0);
	for (auto& s : si) {
		auto &blk = hilbs.hblks[s.first];
		for (size_t i = 0; i < blk.size; ++i) {
			// Calculate Wavefunction
			if (abs(blk.eigvec[s.second*blk.size+i]) < TOL) continue;
			wvfnc[blk.f_ind+i] += pow(blk.eigvec[s.second*blk.size+i],2);
			ulli state = hilbs.Hashback(bindex(s.first,i));
			int Lcnt = 0;
			// This needs to account geometry, Need a cleaner solution?
			for (size_t j = 5*hilbs.cluster->tm_per_site; j < hilbs.cluster->vo_persite; ++j) {
				if ((BIG1<<(j+hilbs.num_corb/2)) & state) Lcnt++;
				if ((BIG1<<(j+hilbs.num_corb+hilbs.num_vorb/2)) & state) Lcnt++;
			}
			if (Lcnt <= ligNum) dLweight[Lcnt] += pow(blk.eigvec[s.second*blk.size+i],2)/si.size();
		}
	}
	if (print) {
		cout << "Ground State composition";
		for (size_t i = 0; i <= ligNum; ++i) {
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
	ligNum = (ligNum < hilbs.num_vh) ? ligNum : hilbs.num_vh;
	vector<bool> firstLh(ligNum+1,true);
	auto all_eig = hilbs.get_all_eigval(false);
	vector<double> distinct = ed::printDistinct(all_eig,0.0,all_eig.size(),is_print);
	for (size_t e = 0; e < distinct.size(); ++e) {
		vector<bindex> gei;
		for (size_t i = 0; i < hilbs.hblks.size(); ++i) {
			if (hilbs.hblks[i].eig == nullptr) continue;
			for (size_t j = 0; j < hilbs.hblks[i].nev; ++j) {
				if (abs(hilbs.hblks[i].eig[j]-distinct[e]) < TOL) {
					gei.push_back({i,j});
				}
			}
		}

		vector<double> dLweight = wvfnc_weight(hilbs,gei,ligNum,false);
		for (size_t i = 0; i <= ligNum; ++i) {
			if (firstLh[i] && dLweight[i] > TOL) { //Cases where different dn states mix
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
	if (firstLh[1]) cout << "Warning: dL energy not set, check effective_delta" << endl;
	if (is_print) cout << "Effective delta: " << dL - d << ", d energy: " << d << ", dL energy: " << dL << endl;
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
	const int gsblk_size = GS.hblks[gbi].size;
	const int exblk_size = EX.hblks[exi].size;
	vector<ulli> gslist = GS.get_hashback_list(gbi);
	vector<ulli> exslist = EX.get_hashback_list(exi);
	#pragma omp parallel
	{
	    std::vector<blapIndex> local_blap;
	    local_blap.reserve(EX.hblks[exi].size*2/omp_get_max_threads());
	    #pragma omp for collapse(2) nowait
	    for (size_t g = 0; g < gsblk_size; g++) {
	        for (size_t e = 0; e < exblk_size; e++) {
	            ulli gs = gslist[g], exs = exslist[e];
	            ulli vh = gs - (gs & exs);
	            if (!ed::is_pw2(vh)) continue;
	            ulli ch = exs - (gs & exs);
	            int coi = EX.orbind(ch);
	            if (coi == -1) continue;
	            int voi = GS.atlist[coi].vind;
	            if (!GS.atlist[voi].contains(vh)) continue;
	            QN chqn = EX.atlist[coi].fast_qn(ch,half_orb,coi);
	            QN vhqn = GS.atlist[voi].fast_qn(vh,half_orb,voi);
	            if (chqn.spin != vhqn.spin || abs(vhqn.ml-chqn.ml) > 1) continue;
	            dcomp blap_val = gaunt(cl,chqn.ml,vl,vhqn.ml)[1] * pow(-1,vhqn.ml-chqn.ml+1)
	                            * GS.Fsign(&vh,gs,1) * EX.Fsign(&ch,exs,1) * proj_pvec(vhqn.ml-chqn.ml,pvec);
	            if (std::norm(blap_val) > 1E-14) {
	                local_blap.emplace_back(g, e, blap_val);
	            }
	        }
	    }

	    // Merge local results back to global vector
	    #pragma omp critical
	    {
	        blap.insert(blap.end(), local_blap.begin(), local_blap.end());
	    }
	}
	return;
}

vecc gen_dipole_state(Hilbert& GS, Hilbert& EX, const PM& pm, const bindex& inds, vecc vec_in, 
						const vector<blapIndex>& blap, bool excite) {
	size_t gbi = inds.first, exi = inds.second;
	if (excite) {
		// Generate excited state with D|g>
		if (vec_in.size() != GS.hblks[gbi].size) {
			cout << "Block size mismatch in gen_dipole_state" << endl;
			exit(1);
		}
		vecc dpvec(EX.hblks[exi].size,0);
		for (auto & b : blap) {
			dpvec[b.e] += vec_in[b.g] * b.blap;
		}
		return dpvec;
	} else { // Generate excited state with D+|g>
		if (vec_in.size() != EX.hblks[exi].size) {
			cout << "Block size mismatch in gen_dipole_state (rev)" << endl;
			exit(1);
		}
		vecc dpvec(GS.hblks[gbi].size,0);
		for (auto & b : blap) {
			dpvec[b.g] += vec_in[b.e] * conj(b.blap);
		}
		return dpvec;
	}
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
			for (size_t i = 0; i < blk.nev; ++i) {
				if (abs(blk.eig[i]-ref_en-peak_en_copy[p]) < 1e-8) 
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

void write_XAS(vecd const& aben, vecd const& intensity, string file_dir, bool exact) {
	multistream mout(true,file_dir,"w");
	mout << setw(15) << "peaks (eV)" << setw(15) << "intensity" << endl;
	for (int i = 0; i < intensity.size(); ++i) {
		if (intensity[i] < TOL and exact) continue;
		mout.set_skip_cout(intensity[i] < PRINT_TOL or !exact);
		mout << setw(15) << aben[i] << setw(15) << intensity[i] << endl;
	}
	return;
}


void XAS(Hilbert& GS, Hilbert& EX, const PM& pm) {
	// Note: different diagonalize routine might yield different results, 
	// if the width of delta function is not small enough.
	double beta = 0, nedos = pm.nedos;//, racah_B = (SC[2]/49) - (5*SC[4]/441);

	if (GS.hblks[0].eig == nullptr) throw runtime_error("Hamiltonian not diagonalized");
	double gs_en = GS.hblks[0].eig[0], ex_en;
	if (!pm.skip_ch_diag) {
		if (EX.hblks[0].eig == nullptr) throw runtime_error("Hamiltonian not diagonalized");
		ex_en = EX.hblks[0].eig[0];
	}
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.nev; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<bindex> gsi; // index for ground state
	// Calculate Partition function
	double Z = 0, emin = pm.ab_range[0], emax = pm.ab_range[1], SDegen = 0;
	for (size_t i = 0; i < GS.hblks.size(); ++i) {
	for (size_t j = 0; j < GS.hblks[i].nev; ++j) {
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
	bool write_output = true;
	if (pm.spec_solver == 4) {
		write_output = false;
		// Lanczos solver
		for (int i = 0; i < nedos; ++i) 
			xas_aben[i] = pm.ab_range[0] + (pm.ab_range[1]-pm.ab_range[0])/nedos*i;
		for (auto &exblk : EX.hblks) {
			size_t last_gs_block = gsi[0].first;
			size_t exblk_ind = &exblk-&EX.hblks[0];
			basis_overlap(GS,EX,bindex(last_gs_block,&exblk-&EX.hblks[0]),blap,pm);
			for (auto &g  : gsi) {
				auto& gsblk = GS.hblks[g.first];
				size_t gblk_size =  GS.hblks[g.first].size;
				// If spin ordered block, check for no spin flip
				if (!GS.SO_on && !EX.SO_on && GS.hblks[g.first].get_sz() != exblk.get_sz()) continue;
				if (g.first != last_gs_block) {
					basis_overlap(GS,EX,bindex(g.first,&exblk-&EX.hblks[0]),blap,pm);
					last_gs_block = g.first;
				}
				// cout << "gsblk: " << g.first << ", exblk: " << &exblk-&EX.hblks[0] << endl;
				vecc gs_vec(gblk_size,0);
				#pragma omp parallel for
				for (int i = 0; i < gblk_size; ++i) 
					gs_vec[i] = dcomp(GS.hblks[g.first].eigvec[g.second*gblk_size+i],0);		
				if (g.first != last_gs_block) {
					basis_overlap(GS,EX,bindex(g.first,exblk_ind),blap,pm);
					last_gs_block = g.first;
				}
				vecc dipole_vec = gen_dipole_state(GS,EX,pm,bindex(g.first,exblk_ind),gs_vec,blap);
				// Perform Lanczos
				ContFracExpan(exblk.ham,dipole_vec,gs_en,xas_aben,xas_int,pm.eps_ab,pm.niterCFE);
			}
		}
		cout << "Finish solving Lanczos!" << endl;
	} else {
		for (auto &exblk : EX.hblks) {
			size_t last_gs_block = gsi[0].first;
			basis_overlap(GS,EX,bindex(last_gs_block,&exblk-&EX.hblks[0]),blap,pm);
			for (auto &g  : gsi) {
				auto& gsblk = GS.hblks[g.first];
				// ed::print_progress((&exblk-&EX.hblks[0])*gsi.size()+(&g-&gsi[0])+1,EX.hblks.size()*gsi.size());
				if (!GS.SO_on && !EX.SO_on && gsblk.get_sz() != exblk.get_sz()) continue; // Spin order blocks
				if (g.first != last_gs_block) {
					basis_overlap(GS,EX,bindex(g.first,&exblk-&EX.hblks[0]),blap,pm);
					last_gs_block = g.first;
				}
				#pragma omp parallel for reduction (vec_double_plus:xas_int) schedule(dynamic)
				for (size_t ei = 0; ei < exblk.nev; ++ei) {
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
		XAS_peak_occupation(GS,EX,vecd({10}),xas_aben,xas_int,gsi,gs_en,"top",write_output);
	}

	// Normalize intensity
	for (auto & val : xas_int) val /= gsi.size();
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";
	cout << "Writing results..." << endl << endl;
	write_XAS(xas_aben,xas_int,"XAS_"+pm.edge+"edge_"+pol_str(pm.pvin)+".txt",write_output);
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
			for (size_t i = 0; i < blk.nev; ++i) {
				if (abs(blk.eig[i]-ref_en-peak_ab_en[p]) < TOL) 
					ab_peaks.push_back(bindex(&blk-&EX.hblks[0],i));
			}
		}
		// Get emission peaks
		double loss_en = pm.eloss ? peak_em_en[p] : peak_ab_en[p]-peak_em_en[p];
		double tolerance = pm.eloss ? TOL : (2*pm.em_energy+pm.ab_range[1]-pm.ab_range[0])/ab_en.size();
		for (auto &blk : GS.hblks) {	
			for (size_t i = 0; i < blk.nev; ++i) {
				if (abs(blk.eig[i]-ref_en-loss_en) < tolerance) 
					em_peaks.push_back(bindex(&blk-&GS.hblks[0],i));
			}
		}
		if (mode == "top" && abs(intensity[indices[p]]) < PRINT_TOL) break;
		cout << "PEAK " << p+1;
		if (mode == "top") cout << ", intensity: " << intensity[indices[p]] << endl;
		else cout << endl;
		cout << "absorption energy: " << peak_ab_en[p] << ", degeneracy: " << ab_peaks.size() << endl;
		if (pm.eloss) cout << "loss energy: ";
		else cout << "emission energy: ";
		cout << peak_em_en[p] << ", degeneracy: " << em_peaks.size() << endl;
		vecd occ_ex = occupation(EX,ab_peaks,false);
		vecd occ_fs = occupation(GS,em_peaks,false);
		occ_ex.insert(occ_ex.end(),occ_fs.begin(),occ_fs.end());
		cout << "--------------------Abs----------Em----" << endl;
		EX.cluster->print_eigstate(occ_ex);
		if (ref_gs) {
			cout << "--------Relative to Ground state-------" << endl;
			for (size_t i = 0; i < occ_ex.size()/2; ++i) {
				occ_ex[i] -= occ_gs[i];
				occ_ex[i+occ_ex.size()/2] -= occ_gs[i];
			}
			EX.cluster->print_eigstate(occ_ex);
		}
		cout << "---------------------------------------" << endl << endl;
	}
	return;
}

void write_RIXS(vecd const& peaks, vecd const& ab, vecd const& em, bool eloss, string file_dir) {
	multistream mout(true,file_dir,"w");
	if (eloss) mout << setw(18) << "absorption (eV)" << setw(18) << "energy loss (eV)" << setw(18) << "intensity" << endl;
	else mout << setw(18) << "absorption (eV)" << setw(18) << "emission (eV)" << setw(18) << "intensity" << endl;
	for (int y = 0; y < em.size(); ++y) {
		for (int x = 0; x < ab.size(); ++x) {
			if (peaks[y*ab.size()+x] < TOL) continue;
			mout.set_skip_cout(peaks[y*ab.size()+x] < PRINT_TOL);
			mout << setw(18) << ab[x];
			mout << setw(18) << em[y];
			mout << setw(18) << peaks[y*ab.size()+x] << endl;
	 	}
	}
	return;
}

void write_kh_RIXS(vecd const& rixsmat, vecd const& eloss, string file_dir) {
	multistream mout(false,file_dir,"w");
	int nedos = eloss.size();
	mout << setw(6) << "energy loss (eV)" << setw(18) << "absorption array" << endl;
	for (int x = 0; x < nedos; ++x) {
	    bool all_zero = std::all_of(rixsmat.begin()+x*nedos, 
	    				rixsmat.begin()+(x+1)*nedos, [](bool elem){return elem < TOL;});
	    if (all_zero) continue;
		mout << eloss[x] << setw(6) << " ";
		for (int y = 0; y < nedos; ++y) {
			mout << rixsmat[x*nedos+y] << " ";
		}
		mout << endl;
	}
	return;
}

void write_iter_RIXS(vecd const& rixs_ab, vecd const& rixs_loss, vecd const& rixs_int, 
						string file_dir, bool write_init) {
	// Write RIXS in absorption, eloss, intensity basis
	int nedos = rixs_ab.size();
	multistream mout;
	if (write_init) {
		mout = multistream(false,file_dir,"w");
		mout << setw(18) << "absorption energy" << setw(18) << "energy loss";
		mout << setw(18) << "intensity" << endl;
	} else mout = multistream(false,file_dir,"a");
	for (int x = 0; x < nedos; ++x) {
		mout << rixs_ab[x] << setw(18) << " ";
		mout << rixs_loss[x] << setw(18) << " ";
		mout << rixs_int[x] << setw(18) << " ";
		mout << endl;
	}
	return;
}

inline void calc_rixs_blap(Hilbert& GS, Hilbert& EX, size_t gsind, size_t exind, 
							vector<blapIndex>& blap_in, vector<blapIndex>& blap_out, const PM& pm) {
	basis_overlap(GS,EX,bindex(gsind,exind),blap_in,pm);
	if (pm.pvin != pm.pvout) {
		basis_overlap(GS,EX,bindex(gsind,exind),blap_out,pm,true);
	} else blap_out = blap_in;
	return;
}

void RIXS(Hilbert& GS, Hilbert& EX, const PM& pm) {
	double beta = 0, hbar = 6.58e-16, nedos = pm.nedos, eloss_min = -2;
	dcomp igamma(0,pm.eps_loss);

	if (GS.hblks[0].eig == nullptr) throw runtime_error("Hamiltonian not diagonalized");
	if (!pm.skip_ch_diag and EX.hblks[0].eig == nullptr) throw runtime_error("Hamiltonian not diagonalized");

	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.nev; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
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
	for (size_t j = 0; j < GS.hblks[i].nev; ++j) {
		if (abs(GS.hblks[i].eig[j]-gs_en) < TOL) {
			gsi.push_back({i,j});
			Z += exp(-beta*GS.hblks[i].eig[j]);
		}
	}} 

	cout << "Calculating cross section..." << endl;
	auto start = chrono::high_resolution_clock::now();
	vecd rixs_em, rixs_ab, rixs_peaks, rixs_peaks_kh, rixs_em_kh;
	vector<blapIndex> blap, blap_in, blap_out;
	bool reuse_blap = false;

	if (pm.spec_solver == 4) {
		// BiCGstab and Lanczos to solve RIXS spectra
		// Solve for each incident energy, not super efficient
		vecc solved_vec;
		// Special case where there are no symmetry, we can use the same basis overlap.
		if (GS.hblks.size()==1 and EX.hblks.size()==1) {
			reuse_blap = true;
			calc_rixs_blap(GS,EX,0,0,blap_in,blap_out,pm);
		}
		if (pm.precond != 0) solved_vec = vecc(EX.hsize,0);
		for (auto &ab_en: pm.inc_e_points) {
			cout << "---------------Solving for incident energy: " << ab_en << "---------------" << endl;
			dcomp z(gs_en+ab_en,-pm.eps_ab);
			vecd rixs_ab_local(nedos,ab_en);
			vecd rixs_em_local(nedos,0);
			vecd rixs_peaks_local(nedos,0);
			// Using previous coverged vector as guess, use a vector with size = Hilbert space
			for (int i = 0; i < nedos; ++i) rixs_em_local[i] = -2 + (pm.em_energy+2)/nedos*i;
			for (auto &exblk : EX.hblks) {
				size_t last_gs_block = gsi[0].first;
				size_t exblk_ind = &exblk-&EX.hblks[0];
				if (!reuse_blap) calc_rixs_blap(GS,EX,last_gs_block,exblk_ind,blap_in,blap_out,pm);
				for (auto &g  : gsi) {
					auto& gsblk = GS.hblks[g.first];
					size_t gblk_size =  GS.hblks[g.first].size;
					// If spin ordered block, check for no spin flip
					if (!GS.SO_on && !EX.SO_on && GS.hblks[g.first].get_sz() != exblk.get_sz()) continue;
					cout << "Block: " << g.first << ", " << exblk_ind << endl;
					vecc gs_vec(gblk_size,0);
					#pragma omp parallel for
					for (int i = 0; i < gblk_size; ++i) 
						gs_vec[i] = dcomp(GS.hblks[g.first].eigvec[g.second*gblk_size+i],0);
					// Recalculate basis state overlap if different blocks
					if (g.first != last_gs_block && !reuse_blap) {
						calc_rixs_blap(GS,EX,g.first,exblk_ind,blap_in,blap_out,pm);
						last_gs_block = g.first;
					}
					vecc dipole_vec = gen_dipole_state(GS,EX,pm,bindex(g.first,exblk_ind),gs_vec,blap_in);
					// Perform BiCGS, solve for intermediate state
					vecc guess_vec; // This is currently unstable...?
					if (pm.precond != 0 && ab_en != pm.inc_e_points[0]) {
						guess_vec = vecc(exblk.size,0);
						std::copy(solved_vec.begin()+exblk.f_ind, 
							solved_vec.begin()+exblk.f_ind+exblk.size, guess_vec.begin());
					}
					vecc midvec = BiCGS(exblk.ham,dipole_vec,z,pm.CG_tol,guess_vec);
					if (pm.precond != 0) std::copy(midvec.begin(), midvec.end(), solved_vec.begin()+exblk.f_ind);
					// De-excitation
					midvec = gen_dipole_state(GS,EX,pm,bindex(g.first,exblk_ind),midvec,blap_out,false);
					int niter_CFE_in = pm.niterCFE;
					if (niter_CFE_in > GS.hblks[g.first].size/100) niter_CFE_in = GS.hblks[g.first].size/100;
					if (niter_CFE_in < 20) niter_CFE_in = 20; 
					cout << "Number of Lanczos Iteration: " << niter_CFE_in << endl;
					ContFracExpan(GS.hblks[g.first].ham,midvec,gs_en,rixs_em_local,rixs_peaks_local,
									pm.eps_loss,niter_CFE_in);
				}
			}
			bool write_init = (ab_en == pm.inc_e_points[0]);
			// Write per absorption to save progress
			write_iter_RIXS(rixs_ab_local,rixs_em_local,rixs_peaks_local,
				"RIXS_"+pm.edge+"edge_"+pol_str(pm.pvin)+"_"+pol_str(pm.pvout)+".txt",write_init);
			cout << "---------------Done---------------" << endl;
		}
	} else {
		// Figure out energy levels for excited states within the spectra range
		rixs_em = vecd(nedos,0);
		rixs_ab = vecd(nedos,0);
		rixs_peaks = vecd(nedos*nedos,0);
		vecd exen;
		for (auto &exblk : EX.hblks) {
			exblk.init_einrange();
			for (int ei = 0; ei < exblk.nev; ++ei) {
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
			for (size_t ei = 0; ei < exblk.nev; ++ei) {
				if (exblk.einrange[ei] == -1) continue;
				exblk.einrange[ei] = ed::binary_search(exen,exblk.eig[ei],std::greater<double>(),
						[&](double a, double b){return abs(a-b) < TOL;});
			}}
		}
		// Calculate K-H frequency range
		int n_min = round((exen[0]-gs_en-2-ab_emin)/(ab_emax-ab_emin)*nedos);
		int n_max = round((exen[exen.size()-1]-gs_en+2-ab_emin)/(ab_emax-ab_emin)*nedos);
		if (n_min < 0) n_min = 0;
		if (n_max > nedos) n_max = nedos;

		// Calculate and store <v|D|i>
		rixs_peaks_kh = vecd(nedos*nedos,0);
		rixs_em_kh = vecd(nedos,0);
		for (int i = 0; i < nedos; ++i) rixs_ab[i] = ab_emin + i*(ab_emax-ab_emin)/nedos;
		vector<dcomp> rixskern(gsi.size()*EX.hsize,0); // This is very memory exhausive???
		vector<bindex> peak_occ_ab, peak_occ_em; // Calculate peak occupation
		cout << "Calculating <v|D|i>" << endl;
		for (auto &exblk : EX.hblks) {
			size_t last_gs_block = gsi[0].first;
			basis_overlap(GS,EX,bindex(last_gs_block,&exblk-&EX.hblks[0]),blap,pm);
			for (auto &g : gsi) {
				auto& gsblk = GS.hblks[g.first];
				// ed::print_progress((&exblk-&EX.hblks[0])*gsi.size()+(&g-&gsi[0])+1,EX.hblks.size()*gsi.size());
				if (!GS.SO_on && !EX.SO_on && gsblk.get_sz() != exblk.get_sz()) continue;
				if (g.first != last_gs_block) {
					basis_overlap(GS,EX,bindex(g.first,&exblk-&EX.hblks[0]),blap,pm);
					last_gs_block = g.first;
				}	
				int gs_num = (&g-&gsi[0]);
				#pragma omp parallel for shared(rixskern) schedule(dynamic)
				for (size_t ei = 0; ei < exblk.nev; ++ei) {
					if (exblk.eig[ei]-gs_en < ab_emin || exblk.eig[ei]-gs_en > ab_emax) continue;
					dcomp csvi = 0;
					for (auto & b : blap) {
						size_t gsind = g.second*gsblk.size+b.g;
						size_t exind = ei*exblk.size+b.e;
						csvi +=  gsblk.eigvec[gsind] * exblk.eigvec[exind] * b.blap;
					}
					if (abs(csvi) > TOL) rixskern[gs_num*EX.hsize+ei+exblk.f_ind] += csvi;
				}
				// Conjugate the <v|D|i> first
				#pragma omp parallel for shared(rixskern) schedule(dynamic)
				for (size_t r = 0; r < rixskern.size(); r++) {
					rixskern[r] = conj(rixskern[r]);
				}
			}
		}

		// NEW IMPLEMENTATION
		// Loop through final states, calculate <f|D|v><v|D|i>, add to spectra
		cout << endl << "Calculating <f|D|v><v|D|i>" << endl;
		double gamma_tol = 6*pm.eps_loss;
		double freq_step = (ab_emax-ab_emin)/nedos;
		for (auto &fsblk : GS.hblks) {
		for (auto &exblk : EX.hblks) {
			if (!GS.SO_on && !EX.SO_on && fsblk.get_sz() != exblk.get_sz()) continue; // Spin order blocks
			cout << "FS blk: " << &fsblk-&GS.hblks[0] << ", EX blk: " << &exblk-&EX.hblks[0] << endl;
			basis_overlap(GS,EX,bindex(&fsblk-&GS.hblks[0],&exblk-&EX.hblks[0]),blap,pm,true);
			for (size_t fi = 0; fi < fsblk.nev; ++fi) {
				// ed::print_progress((double)fi+1,(double)fsblk.nev);
				vector<dcomp> fDv = vector<dcomp>(exblk.nev,0);
				double fs_en = fsblk.eig[fi];
				if (fs_en-gs_en > pm.em_energy) continue;
				// Calculate <f|D|v> for all v		
				#pragma omp parallel for shared(fDv) schedule(dynamic)
				for (size_t ei = 0; ei < exblk.nev; ++ei) {
					if (exblk.einrange[ei] == -1) continue;
					dcomp csvf = 0;
					for (size_t b = 0; b < blap.size(); ++b) {
						size_t fsind = fi*fsblk.size+blap[b].g;
						size_t exind = ei*exblk.size+blap[b].e;
						csvf +=  fsblk.eigvec[fsind] * exblk.eigvec[exind] * blap[b].blap;
					}
					fDv[ei] = csvf;
				}
				// Calculate sum <f|D|v><v|D|i> for a pair of f,i
				for (auto& g : gsi) {
					int gs_num = &g-&gsi[0];
					// precompute <f|D|v><v|D|i> for all v		
					vector<dcomp> fDvvDi = vector<dcomp>(exblk.nev,0);
					#pragma omp parallel for shared(fDvvDi,fDv,rixskern) schedule(dynamic)
					for (size_t ei = 0; ei < exblk.nev; ++ei) {
						// RIXSKERN already conjugated
						fDvvDi[ei] += fDv[ei] * rixskern[gs_num*EX.hsize+ei+exblk.f_ind];
					}
					// Sweep through absorption frequency
					if (pm.spec_solver == 2 || pm.spec_solver == 3) {
						int eloss_ind = floor((fs_en-gs_en-eloss_min)/(pm.em_energy-eloss_min)*nedos);
						rixs_em_kh[eloss_ind] = fs_en-gs_en;
						#pragma omp parallel for shared(rixs_peaks_kh)
						for (size_t n = n_min; n < n_max; ++n) {
							dcomp intensity = 0;
							double omega_in = ab_emin + n*freq_step;
							// #pragma omp parallel for reduction(+:intensity)
							for (size_t ei = 0; ei < exblk.nev; ++ei) {
								if (abs(fDvvDi[ei]) < TOL) continue;
								// if (abs(omega_in-(exblk.eig[ei]-gs_en))>gamma_tol) continue; 
								intensity += fDvvDi[ei]/(omega_in-(exblk.eig[ei]-gs_en)+igamma);
							}
							rixs_peaks_kh[eloss_ind*nedos+n] += exp(-beta*gs_en)*pow(abs(intensity),2);
						}
					}
					// Exact solution
					if (pm.spec_solver == 1 || pm.spec_solver == 3) {
						int eloss_ind = floor((fs_en-gs_en)/(pm.em_energy)*nedos);
						rixs_em[eloss_ind] = fs_en-gs_en;
						vector<dcomp> csum(exen.size(),0);
						#pragma omp parallel for shared(fDvvDi,csum)
						for (size_t ei = 0; ei < exblk.nev; ++ei) {
							if (exblk.einrange[ei] == -1) continue;
							csum[exblk.einrange[ei]] += fDvvDi[ei];
						}
						#pragma omp parallel for shared(rixs_ab,rixs_peaks)
						for (size_t e = 0; e < exen.size(); e++) {
							if (abs(csum[e]) < TOL) continue;
							double ab_en = exen[e]-gs_en;
							int abind = floor((ab_en-ab_emin)/(ab_emax-ab_emin)*nedos);
							rixs_ab[abind] = ab_en;
							rixs_peaks.at(eloss_ind*nedos+abind) += exp(-beta*gs_en)*pow(abs(csum[e]),2);
						}
					}
				}
			}
		}}
	}

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";

	cout << "Writing results..." << endl;
	if (pm.spec_solver == 2 || pm.spec_solver == 3)
		// ed::write_vec(rixs_peaks_kh,nedos,nedos,"RIXS_"+pm.edge+"edge_"+pol_str(pm.pvin)+"_"+pol_str(pm.pvout)+"-full.csv");
		cout << "RIXS matrix wrote to RIXS_"+string(pm.edge)+"edge_"+pol_str(pm.pvin)+"_"+pol_str(pm.pvout)+"-kh.txt" << endl;
		for (auto & val : rixs_peaks_kh) val /= gsi.size();
		write_kh_RIXS(rixs_peaks_kh,rixs_em_kh,"RIXS_"+pm.edge+"edge_"+pol_str(pm.pvin)+"_"+pol_str(pm.pvout)+"-kh.txt");

	if (pm.spec_solver == 1 || pm.spec_solver == 3) {
		cout << endl;
		for (auto & val : rixs_peaks) val /= gsi.size();
		RIXS_peak_occupation(GS,EX,vecd({20}),rixs_ab,rixs_em,rixs_peaks,gsi,pm,gs_en,"top",true);
		write_RIXS(rixs_peaks,rixs_ab,rixs_em,pm.eloss,"RIXS_"+pm.edge+"edge_"+pol_str(pm.pvin)+"_"+pol_str(pm.pvout)+".txt");
	} 

	return;
}