#include <complex>
#include <chrono>
#include <omp.h>
#include "multiplet.hpp"
#include "photon.hpp"

using namespace std;
// Contains photon based spectroscopy

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

void calc_ham(Hilbert& hilbs, vector<double*>& SC, double* FG, double* CF, double const& SO) {
	// Calculate Hamiltonian of the hilbert space
	calc_coulomb(hilbs,SC); 
	if (hilbs.SO_on) calc_SO(hilbs,SO);
	if (hilbs.CF_on) calc_CF(hilbs,CF);
	if (hilbs.CV_on) calc_CV(hilbs,FG);
	if (hilbs.HYB_on) calc_HYB(hilbs,SC);
	return;
}

void occupation(Hilbert& hilbs, const vector<bindex>& si) {
	// Calculation occupation of orbitals, only valid when matrix is diagonalized
	for (auto& blk : hilbs.hblks) if (blk.eigvec == nullptr) throw runtime_error("matrix not diagonalized for occupation");
	if (hilbs.coord != "sqpl") throw runtime_error("non square planar geometry occupation not coded yet");

	int nvo = hilbs.num_vorb/2, nco = hilbs.num_corb/2;
	vecc U(nvo*nvo,0);
	U[0*nvo+0] = U[0*nvo+4] = {1/sqrt(2),0};
	U[1*nvo+2] = {1,0};
	U[2*nvo+0] = {0,1/sqrt(2)};
	U[2*nvo+4] = {0,-1/sqrt(2)};
	U[3*nvo+1] = {1/sqrt(2),0};
	U[3*nvo+3] = {-1/sqrt(2),0};
	U[4*nvo+1] = U[4*nvo+3] = {0,1/sqrt(2)};
	U[5*nvo+5] = U[5*nvo+7] = U[8*nvo+8] = U[8*nvo+10] = {1/sqrt(2),0};
	U[6*nvo+5] = U[9*nvo+8] = {0,1/sqrt(2)};
	U[6*nvo+7] = U[9*nvo+10] = {0,-1/sqrt(2)};
	U[7*nvo+6] = U[10*nvo+9] = {1,0};

	vector<double> occ_totc(nvo,0);
	for (auto &s  : si) {
		Block& blk = hilbs.hblks[s.first];
		// Pick an operator
		for (size_t c = 0; c < nvo; ++c) {
			vector<complex<double>> wvfnc(ed::choose(2*nvo,hilbs.num_vh-1),0);
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
						op_state -= op;
						int ind = 0, cnt = 0;
						// Figure out hash index for wavefunction
						for (size_t i = 0; i < nvo; ++i)
							if (1 << (i+nco) & op_state) ind += ed::choose(i,++cnt);
						for (size_t i = nvo; i < 2*nvo; ++i)
							if (1 << (i+2*nco) & op_state) ind += ed::choose(i,++cnt);
						// Figure out fermion sign
						op_state |= op;
						int fsgn = pow(-1,ed::count_bits(op_state/op));
						wvfnc[ind] += fsgn*sj*U[c*nvo+ci];
					}
				}
			}
			for (auto &w : wvfnc) occ_totc[c] += abs(pow(w,2))/si.size();
		}
	}

	double tot_h = 0;
	int p = 5;
	for (auto &h : occ_totc) tot_h += h;
	cout << "Number of Holes: " << tot_h << endl;
	cout << "dx2: " << fixed << setprecision(p) << occ_totc[0] << endl;
	cout << "dz2: " << fixed << setprecision(p) << occ_totc[1] << endl;
	cout << "dxy: " << fixed << setprecision(p) << occ_totc[2] << endl;
	cout << "dxz: " << fixed << setprecision(p) << occ_totc[3] << endl;
	cout << "dyz: " << fixed << setprecision(p) << occ_totc[4] << endl;
	cout << "pxx: " << fixed << setprecision(p) << occ_totc[5] << endl;
	cout << "pxy: " << fixed << setprecision(p) << occ_totc[6] << endl;
	cout << "pxz: " << fixed << setprecision(p) << occ_totc[7] << endl;
	cout << "pyx: " << fixed << setprecision(p) << occ_totc[8] << endl;
	cout << "pyy: " << fixed << setprecision(p) << occ_totc[9] << endl;
	cout << "pyz: " << fixed << setprecision(p) << occ_totc[10] << endl;
	return;
}

void write_XAS(vecd const& aben, vecd const& intensity, string file_dir, bool print) {
	std::ofstream xasfile;
	xasfile.open(file_dir);
	cout << fixed << setprecision(5);
	xasfile << setw(15) << "peaks (eV)" << setw(15) << "intensity" << endl;
	if (print) cout << setw(15) << "peaks" << setw(15) << "intensity" << endl;
	for (int i = 0; i < intensity.size(); ++i) {
		if (intensity[i] < TOL) continue;
		xasfile << setw(15) << aben[i] << setw(15) << intensity[i] << endl;
		if (print) cout << setw(15) << aben[i] << setw(15) << intensity[i] << endl;
	}
	xasfile.close();
}

void basis_overlap(Hilbert& GS, Hilbert& EX, bindex inds, vector<dcomp>& blap, 
					const PM& pm, bool pvout) {
	// Calculate basis state overlap base on the indices of the blocks
	blap.clear();
	size_t gbi = inds.first, exi = inds.second;
	blap = vector<dcomp>(GS.hblks[gbi].size*EX.hblks[exi].size,0);
	vecd pvec = pvout ? pm.pvout : pm.pvin;
	int cl = GS.atlist[0].l, vl = GS.atlist[GS.atlist[0].vind].l;
	int half_orb = (EX.num_vorb+EX.num_corb)/2;
	// #pragma omp parallel for
	for (size_t g = 0; g < GS.hblks[gbi].size; g++) {
		for (size_t e = 0; e < EX.hblks[exi].size; e++) {
			ulli gs = GS.Hashback(bindex(gbi,g)), exs = EX.Hashback(bindex(exi,e));
			ulli ch = exs - (gs & exs), vh = gs - (gs & exs);
			int coi = EX.orbind(ch), voi = GS.atlist[coi].vind;
			if (!ed::is_pw2(vh) || !GS.atlist[voi].contains(vh)) continue;
			QN chqn = EX.atlist[coi].fast_qn(ch,half_orb,coi);
			QN vhqn = GS.atlist[voi].fast_qn(vh,half_orb,voi);
			if (chqn.spin != vhqn.spin || abs(vhqn.ml-chqn.ml) > 1) continue;
			blap[g*EX.hblks[exi].size+e] = gaunt(cl,chqn.ml,vl,vhqn.ml)[1] 
				* GS.Fsign(&vhqn,gs,1) * EX.Fsign(&chqn,exs,1) * proj_pvec(vhqn.ml-chqn.ml,pvec);
		}
	}
	return;
}



void XAS(string input_dir, vector<double*>& SC, double* FG, double* CF, double SO, 
			bool HYB, int nedos, const PM& pm) {
	// Note: different diagonalize routine might yield different results, if the width of delta function is not small enough.
	double beta = 0;//, racah_B = (SC[2]/49) - (5*SC[4]/441);

	cout << "Reading files..." << endl;
	auto start = chrono::high_resolution_clock::now();
	Hilbert GS(input_dir,SC,FG,CF,SO,HYB,pm.edge,false);
	Hilbert EX(input_dir,SC,FG,CF,SO,HYB,pm.edge,true);
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Assembling Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	calc_ham(GS,SC,FG,CF,SO);
	calc_ham(EX,SC,FG,CF,SO);
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Diagonalizing Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	#pragma omp parallel for
	for (size_t g = 0; g < GS.hblks.size(); g++) GS.hblks[g].diag_dsyev();
	#pragma omp parallel for
	for (size_t e = 0; e < EX.hblks.size(); e++) EX.hblks[e].diag_dsyev();
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<bindex> gsi,exi; // index for ground state and excited states
	// Calculate Partition function
	double Z = 0, emin = -25.0, emax = 25.0;
	for (size_t i = 0; i < GS.hblks.size(); ++i) {
	for (size_t j = 0; j < GS.hblks[i].size; ++j) {
		if (abs(GS.hblks[i].eig[j]-gs_en) < TOL) {
			Z += exp(-beta*GS.hblks[i].eig[j]);
			gsi.push_back({i,j});
		}
	}}
	// cout << "Calculating Occupation (with degeneracy): " << gsi.size() << endl;
	// occupation(GS,gsi);

	cout << "Calculating cross section..." << endl;
	start = chrono::high_resolution_clock::now();

	vector<dcomp> blap;
	vecd xas_aben(nedos,0), xas_int(nedos,0);
	vector<pair<size_t,size_t>> peak_occ;
	for (auto &exblk : EX.hblks) {	
		size_t last_gs_block = gsi[0].first;
		basis_overlap(GS,EX,bindex(last_gs_block,&exblk-&EX.hblks[0]),blap,pm);
		for (auto &g  : gsi) {
			Block& gsblk = GS.hblks[g.first];
			if (g.first != last_gs_block) {
				basis_overlap(GS,EX,bindex(g.first,&exblk-&EX.hblks[0]),blap,pm);
				last_gs_block = g.first;
			}
			#pragma omp parallel for
			for (size_t ei = 0; ei < exblk.size; ++ei) {
				if (exblk.eig[ei]-gs_en < emin || exblk.eig[ei]-gs_en > emax) continue;
				dcomp cs = 0;
				for (size_t gj = 0; gj < gsblk.size; ++gj) {
					if (abs(gsblk.eigvec[g.second*gsblk.size+gj]) < TOL) continue;
					for (size_t ej = 0; ej < exblk.size; ++ej) {
						if (abs(exblk.eigvec[ei*exblk.size+ej]) < TOL && 
							abs(blap[gj*exblk.size+ej]) < TOL) continue;
						cs +=  gsblk.eigvec[g.second*gsblk.size+gj] * exblk.eigvec[ei*exblk.size+ej]
							   * blap[gj*exblk.size+ej];
					}
				}
				if (abs(cs) > TOL) { // Is this arbritrary?????	
					cout << "cross section: " << cs << endl;
					size_t peak_pos = round((exblk.eig[ei]-gs_en-emin)/((emax-emin)/nedos));
					xas_aben[peak_pos] = exblk.eig[ei]-gs_en;
					// Peaks position might be different if delta function is too wide, we can fix this by taking the smaller value
					xas_int[peak_pos] += exp(-beta*gs_en)*pow(abs(cs),2);
				}
				// Find occupation of peak here
				// if (abs(exblk.eig[ei]-gs_en+2.05315) < 0.01) {
				// 	// peak_occ.clear();
				// 	cout << "exblk pos: " << &exblk-&EX.hblks[0] << ", " << ei << ", pos: " << exblk.eig[ei]-gs_en << ", intensity: " << exp(-beta*gs_en)*pow(abs(cs),2) << endl;
				// 	peak_occ.push_back(bindex(&exblk-&EX.hblks[0],ei));
				// 	occupation(EX,vector<bindex>({bindex(&exblk-&EX.hblks[0],ei)}));
				// }
				// Find occupation of peak here
			}
		}
	}
	// cout << "number of peaks: " << peak_occ.size() << endl;
	// occupation(EX,peak_occ);

	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";
	cout << "Writing results..." << endl << endl;
	write_XAS(xas_aben,xas_int,"xas_peaks.txt");

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
				cout << setw(18) << ab[x];
				cout << setw(18) << em[y];
				cout << setw(18) << peaks[y*ab.size()+x] << endl;
			}
	 	}
	}
	rixsfile.close();
	return;
}

void RIXS(string input_dir, vector<double*>& SC, double* FG, double* CF, double SO, 
			bool HYB, int nedos, const PM& pm) {

	bool sf = true, nsf = true; // Spin Flip & No Spin Flip
	double beta = 0, hbar = 6.58e-16;
	// double Racah_B = (SC[2]/49) - (5*SC[4]/441);
	dcomp igamma(0,1*0.1);

	cout << "Reading files..." << endl;
	auto start = chrono::high_resolution_clock::now();
	Hilbert GS(input_dir,SC,FG,CF,SO,HYB,pm.edge,false);
	Hilbert EX(input_dir,SC,FG,CF,SO,HYB,pm.edge,true);
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Assembling Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	calc_ham(GS,SC,FG,CF,SO);
	calc_ham(EX,SC,FG,CF,SO);
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Diagonalizing Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	#pragma omp parallel for
	for (size_t g = 0; g < GS.hblks.size(); g++) GS.hblks[g].diag_dsyev();
	#pragma omp parallel for
	for (size_t e = 0; e < EX.hblks.size(); e++) EX.hblks[e].diag_dsyev();
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<bindex> gsi,exi; // index for ground state and excited states
	// Calculate Partition function
	double Z = 0, ab_emin = -25.0, ab_emax = 25.0, elm_min, elm_max;
	if (pm.eloss) { // Energy Loss
		elm_max = 15;
		elm_min = 0;
	} else { // Emission Energy
		elm_max = 15;
		elm_min = -15;
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

	cout << "ground state degeneracy: " << gsi.size() << ", excited state levels: " << exen.size() << endl;
	cout << "gsen: " << gs_en << ", exen: ";
	for (auto & e : exen) cout << e << ", ";
	cout << endl << "index: ";
	
	if (exen.size() > 1) { // Skip this if there are less than 1 excited state
		for (auto &exblk : EX.hblks) {
		for (size_t ei = 0; ei < exblk.size; ++ei) {
			if (exblk.einrange[ei] == -1) continue;
			exblk.einrange[ei] = ed::binary_search(exen,exblk.eig[ei]);
		}}
	}

	cout << "Calculating cross section..." << endl;
	start = chrono::high_resolution_clock::now();

	// Calculate and store <v|D|i>
	vecd rixs_em(nedos,0),rixs_ab(nedos,0),rixs_loss(nedos,0),rixs_peaks(nedos*nedos,0);
	for (int i = 0; i < nedos; ++i) rixs_ab[i] = ab_emin + i*(ab_emax-ab_emin)/nedos;
	vector<dcomp> blap, rixskern(gsi.size()*EX.hsize,0);
	vector<bindex> peak_occ_ab, peak_occ_em; // Calculate peak occupation
	for (auto &exblk : EX.hblks) {
		size_t last_gs_block = gsi[0].first;
		basis_overlap(GS,EX,bindex(last_gs_block,&exblk-&EX.hblks[0]),blap,pm);
		for (auto &g : gsi) {
			Block& gsblk = GS.hblks[g.first];
			if (g.first != last_gs_block) {
				basis_overlap(GS,EX,bindex(g.first,&exblk-&EX.hblks[0]),blap,pm);
				last_gs_block = g.first;
			}	
			#pragma omp parallel for
			for (size_t ei = 0; ei < exblk.size; ++ei) {
				if (exblk.eig[ei]-gs_en < ab_emin || exblk.eig[ei]-gs_en > ab_emax) continue;
				dcomp csvi = 0;
				for (size_t gj = 0; gj < gsblk.size; ++gj) {
					if (abs(gsblk.eigvec[g.second*gsblk.size+gj]) < TOL) continue;
					for (size_t ej = 0; ej < exblk.size; ++ej) {
						if (abs(exblk.eigvec[ei*exblk.size+ej]) < TOL 
							&& abs(blap[gj*exblk.size+ej]) < TOL) continue;
						csvi +=  gsblk.eigvec[g.second*gsblk.size+gj] * exblk.eigvec[ei*exblk.size+ej]
								 * blap[gj*exblk.size+ej];
					}
				}
				if (abs(csvi) > TOL) rixskern[(&g-&gsi[0])*EX.hsize+ei+exblk.f_ind] += csvi;
			}
		}
	}

	// Loop through final states, calculate <f|D|v><v|D|i>, add to spectra
	int gscnt = 0;
	bool is_gs = false;
	for (auto &fsblk : GS.hblks) {
	for (auto &exblk : EX.hblks) {
		basis_overlap(GS,EX,bindex(&fsblk-&GS.hblks[0],&exblk-&EX.hblks[0]),blap,pm,true);
		for (size_t fi = 0; fi < fsblk.size; ++fi) {
			if (pm.eloss && fsblk.eig[fi]-gs_en > elm_max) continue;
			vector<dcomp> csum(gsi.size()*exen.size(),0);
			if (abs(fsblk.eig[fi]-gs_en) < TOL) {
				is_gs = true;
				gscnt++;
			} else is_gs = false;
			#pragma omp parallel for
			for (size_t ei = 0; ei < exblk.size; ++ei) {
				if (pm.eloss) {if (exblk.eig[ei]-gs_en < ab_emin || exblk.eig[ei]-gs_en > ab_emax) continue;}
				else {if (exblk.eig[ei]-fsblk.eig[fi] > elm_max || exblk.eig[ei]-fsblk.eig[fi] < elm_min) continue;}
				dcomp csvf = 0;
				if (is_gs && pm.pvin == pm.pvout) csvf = rixskern[(gscnt-1)*EX.hsize+ei+exblk.f_ind];
				else {
					for (size_t fj = 0; fj < fsblk.size; ++fj) {
						if (abs(fsblk.eigvec[fi*fsblk.size+fj]) < TOL) continue;
						for (size_t ej = 0; ej < exblk.size; ++ej) {
							if (abs(exblk.eigvec[ei*exblk.size+ej]) < TOL || 
								abs(blap[fj*exblk.size+ej]) < TOL) continue;
							csvf +=  fsblk.eigvec[fi*fsblk.size+fj] * exblk.eigvec[ei*exblk.size+ej]
									 * blap[fj*exblk.size+ej];
						}
					}
				}
				if (abs(csvf) < TOL) continue;
				for (size_t gi = 0; gi < gsi.size(); ++gi) {
					if (abs(rixskern[gi*EX.hsize+ei+exblk.f_ind]) < TOL) continue;
					if (GS.hblks[gsi[gi].first].get_sz() != fsblk.get_sz()) {;} // Spin flip
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
					int elmind, abind = round((exen[e]-gs_en-ab_emin)/(ab_emax-ab_emin)*nedos);
					rixs_ab[abind] = exen[e]-gs_en;
					if (pm.eloss) {
						elmind = round((fsblk.eig[fi]-gs_en)/elm_max*nedos);
						rixs_em[elmind] = fsblk.eig[fi]-gs_en; //exen[e]-fsblk.eig[fi];
					} else {
						elmind = round((exen[e]-fsblk.eig[fi]-elm_min)/(elm_max-elm_min)*nedos);
						rixs_em[elmind] = exen[e]-fsblk.eig[fi];
					}
					rixs_peaks[elmind*nedos+abind] += exp(-beta*gs_en)*pow(abs(csum[g*exen.size()+e]),2);

					// Find occupation of f/v here
					if (exen[e]-gs_en > -0.5 && exen[e]-gs_en < -0.4) { // Absorptioin
						if (exen[e]-fsblk.eig[fi] > -2.4 && exen[e]-fsblk.eig[fi] < -2.3) { // Emission
							cout << "calculate peak" <<  endl;
							bindex fvind = bindex(&fsblk-&GS.hblks[0],fi);
							if (peak_occ_em.size() != 0 && peak_occ_em.back() == fvind) continue;
							cout << fvind.first << "," << fvind.second << endl;
							peak_occ_em.push_back(fvind);
						}
					} 
					// Find occupation of f/v here
				}
			}
		}
	}}

	cout << "ground state number:" << peak_occ_ab.size() << endl;
	occupation(GS,gsi);
	cout << "final state number:" << peak_occ_em.size() << endl;
	occupation(GS,peak_occ_em);


	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";
	cout << "Writing results..." << endl << endl;
	write_RIXS(rixs_peaks,rixs_ab,rixs_em,"rixs_peaks.txt",pm.eloss);
	return;
}