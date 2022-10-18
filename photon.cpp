#include <complex>
#include <chrono>
#include <bitset>
#include "multiplet.hpp"
#include "photon.hpp"

using namespace std;
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
		std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
		initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

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
	if (hilbs.coord != "sqpl") {
		cout << "non square planar geometry occupation not coded yet" << endl;
		return;
	}

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
	int p = 5;
	for (auto &h : occ_totc) tot_h += h;
	cout << "Number of Holes: " << fixed << setprecision(p) << tot_h << endl;
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
		cout << fixed << setprecision(5);
		cout << "Ground State composition";
		for (size_t i = 0; i < ligNum; ++i) {
			cout << ", d" << 10-hilbs.num_vh+i << "L: " << dLweight[i];
		}
		cout << endl;
	}
	return dLweight;
}

double effective_delta(Hilbert& hilbs, int ligNum) {
	// Check effective delta, needs t = 0
	double d, dL;
	vector<bool> firstLh(ligNum,true);
	vector<double> distinct = ed::printDistinct(hilbs.get_all_eigval(),0.0,hilbs.hsize);
	for (size_t e = 0; e < distinct.size(); ++e) {
		vector<bindex> gei;
		for (size_t i = 0; i < hilbs.hblks.size(); ++i) {
		for (size_t j = 0; j < hilbs.hblks[i].size; ++j) {
			if (abs(hilbs.hblks[i].eig[j]-distinct[e]) < TOL) {
				gei.push_back({i,j});
			}
		}}
		vector<double> dLweight = wvfnc_weight(hilbs,gei,ligNum,false);
		for (size_t i = 0; i < ligNum; ++i) {
			if (abs(dLweight[i]-1) < TOL && firstLh[i]) {
				if (i == 0) d = distinct[e];
				if (i == 1) dL = distinct[e];
				firstLh[i] = false;
				cout << "Energy: " << distinct[e] << ", degeneracy: " << gei.size() << endl;
				wvfnc_weight(hilbs,gei,ligNum,true);
			}
		}
	}
	cout << "Effective delta: " << dL - d << endl;
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

void peak_occupation(Hilbert& hilbs, vecd const& peak_en, vecd const& energy, 
						vecd const& intensity, double ref_en, string mode) {
	// Calculate peak occupation using energy
	// If mode = top, it will use the peak_en as the top n number of peaks
	// TODO: Range of peak energy
	vecd peak_en_copy;
	if (mode == "list") peak_en_copy = peak_en;
	else if (mode == "top") {
		size_t top = (size_t)peak_en[0];
		vecd intensity_copy = intensity;
		vector<size_t> indices = ed::top_n(intensity_copy,top);
		peak_en_copy.resize(top);
		for (size_t i = 0; i < top; ++i) peak_en_copy[i] = energy[indices[i]];
	}
	for (size_t p = 0; p < peak_en_copy.size(); p++) {
		vector<bindex> peaks;
		for (auto &blk : hilbs.hblks) {	
			for (size_t i = 0; i < blk.size; ++i) {
				if (abs(blk.eig[i]-ref_en-peak_en_copy[p]) < 1e-5) 
					peaks.push_back(bindex(&blk-&hilbs.hblks[0],i));
			}
		}
		cout << "peak " << p << " energy: " << peak_en_copy[p];
		cout << ", degeneracy: " << peaks.size() << endl;
		occupation(hilbs,peaks);
	}
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
				* GS.Fsign(&vh,gs,1) * EX.Fsign(&ch,exs,1) * proj_pvec(vhqn.ml-chqn.ml,pvec);
		}
	}
	return;
}

void XAS(string input_dir, vector<double*>& SC, double* FG, double* CF, double SO, 
			bool HYB, int nedos, const PM& pm) {
	// Note: different diagonalize routine might yield different results, 
	// if the width of delta function is not small enough.
	// TODO: in small Hilbert space, parallelize won't yield correct result, rewrite it
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
	// Different SC for excited state
	double SC2[5]{6.054,0,10.7368,0,6.7152};
	SC[1] = SC2;
	calc_ham(EX,SC,FG,CF,SO);
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Diagonalizing Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	// #pragma omp parallel for schedule(static) num_threads(4)
	for (size_t g = 0; g < GS.hblks.size(); g++) {
		start = chrono::high_resolution_clock::now();
		cout << "Number of threads: " << mkl_get_max_threads() << endl;
		cout << "Diagonalizing block number: " << g << ", matrix size: " << GS.hblks[g].size << endl;
		GS.hblks[g].diagonalize();
		for (size_t i = 0; i < pow(GS.hblks[g].size,2); ++i) 
			if (abs(GS.hblks[g].eigvec[i]) < TOL) GS.hblks[g].eigvec[i] = 0;
		stop = chrono::high_resolution_clock::now();
		duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
		cout << "Run time = " << duration.count() << " ms\n" << endl;
	}

	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<bindex> gsi; // index for ground state
	// Calculate Partition function
	double Z = 0, emin = -25.0, emax = 25.0, SDegen = 0;
	for (size_t i = 0; i < GS.hblks.size(); ++i) {
	for (size_t j = 0; j < GS.hblks[i].size; ++j) {
		if (abs(GS.hblks[i].eig[j]-gs_en) < TOL) {
			Z += exp(-beta*GS.hblks[i].eig[j]);
			if (abs(GS.hblks[i].get_sz()) >= SDegen) SDegen = abs(GS.hblks[i].get_sz());
			gsi.push_back({i,j});
		}
	}}
	cout << "Grounds State Spin Quantum Number: " << SDegen << endl;
	cout << "Calculating Occupation (with degeneracy): " << gsi.size() << endl;
	occupation(GS,gsi);
	// effective_delta(GS,3);
	// return;

	// #pragma omp parallel for
	for (size_t e = 0; e < EX.hblks.size(); e++) {
		start = chrono::high_resolution_clock::now();
		cout << "Number of threads: " << mkl_get_max_threads() << endl;
		cout << "Diagonalizing block number: " << e << ", matrix size: " << EX.hblks[e].size << endl;
		EX.hblks[e].diagonalize();
		for (size_t i = 0; i < pow(EX.hblks[e].size,2); ++i) 
			if (abs(EX.hblks[e].eigvec[i]) < TOL) EX.hblks[e].eigvec[i] = 0;
		stop = chrono::high_resolution_clock::now();
		duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
		cout << "Run time = " << duration.count() << " ms\n" << endl;
	}
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Calculating cross section..." << endl;
	start = chrono::high_resolution_clock::now();
	vector<dcomp> blap;
	vecd xas_aben(nedos,0), xas_int(nedos,0);
	for (auto &exblk : EX.hblks) {	
		size_t last_gs_block = gsi[0].first;
		basis_overlap(GS,EX,bindex(last_gs_block,&exblk-&EX.hblks[0]),blap,pm);
		for (auto &g  : gsi) {
			Block& gsblk = GS.hblks[g.first];
			if (g.first != last_gs_block) {
				basis_overlap(GS,EX,bindex(g.first,&exblk-&EX.hblks[0]),blap,pm);
				last_gs_block = g.first;
			}
			cout << "GS blk number: " << &gsblk-&GS.hblks[0] << ", EX blk number: " << &exblk-&EX.hblks[0] << endl;
			#pragma omp parallel for reduction (vec_double_plus:xas_int)
			for (size_t ei = 0; ei < exblk.size; ++ei) {
				if (exblk.eig[ei]-gs_en < emin || exblk.eig[ei]-gs_en > emax) continue;
				dcomp cs = 0;
				for (size_t gj = 0; gj < gsblk.size; ++gj) {
					size_t gsind = g.second*gsblk.size+gj;
					if (gsblk.eigvec[gsind] == 0) continue;
					for (size_t ej = 0; ej < exblk.size; ++ej) {
						size_t exind = ei*exblk.size+ej, blapind = gj*exblk.size+ej;
						if (exblk.eigvec[exind] == 0 || blap[blapind] == 0.0) continue;
						cs +=  gsblk.eigvec[gsind] * exblk.eigvec[exind] * blap[blapind];
					}
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

	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";
	cout << "Writing results..." << endl << endl;
	write_XAS(xas_aben,xas_int,"xas_peaks.txt");

	// Calculate peak intensity
	peak_occupation(EX,vecd({2}),xas_aben,xas_int,gs_en,"top");

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
	for (size_t g = 0; g < GS.hblks.size(); g++) {
		GS.hblks[g].diagonalize();
		for (size_t i = 0; i < pow(GS.hblks[g].size,2); ++i) 
			if (abs(GS.hblks[g].eigvec[i]) < TOL) GS.hblks[g].eigvec[i] = 0;
	}
	#pragma omp parallel for
	for (size_t e = 0; e < EX.hblks.size(); e++) {
		EX.hblks[e].diagonalize();
		for (size_t i = 0; i < pow(EX.hblks[e].size,2); ++i) 
			if (abs(EX.hblks[e].eigvec[i]) < TOL) EX.hblks[e].eigvec[i] = 0;
	}
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
			exblk.einrange[ei] = ed::binary_search(exen,exblk.eig[ei],std::greater<int>(),
					[&](double a, double b){return abs(a-b) < TOL;});
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
					size_t gsind = g.second*gsblk.size+gj;
					if (gsblk.eigvec[gsind] == 0) continue;
					for (size_t ej = 0; ej < exblk.size; ++ej) {
						size_t exind = ei*exblk.size+ej, blapind = gj*exblk.size+ej;
						if (exblk.eigvec[exind] == 0 || blap[blapind] == 0.0) continue;
						csvi +=  gsblk.eigvec[gsind] * exblk.eigvec[exind] * blap[blapind];
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
		if (exblk.get_sz() != fsblk.get_sz()) {;} // Spin flip
		basis_overlap(GS,EX,bindex(&fsblk-&GS.hblks[0],&exblk-&EX.hblks[0]),blap,pm,true);
		cout << "FS blk: " << &fsblk-&GS.hblks[0] << ", EX blk: " << &exblk-&EX.hblks[0] << endl;
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
						size_t fsind = fi*fsblk.size+fj;
						if (fsblk.eigvec[fsind] == 0.0) continue;
						for (size_t ej = 0; ej < exblk.size; ++ej) {
							size_t exind = ei*exblk.size+ej, blapind = fj*exblk.size+ej;
							if (exblk.eigvec[exind] == 0.0 || blap[blapind] == 0.0) continue;
							csvf +=  fsblk.eigvec[fsind] * exblk.eigvec[exind] * blap[blapind];
						}
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
							if (!peak_occ_em.empty() && peak_occ_em.back() == fvind) continue;
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