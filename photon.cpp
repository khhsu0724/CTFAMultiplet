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

void occupation(Hilbert& hilbs) {
	// Calculation occupation of orbitals, only valid when matrix is diagonalized
	for (auto& blk : hilbs.hblks) if (blk.eigvec == nullptr) throw runtime_error("matrix not diagonalized for occupation");
	if (hilbs.coord != "sqpl") throw runtime_error("non square planar geometry occupation not coded yet");
	double gs_en = hilbs.hblks[0].eig[0];
	for (auto &gsb : hilbs.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<pair<size_t,size_t>> gsi;
	for (size_t i = 0; i < hilbs.hblks.size(); ++i) {
	for (size_t j = 0; j < hilbs.hblks[i].size; ++j) {
		if (abs(hilbs.hblks[i].eig[j]-gs_en) < TOL) gsi.push_back({i,j});
	}}

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

	cout << "Ground state: " << gs_en << ", degeneracy: " << gsi.size() << endl;
	vector<double> occ_totc(nvo,0);
	for (auto &g  : gsi) {
		Block& gsblk = hilbs.hblks[g.first];
		// Pick an operator
		for (size_t c = 0; c < nvo; ++c) {
			vector<complex<double>> wvfnc(ed::choose(2*nvo,hilbs.num_vh-1),0);
			for (size_t ci = 0; ci < nvo; ++ci) {
				if (abs(U[c*nvo+ci]) < TOL) continue;
				// #pragma omp parallel for
				for (size_t gj = 0; gj < gsblk.size; ++gj) {
					double gsgj = gsblk.eigvec[g.second*gsblk.size+gj];
					if (abs(gsgj) < TOL) continue;
					// Operate on spin up and spin down
					for (int sud = nco; sud <= 2*nco+nvo; sud += nco+nvo) {
						ulli op_state = hilbs.Hashback(gj), op;
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
						wvfnc[ind] += fsgn*gsgj*U[c*nvo+ci];
					}
				}
			}
			for (auto &w : wvfnc) occ_totc[c] += abs(pow(w,2))/gsi.size();
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
	xasfile << setw(15) << "peaks (eV)" << setw(15) << "intensity" << endl;
	if (print) cout << setw(15) << "peaks" << setw(15) << "intensity" << endl;
	for (int i = 0; i < intensity.size(); ++i) {
		if (intensity[i] < TOL) continue;
		xasfile << setw(15) << aben[i] << setw(15) << intensity[i] << endl;
		if (print) cout << setw(15) << aben[i] << setw(15) << intensity[i] << endl;
	}
	xasfile.close();
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
	for (auto & gsblk : GS.hblks) gsblk.diag_dsyev();
	// cout << "Calculating Occupation: " << endl;
	// occupation(GS);
	// return;


	for (auto & exblk : EX.hblks) exblk.diag_dsyev();
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<pair<size_t,size_t>> gsi,exi; // index for ground state and excited states
	// Calculate Partition function
	double Z = 0, emin = -25.0, emax = 25.0;
	for (size_t i = 0; i < GS.hblks.size(); ++i) {
	for (size_t j = 0; j < GS.hblks[i].size; ++j) {
		if (abs(GS.hblks[i].eig[j]-gs_en) < TOL) {
			Z += exp(-beta*GS.hblks[i].eig[j]);
			gsi.push_back({i,j});
		}
	}}
	
	vecd xas_aben(nedos,0), xas_int(nedos,0);
	vector<dcomp> blap(GS.hsize*EX.hsize,0);
	int half_orb = (EX.num_vorb+EX.num_corb)/2;

	cout << "Calculating basis overlap..." << endl;
	start = chrono::high_resolution_clock::now();
	// Calculating basis state overlap
	int cl = GS.atlist[0].l, vl =  GS.atlist[GS.atlist[0].vind].l;
	#pragma omp parallel for
	for (size_t g = 0; g < GS.hsize; g++) {
		for (size_t e = 0; e < EX.hsize; e++) {
			ulli gs = GS.Hashback(g), exs = EX.Hashback(e);
			ulli ch = exs - (gs & exs), vh = gs - (gs & exs);
			int coi = EX.orbind(ch), voi = GS.atlist[coi].vind;
			if (!ed::is_pw2(vh) || !GS.atlist[voi].contains(vh)) continue;
			QN chqn = EX.atlist[coi].fast_qn(ch,half_orb,coi);
			QN vhqn = GS.atlist[voi].fast_qn(vh,half_orb,voi);
			if (chqn.spin != vhqn.spin || abs(vhqn.ml-chqn.ml) > 1) continue;
			blap[g*EX.hsize+e] = gaunt(cl,chqn.ml,vl,vhqn.ml)[1] * GS.Fsign(&vhqn,gs,1)
				* EX.Fsign(&chqn,exs,1) * proj_pvec(vhqn.ml-chqn.ml,pm.pvin);
		}
	}

	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";

	cout << "Calculating cross section..." << endl;
	start = chrono::high_resolution_clock::now();

	for (auto &g  : gsi) {
		Block& gsblk = GS.hblks[g.first];
		for (auto &exblk : EX.hblks) {	
		// #pragma omp parallel for shared(xas_aben,xas_int)
		for (size_t ei = 0; ei < exblk.size; ++ei) {
			if (exblk.eig[ei]-gs_en < emin || exblk.eig[ei]-gs_en > emax) continue;
			dcomp cs = 0;
			for (size_t gj = 0; gj < gsblk.size; ++gj) {
				if (abs(gsblk.eigvec[g.second*gsblk.size+gj]) < TOL) continue;
				for (size_t ej = 0; ej < exblk.size; ++ej) {
					if (abs(exblk.eigvec[ei*exblk.size+ej]) < TOL && abs(blap[gj*EX.hsize+ej]) < TOL) continue;
					cs +=  gsblk.eigvec[g.second*gsblk.size+gj] * exblk.eigvec[ei*exblk.size+ej]
						   * blap[gj*EX.hsize+ej];
				}
			}
			if (abs(cs) > TOL) {
				size_t peak_pos = round((exblk.eig[ei]-gs_en-emin)/((emax-emin)/nedos));
				xas_aben[peak_pos] = exblk.eig[ei]-gs_en;
				// Peaks position might be different if delta function is too wide, we can fix this by taking the smaller value
				xas_int[peak_pos] += exp(-beta*gs_en)*pow(abs(cs),2);
			}
		}}
	}

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
	
	// DEBUG: Spin Flip
	int gscnt = 0, excnt = 0;
	vecd gbss(GS.hsize,0), gsts(GS.hsize,0);
	cout << "Number of core orbitals: " << GS.num_corb << endl;
	for (int i = 0; i < GS.hsize; ++i) {
		// sum up all atoms qn.spin in gbss
		ulli s = GS.Hashback(i);
		for (int j = GS.num_corb/2; j < (GS.num_vorb+GS.num_corb)/2; ++j) if (s & 1<<j) gbss[i] -= 0.5;
		gbss[i] = 2 * gbss[i] + 0.5 * GS.num_vh;
		// cout << i << ", state: " << bitset<28>(s) << ", spin: " << gbss[i] << endl;
	}
	// cout << "non zero of ham" << endl;
	// for (int i = 0; i < GS.hsize; ++i) {
	// 	for (int j = 0; j < GS.hsize; ++j) {
	// 		// cout << gbss[i] << "," << gbss[j] << "," << GS.hblks[0].ham[i+GS.hsize*j] << endl;
	// 		if (gbss[i] != gbss[j] && GS.hblks[0].ham[i+GS.hsize*j] != 0) {
	// 			cout << "i: " << bitset<16>(GS.Hashback(i)) << ", j: " << bitset<16>(GS.Hashback(j)) << ", entry: " << GS.hblks[0].ham[i+GS.hsize*j] << endl; 
	// 		}
	// 	}
	// }
	// // DEBUG

	cout << "Diagonalizing Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	for (auto & gsblk : GS.hblks) gsblk.diag_dsyev();
	for (auto & exblk : EX.hblks) exblk.diag_dsyev();
	stop = chrono::high_resolution_clock::now();
	// ed::write_mat(GS.hblks[0].eigvec,GS.hsize,GS.hsize,"./GS_eigvec.txt");
	// ed::write_mat(EX.hblks[0].eigvec,EX.hsize,EX.hsize,"./EX_eigvec.txt");
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;


	// DEBUG
	// cout << "non zero of eigvec" << endl;
	// for (int i = 0; i < GS.hsize; ++i) {
	// 	for (int j = 0; j < GS.hsize; ++j) {
	// 		// cout << gbss[i] << "," << gbss[j] << "," << GS.hblks[0].ham[i+GS.hsize*j] << endl;
	// 		if (gbss[i] != gbss[j] && GS.hblks[0].eigvec[i+GS.hsize*j] != 0) {
	// 			cout << "i: " << bitset<16>(GS.Hashback(i)) << ", j: " << bitset<16>(GS.Hashback(j)) << ", entry: " << GS.hblks[0].eigvec[i+GS.hsize*j] << endl; 
	// 		}
	// 	}
	// }
	// DEBUG

	cout << "Ground state eigenvalues" << endl;
	ed::printDistinct(GS.hblks[0].eig,GS.hblks[0].size);
	cout << "Core Hole state eigenvalues" << endl;
	ed::printDistinct(EX.hblks[0].eig,EX.hblks[0].size);


	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<pair<size_t,size_t>> gsi,exi; // index for ground state and excited states
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
	// double emin = ab_emin > em_emin ? ab_emin + gs_en : em_emin + gs_en;
	// double emax = ab_emax > em_emax ? ab_emax + gs_en : em_emax + gs_en;
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

	for (size_t i = 0; i < EX.hblks[0].size; ++i) cout << EX.hblks[0].einrange[i] << ", ";
	cout << endl << endl;

	// Figure out total spin of ground and final state

	for (auto &gsblk : GS.hblks) {
	for (size_t gi = 0; gi < gsblk.size; ++gi) {
		for (size_t gj = 0; gj < gsblk.size; ++gj) {
			// Loop through eigenvectors to find out total spin of each state Sz
			if (abs(gsblk.eigvec[gi*gsblk.size+gj]) > TOL)  {
				gsts[gi] += pow(gsblk.eigvec[gi*gsblk.size+gj],2) * gbss[gj];
				// cout << "index: " << gj << ", eigval: " << gsblk.eigvec[gi*gsblk.size+gj];
				// cout << ", state: " << bitset<16>(GS.Hashback(gj));
				// cout << ", spin: " << gbss[gj] << endl;
			}
		}
		// cout << "Total spin of this state: " << gsts[gi] << ", eigenerg: " << gsblk.eig[gi] << endl << endl;
		gscnt++;
	}}
	gscnt = 0;
	// return;

	vecd rixs_em(nedos,0),rixs_ab(nedos,0),rixs_loss(nedos,0),rixs_peaks(nedos*nedos,0);
	for (int i = 0; i < nedos; ++i) rixs_ab[i] = ab_emin + i*(ab_emax-ab_emin)/nedos;
	vector<dcomp> ivlap(GS.hsize*EX.hsize,0), fvlap(GS.hsize*EX.hsize,0);
	int half_orb = (EX.num_vorb+EX.num_corb)/2;

	cout << "Calculating cross section..." << endl;
	start = chrono::high_resolution_clock::now();
	// Calculating basis state overlap
	int cl = GS.atlist[0].l, vl =  GS.atlist[GS.atlist[0].vind].l;
	// #pragma omp parallel for
	for (size_t g = 0; g < GS.hsize; g++) {
		for (size_t e = 0; e < EX.hsize; e++) {
			ulli gs = GS.Hashback(g), exs = EX.Hashback(e);
			ulli ch = exs - (gs & exs), vh = gs - (gs & exs);
			int coi = EX.orbind(ch), voi = GS.atlist[coi].vind;
			if (!ed::is_pw2(vh) || !GS.atlist[voi].contains(vh)) continue;
			QN chqn = EX.atlist[coi].fast_qn(ch,half_orb,coi);
			QN vhqn = GS.atlist[voi].fast_qn(vh,half_orb,voi);
			if (chqn.spin != vhqn.spin || abs(vhqn.ml-chqn.ml) > 1) continue;
			double overlap = gaunt(cl,chqn.ml,vl,vhqn.ml)[1] * GS.Fsign(&vhqn,gs,1) 
								 * EX.Fsign(&chqn,exs,1);
			ivlap[g*EX.hsize+e] = overlap * proj_pvec(vhqn.ml-chqn.ml,pm.pvin); // Phase out
			fvlap[g*EX.hsize+e] = overlap * proj_pvec(vhqn.ml-chqn.ml,pm.pvout);
		}
	}
	// Brute force test
	// for (auto &g : gsi) {
	// 	// cout << "Brute force" << endl;
	// 	Block& gsblk = GS.hblks[g.first];
	// 	for (auto &fsblk : GS.hblks) {
	// 	for (int fi = 0; fi < fsblk.size; ++fi) {
	// 		vector<dcomp> csum(exen.size(),0);
	// 		for (auto &exblk : EX.hblks) {
	// 		for (int ei = 0; ei < exblk.size; ++ei) {
	// 			if (exblk.eig[ei]-gs_en < ab_emin || exblk.eig[ei]-gs_en > ab_emax) continue;
	// 			if (exblk.eig[ei]-fsblk.eig[fi] < em_emin || exblk.eig[ei]-fsblk.eig[fi] > em_emax) continue;
	// 			dcomp csvi(0,0), csvf(0,0);
	// 			for (int gj = 0; gj < gsblk.size; ++gj) {
	// 				if (abs(gsblk.eigvec[g.second*gsblk.size+gj]) < TOL) continue;
	// 				for (int ej = 0; ej < exblk.size; ++ej) {
	// 					if (abs(exblk.eigvec[ei*exblk.size+ej]) < TOL/* && abs(ivlap[gj*EX.hsize+ej]) < TOL*/) continue;
	// 					csvi +=  gsblk.eigvec[g.second*gsblk.size+gj] * exblk.eigvec[ei*exblk.size+ej]
	// 							 * ivlap[gj*EX.hsize+ej];
	// 				}
	// 			}
	// 			for (int fj = 0; fj < fsblk.size; ++fj) {
	// 				if (abs(fsblk.eigvec[fi*fsblk.size+fj]) < TOL) continue;
	// 				for (int ej = 0; ej < exblk.size; ++ej) {
	// 					if (abs(exblk.eigvec[ei*exblk.size+ej]) < TOL && abs(fvlap[fj*EX.hsize+ej]) < TOL) continue;
	// 					csvf +=  fsblk.eigvec[fi*fsblk.size+fj] * exblk.eigvec[ei*exblk.size+ej]
	// 							 * fvlap[fj*EX.hsize+ej];
	// 				}
	// 			}
	// 			// cout << "index: " << exblk.eind[ei] << ", energy: " << exblk.eig[ei] << endl;
	// 			csum[exblk.einrange[ei]] += conj(csvf) * csvi;
	// 			// if (abs(conj(csvf)*csvi) > TOL) cout << "i: " << g.second << ", f: " << fi << ", v:"  << ei << ", csvi: " << csvi << ", csvf: " << csvf <<  ", cross section: " << conj(csvf) * csvi << ", energy: " << exen[exblk.eind[ei]] << endl;
	// 		}}
	// 		// Add to spectra
	// 		for (int e = 0; e < exen.size(); e++) {
	// 			if (abs(csum[e]) < TOL) continue;
	// 			int abind = round((exen[e]-gs_en-ab_emin)/(ab_emax-ab_emin)*nedos);
	// 			int emind = round((exen[e]-fsblk.eig[fi]-em_emin)/(em_emax-em_emin)*nedos);
	// 			rixs_ab[abind] = exen[e]-gs_en;
	// 			rixs_em[emind] = exen[e]-fsblk.eig[fi];
	// 			// cout << "gi: " << g.second << ", fi: " << fi << ", gse: " << gs_en << ", ve: " << exen[e] << ", fe: " << fsblk.eig[fi] << ", cs: " << csum[e] << endl;
	// 			rixs_peaks[emind*nedos+abind] += exp(-beta*gs_en)*pow(abs(csum[e]),2);
	// 		}
	// 	}}
	// }

	// stop = chrono::high_resolution_clock::now();
	// duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	// cout << "Run time = " << duration.count() << " ms\n";
	// cout << "Writing results..." << endl << endl;
	// // ed::write_vec(rixs_peaks,nedos,nedos,"rsixs_mat.txt");
	// write_RIXS(rixs_peaks,rixs_ab,rixs_em,"rixs_peaks.txt");
	// return;

	// Calculate and store <v|D|i> and elastic scattering line
	vector<dcomp> rixskern(gsi.size()*EX.hsize,0);
	gscnt = 0, excnt = 0;
	for (auto &g : gsi) {
		Block& gsblk = GS.hblks[g.first];	
		// vector<dcomp> csum(exen.size(),0);
		for (auto &exblk : EX.hblks) {
		// #pragma omp parallel for
		for (size_t ei = 0; ei < exblk.size; ++ei) {
			if (exblk.eig[ei]-gs_en < ab_emin || exblk.eig[ei]-gs_en > ab_emax) continue;
			dcomp csvi = 0;
			for (size_t gj = 0; gj < gsblk.size; ++gj) {
				if (abs(gsblk.eigvec[g.second*gsblk.size+gj]) < TOL) continue;
				for (size_t ej = 0; ej < exblk.size; ++ej) {
					if (abs(exblk.eigvec[ei*exblk.size+ej]) < TOL && abs(ivlap[gj*EX.hsize+ej]) < TOL) continue;
					csvi +=  gsblk.eigvec[g.second*gsblk.size+gj] * exblk.eigvec[ei*exblk.size+ej]
							 * ivlap[gj*EX.hsize+ej];
				}
			}
			if (abs(csvi) > TOL) rixskern[gscnt*EX.hsize+ei] += csvi;
		}}
		gscnt++;
	}

	// Loop through final states, calculate <f|D|v><v|D|i>, add to spectra
	gscnt = 0;
	bool is_gs = false;
	for (auto &fsblk : GS.hblks) {
	for (size_t fi = 0; fi < fsblk.size; ++fi) {
		// Check if spin flips, this can be done by checking blocks instead of calculating
		excnt = 0;
		vector<dcomp> csum(gsi.size()*exen.size(),0);
		if (pm.eloss && fsblk.eig[fi]-gs_en > elm_max) continue;
		if (abs(fsblk.eig[fi]-gs_en) < TOL) {
			is_gs = true;
			gscnt++;
		} else is_gs = false;
		for (auto &exblk : EX.hblks) {
		// #pragma omp parallel for
		for (size_t ei = 0; ei < exblk.size; ++ei) {
			if (pm.eloss) {if (exblk.eig[ei]-gs_en < ab_emin || exblk.eig[ei]-gs_en > ab_emax) continue;}
			else {if (exblk.eig[ei]-fsblk.eig[fi] > elm_max || exblk.eig[ei]-fsblk.eig[fi] < elm_min) continue;}
			dcomp csvf = 0;
			if (is_gs && pm.pvin == pm.pvout) csvf = rixskern[(gscnt-1)*EX.hsize+ei];// Wrong!!!!!!!
			else {
				for (size_t fj = 0; fj < fsblk.size; ++fj) {
					if (abs(fsblk.eigvec[fi*fsblk.size+fj]) < TOL) continue;
					for (size_t ej = 0; ej < exblk.size; ++ej) {
						if (abs(exblk.eigvec[ei*exblk.size+ej]) < TOL || abs(fvlap[fj*EX.hsize+ej]) < TOL) continue;
						csvf +=  fsblk.eigvec[fi*fsblk.size+fj] * exblk.eigvec[ei*exblk.size+ej]
								 * fvlap[fj*EX.hsize+ej];
					}
				}
			}
			if (abs(csvf) < TOL) continue;
			for (size_t gi = 0; gi < gsi.size(); ++gi) {
				// Probably have to calculate spin flip here
				if (abs(rixskern[gi*EX.hsize+ei]) < TOL) continue;
				/*if (abs(gsts[gi]-gsts[fi]) < TOL)*/ // Spin Flip
				csum[gi*exen.size()+exblk.einrange[ei]] += conj(csvf) * rixskern[gi*EX.hsize+ei];
				// cout << "gi: " << gi << ", gi spin: " << gsts[gi] << ", fi: " << fi << ", fi spin: " << gsts[fi] << ", spin diff: " << gsts[gi]-gsts[fi] << endl;
				// cout << "eind: " << exblk.eind[ei] << ", eigen: " << exblk.eig[ei] << ", en: " << conj(csvf) * rixskern[gi*EX.hsize+ei] << endl;
			}
		}}

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
		// #pragma omp parallel for
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
				// cout << "abe: " << exen[e]-gs_en-ab_emin << ", eme: " << exen[e]-fsblk.eig[fi]-em_emin << endl;
				// cout << "abind: " << abind << ", emind: " << emind << endl;
				// cout << "absorption: " << exen[e]-gs_en << ", loss: " << fsblk.eig[fi]-gs_en << ", cs: " << exp(-beta*gs_en)*pow(abs(csum[g*exen.size()+e]),2) << endl;
				// for (int e2 = 0; e2 < exen.size(); e2++) {
				// 	if (abs(csum[g*exen.size()+e2]) < TOL) continue;
				// 	cs += (csum[g*exen.size()+e2])/(exen.at(e)-exen.at(e2)-igamma);
				// }
				// rixs_peaks[emind*nedos+abind] += exp(-beta*gs_en)*pow(abs(cs),2);
			}
		}
	}}
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";
	cout << "Writing results..." << endl << endl;
	// ed::write_vec(rixs_peaks,nedos,nedos,"rixs_mat.txt");
	write_RIXS(rixs_peaks,rixs_ab,rixs_em,"rixs_peaks.txt",pm.eloss);
	return;
}