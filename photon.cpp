#include <complex>
#include <chrono>
#include "multiplet.hpp"
#include "photon.hpp"

using namespace std;
typedef complex<double> dcomp;
#define LIM = 1E-7;

// Contains photon based spectroscopy

dcomp proj_pvec(int ml, vecd& pvec) {
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

void calc_ham(Hilbert& hilbs, double* SC, double* FG, double const& tenDQ, double const& lambda) {
	// Calculate Hamiltonian of the hilbert space
	int nd = 0;
	for (auto &at:hilbs.atlist) if(at.is_val && at.l == 2) nd = at.num_h - hilbs.is_ex;
	double del = 3+(nd-1)*(SC[0]-SC[2]*2.0/63-SC[4]*2.0/63); // ep = 3 eV
	// cout << "nd: " << nd << ", del: " << del << endl;
	calc_coulomb(hilbs,SC); 
	if (hilbs.SO) calc_SO(hilbs, lambda);
	if (hilbs.CF) calc_CF(hilbs, del, tenDQ, 6);
	if (hilbs.is_ex && hilbs.CV) calc_CV(hilbs, FG);
	return;
}

void write_XAS(vecd const& peaks, vecd const& y, string file_dir, bool print) {
	std::ofstream xasfile;
	xasfile.open(file_dir);
	xasfile << setw(15) << "peaks" << setw(15) << "intensity" << endl;
	if (print) cout << setw(15) << "peaks" << setw(15) << "intensity" << endl;
	for (int i = 0; i < y.size(); ++i) {
		if (y[i] != 0) {
			xasfile << setw(15) << peaks[i] << setw(15) << y[i] << endl;
			if (print) cout << setw(15) << peaks[i] << setw(15) << y[i] << endl;
 		}
	}
	xasfile.close();
}

void XAS(double* SC, double* FG, double tenDQ, double lambda, vecd& pvec, int nedos) {
	// Scatters from p6d9 to p5d10
	double beta = 0, racah_B = (SC[2]/49) - (5*SC[4]/441);

	cout << "polarization: ";
	for (auto &v:pvec) cout << v << " ";
	cout << endl << endl;

	cout << "Reading files..." << endl;
	auto start = chrono::high_resolution_clock::now();
	Hilbert GS("./INPUT",false);
	Hilbert EX("./INPUT",true);
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Assembling Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	calc_ham(GS,SC,FG,tenDQ,lambda);
	calc_ham(EX,SC,FG,tenDQ,lambda);
	ed::write_mat(EX.hblks[0].ham,EX.hblks[0].size,EX.hblks[0].size,"./exham.txt");
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Diagonalizing Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	for (auto & gsblk : GS.hblks) gsblk.diag_dsyev();
	// CAN PROBABLY DISCARD UNEEDED STATES HERE?
	for (auto & exblk : EX.hblks) exblk.diag_dsyev();

	ed::write_mat(GS.hblks[0].eigvec,GS.hblks[0].size,GS.hblks[0].size,"./gsev.txt");
	ed::write_mat(EX.hblks[0].eigvec,EX.hblks[0].size,EX.hblks[0].size,"./exev.txt");

	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;


	cout << "gs hsize: " << GS.hblks[0].size << ", ex hsize: " << EX.hblks[0].size << endl;

	// Calculate XAS matrix element, can be phased out
	dcomp* xas_mat = new dcomp[GS.hblks[0].size*EX.hblks[0].size]{0};
	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<pair<int,int>> gsi,exi; // index for ground state and excited states
	// Calculate Partition function
	double Z = 0, emin = -25.0, emax = 25.0;
	for (size_t i = 0; i < GS.hblks.size(); ++i) {
		for (size_t j = 0; j < GS.hblks[i].size; ++j) {
			Z += exp(-beta*GS.hblks[i].eig[j]);
			if (abs(GS.hblks[i].eig[j]-gs_en) < 1e-4) gsi.push_back({i,j});
		}
	} 

	cout << "GS energy: " << gs_en << ", degeneracy: " << gsi.size() << endl;
	cout << "partition function: " << Z << endl << endl;

	vecd xas_x,xas_y,xas_peaks;
	vector<dcomp> blap(GS.hsize*EX.hsize,0);
	for (double i = emin; i < emax; i += (emax-emin)/nedos) xas_x.push_back(i);
	xas_y.resize(xas_x.size(),0.0);
	xas_peaks.resize(xas_x.size(),0.0);

	int num_transitions = 0, half_orb = (EX.num_vorb+EX.num_corb)/2;

	cout << "Calculating cross section..." << endl;
	start = chrono::high_resolution_clock::now();

	// Calculating basis state overlap
	// Need to add in polarization
	for (size_t g = 0; g < GS.hsize; g++) {
		for (size_t e = 0; e < EX.hsize; e++) {
			ulli gs = GS.Hashback(g), exs = EX.Hashback(e);
			ulli ch = exs - (gs & exs), vh = gs - (gs & exs);
			int coi = EX.orbind(ch), voi = GS.atlist[coi].vind;
			if (!ed::is_pw2(vh) || !GS.atlist[voi].contains(vh)) continue;
			QN chqn = EX.atlist[coi].fast_qn(ch,half_orb,coi);
			QN vhqn = GS.atlist[voi].fast_qn(vh,half_orb,voi);
			if (chqn.spin != vhqn.spin || abs(vhqn.ml-chqn.ml) > 1) continue;
			blap[g*EX.hsize+e] = gaunt(1,chqn.ml,2,vhqn.ml)[1] * proj_pvec(vhqn.ml-chqn.ml,pvec)
								* GS.Fsign(&vhqn,gs,1) * EX.Fsign(&chqn,exs,1);
		}
	}

	for (auto &g : gsi) {
		Block gsblk = GS.hblks[g.first];
		for (auto &exblk : EX.hblks) {
		for (int ei = 0; ei < exblk.size; ++ei) {
			if (exblk.eig[ei]-gs_en < emin || exblk.eig[ei]-gs_en > emax) continue;
			dcomp cs = 0;
			for (int gj = 0; gj < gsblk.size; ++gj) {
				if (abs(gsblk.eigvec[g.second*gsblk.size+gj]) < 1e-7) continue;
				for (int ej = 0; ej < exblk.size; ++ej) {
					if (abs(exblk.eigvec[ei*exblk.size+ej]) < 1e-7) continue;
					cs +=  gsblk.eigvec[g.second*gsblk.size+gj] * exblk.eigvec[ei*exblk.size+ej]
						   * blap[gj*EX.hsize+ej];
				}
			}
			if (abs(cs) > 1e-7) {
				num_transitions += 1;
				xas_mat[g.second+GS.hblks[0].size*ei] += cs;
				xas_peaks[(int)floor((exblk.eig[ei]-gs_en-emin)/((emax-emin)/nedos))] = exblk.eig[ei]-gs_en;
				xas_y[(int)floor((exblk.eig[ei]-gs_en-emin)/((emax-emin)/nedos))] += exp(-beta*gs_en)*pow(abs(cs),2);
			}
		}}
	}

	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n";

	cout << "Number of Transitions: " << num_transitions << endl << endl;
	cout << "Writing results..." << endl << endl;
	// ed::write_mat(xas_mat,GS.hblks[0].size,EX.hblks[0].size,"./xas.txt");
	write_XAS(xas_peaks, xas_y, "xas_peaks.txt");

	return;
}
