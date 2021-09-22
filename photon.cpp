#include <complex>
#include <ctime>
#include "multiplet.hpp"
#include "photon.hpp"

using namespace std;
typedef complex<double> dcomp;
#define LIM = 1E-7;

// Contains photon based spectroscopy

double* get_eigvec(int index, int size, double* eigvec_mat) {
	double* eigvec = new double[size];
	for (int i = 0; i < size; ++i) {
		eigvec[i] = eigvec_mat[index*size+i];
	}
	return eigvec;
}

dcomp proj_pvec(int ml, vector<double>& pvec) {
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
	cout << "nd: " << nd << ", del: " << del << endl;
	calc_coulomb(hilbs,SC);
	if (hilbs.SO) calc_SO(hilbs, lambda);
	if (hilbs.CF) calc_CF(hilbs, del, tenDQ, 6);
	if (hilbs.is_ex && hilbs.CV) calc_CV(hilbs, FG);
	return;
}

void XAS(double* SC, double* FG, double tenDQ, double lambda, vector<double>& pvec, int nedos) {
	// Scatters from p6d9 to p5d10
	double beta = 0, racah_B = (SC[2]/49) - (5*SC[4]/441);
	// ed::norm_vec(pvec);

	cout << "polarization: ";
	for (auto &v:pvec) cout << v << " ";
	cout << endl;


	Hilbert GS("./INPUT",false);
	Hilbert EX("./INPUT",true);
	calc_ham(GS,SC,FG,tenDQ,lambda);
	calc_ham(EX,SC,FG,tenDQ,lambda);
	for (auto & gsblk : GS.hblks) gsblk.diag_dsyev();
	// CAN PROBABLY DISCARD UNEEDED STATES HERE?
	for (auto & exblk : EX.hblks) exblk.diag_dsyev();

	ed::write_mat(GS.hblks[0].eigvec,GS.hblks[0].size,GS.hblks[0].size,"./gsev.txt");
	ed::write_mat(EX.hblks[0].eigvec,EX.hblks[0].size,EX.hblks[0].size,"./exev.txt");

	cout << "gs hsize: " << GS.hblks[0].size << ", ex hsize: " << EX.hblks[0].size << endl;
	cout << "ground state eigenstates: ";
	ed::printDistinct(GS.hblks[0].eig,GS.hblks[0].size);
	cout << endl << "excited state eigenstates: ";
	ed::printDistinct(EX.hblks[0].eig,EX.hblks[0].size);
	cout << endl;

	// Calculate XAS matrix element, can be phased out
	dcomp* xas_mat = new dcomp[GS.hblks[0].size*EX.hblks[0].size]{0};
	double gs_en = GS.hblks[0].eig[0];
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	vector<pair<int,int>> gsi,exi; // index for ground state and excited states
	// Calculate Partition function
	double Z = 0, emin = -15.0, emax = 15.0;
	for (size_t i = 0; i < GS.hblks.size(); ++i) {
		for (size_t j = 0; j < GS.hblks[i].size; ++j) {
			Z += exp(-beta*GS.hblks[i].eig[j]);
			if (abs(GS.hblks[i].eig[j]-gs_en) < 1e-7) gsi.push_back({i,j});
		}
	} 
	cout << "GS energy: " << gs_en << ", degeneracy: " << gsi.size() << endl;

	cout << "partition function: " << Z << endl;
	vector<double> xas_x,xas_y, xas_temp;
	for (double i = emin; i < emax; i += (emax-emin)/nedos) xas_x.push_back(i);
	xas_y.resize(xas_x.size(),0.0);
	xas_temp.resize(xas_x.size(),0.0);

	int num_transitions = 0, half_orb = (EX.num_vorb+EX.num_corb)/2;

	// for (auto &g : gsi) {
	// 	cout << "Eigen vectors" << endl;
	// 	cout << "g first: " << g.first << ", g second: " << g.second << endl;
	// 	Block gsblk = GS.hblks[g.first];
	// 	for (int gj = 0; gj < gsblk.size; ++gj) {
	// 		if (abs(gsblk.eigvec[g.second*gsblk.size+gj]) < 1e-7) continue;
	// 		ulli gs = GS.Hashback(gj);
	// 		cout << "gj: " << gj << ", basis state: " << gs << ", " << bitset<16>(gs) << ", value: " << gsblk.eigvec[g.second*gsblk.size+gj] << endl;
	// 	}
	// }

	// for (int gj = 0; gj < GS.hblks[0].size; ++gj) {
	// 	for (int ej = 0; ej < EX.hblks[0].size; ++ej) {
	// 		ulli gs = GS.Hashback(gj), exs = EX.Hashback(ej), ch = exs - (gs & exs), vh = gs - (gs & exs);
	// 		int coi = EX.orbind(ch), voi = GS.atlist[coi].vind;
	// 		if (!ed::is_pw2(vh) || !GS.atlist[voi].contains(vh)) continue;
	// 		QN chqn = EX.atlist[coi].get_qn(ch,half_orb,coi);
	// 		QN vhqn = GS.atlist[voi].get_qn(vh,half_orb,voi);
	// 		if (chqn.spin != vhqn.spin || abs(vhqn.ml-chqn.ml) > 1) continue;
	// 		if (gj == 0 && ej == 10) {
	// 			cout << "ex: " << bitset<16>(exs) << ", gs: " << bitset<16>(gs) << endl;
	// 			cout << "ch: " << bitset<16>(ch) << ", vh: " << bitset<16>(vh) << endl;
	// 			cout << "Gaunt coeff: " << gaunt(1,chqn.ml,2,vhqn.ml)[1] << endl;
	// 			cout << "Projection Operator: " << proj_pvec(vhqn.ml-chqn.ml,pvec) << endl;
	// 			cout << "Fermion Sign: " << GS.Fsign(&vhqn,gs,1) * EX.Fsign(&chqn,exs,1) << endl;
	// 			cout << "result: " << gaunt(1,chqn.ml,2,vhqn.ml)[1] * proj_pvec(vhqn.ml-chqn.ml,pvec) * GS.Fsign(&vhqn,gs,1) * EX.Fsign(&chqn,exs,1);
	// 		}
	// 		xas_mat[gj+GS.hblks[0].size*ej] = gaunt(1,chqn.ml,2,vhqn.ml)[1] * proj_pvec(vhqn.ml-chqn.ml,pvec) * GS.Fsign(&vhqn,gs,1) * EX.Fsign(&chqn,exs,1);
	// 	}
	// }
	// ed::write_mat(xas_mat,GS.hblks[0].size,EX.hblks[0].size,"./xas.txt");
	// gsi = vector<pair<int,int>>{{0,23}};
	// cout << "bitstate: " << bitset<16>(GS.Hashback(23)) << endl;
	for (auto &g : gsi) {
		Block gsblk = GS.hblks[g.first];
		for (auto &exblk : EX.hblks) {
		for (int ei = 0; ei < EX.hblks[0].size; ++ei) {
			// cout << endl << "energy diff: " << exblk.eig[ei]-gs_en;
			if (exblk.eig[ei]-gs_en < emin || exblk.eig[ei]-gs_en > emax) continue;
			dcomp cs = 0;
			for (int gj = 0; gj < gsblk.size; ++gj) {
				if (abs(gsblk.eigvec[g.second*gsblk.size+gj]) < 1e-7) continue;
				ulli gs = GS.Hashback(gj);
				for (int ej = 0; ej < exblk.size; ++ej) {
					if (abs(exblk.eigvec[ei*exblk.size+ej]) < 1e-7) continue;
					ulli exs = EX.Hashback(ej), ch = exs - (gs & exs), vh = gs - (gs & exs);
					int coi = EX.orbind(ch), voi = GS.atlist[coi].vind;
					if (!ed::is_pw2(vh) || !GS.atlist[voi].contains(vh)) continue;
					QN chqn = EX.atlist[coi].get_qn(ch,half_orb,coi);
					QN vhqn = GS.atlist[voi].get_qn(vh,half_orb,voi);
					if (chqn.spin != vhqn.spin || abs(vhqn.ml-chqn.ml) > 1) continue;  	// Selection rule

					
					// Debug Block
					// cout << endl << "one of the cross section" << endl;
					// // // cout << endl << "gs: " << gj << ", " << bitset<16>(gs) << ", ex: " << ej << ", " << bitset<16>(exs) << endl;
					// // // cout << "GS Energy: " << gsblk.eig[g.second] << ", EX energy: " << exblk.eig[ei] << endl;
					// cout << "gs Eigenvector: " << gsblk.eigvec[g.second*gsblk.size+gj] << ", state: " << bitset<16>(gs) << endl;
					// cout << "ex Eigenvector: " << exblk.eigvec[ei*exblk.size+ej] << ", state: " << bitset<16>(exs) << endl;
					// cout << "ex: " << bitset<16>(exs) << ", gs: " << bitset<16>(gs) << endl;
					// cout << "ch: " << bitset<16>(ch) << ", vh: " << bitset<16>(vh) << endl;
					// // cout << "Fermion Sign: " << GS.Fsign(&vhqn,gs,1) << ", gaunt coeff: " << gaunt(1,chqn.ml,2,vhqn.ml)[1] << ", proj ml: " << vhqn.ml-chqn.ml << endl;
					// // cout << "ch ml: " << chqn.ml << ", ch spin: " << chqn.spin << ", vh ml: " << vhqn.ml << ", vh spin: " << vhqn.spin << ", voi: "<< voi << endl;
					// cout << "dcs: " << gsblk.eigvec[g.second*gsblk.size+gj] * exblk.eigvec[ei*exblk.size+ej]
					// 		* gaunt(1,chqn.ml,2,vhqn.ml)[1] * proj_pvec(vhqn.ml-chqn.ml,pvec)
					// 		* GS.Fsign(&vhqn,gs,1) * EX.Fsign(&chqn,exs,1) << endl;
					// Debug Block

					cs += gsblk.eigvec[g.second*gsblk.size+gj] * exblk.eigvec[ei*exblk.size+ej]
							* gaunt(1,chqn.ml,2,vhqn.ml)[1] * proj_pvec(vhqn.ml-chqn.ml,pvec)
							* GS.Fsign(&vhqn,gs,1) * EX.Fsign(&chqn,exs,1);// Something funny with the fermion sign
				}
			}
			if (abs(cs) > 1e-7) {
				num_transitions += 1;
				xas_mat[g.second+GS.hblks[0].size*ei] += cs;
				xas_temp[(int)floor((exblk.eig[ei]-gs_en-emin)/((emax-emin)/nedos))] = exblk.eig[ei]-gs_en;
				xas_y[(int)floor((exblk.eig[ei]-gs_en-emin)/((emax-emin)/nedos))] += exp(-beta*gs_en)*pow(abs(cs),2);///Z;
			}
		}}
	}

	cout << "Number of Transitions: " << num_transitions << endl;

	// Normalize XAS intensity
	ed::write_mat(xas_mat,GS.hblks[0].size,EX.hblks[0].size,"./xas.txt");
	// for (auto &x : xas_x) cout << setw(5) << x << " ";
	// cout << endl;
	// cout.precision(3);
	for (int i = 0; i < xas_y.size(); ++i) {
		if (xas_y[i] != 0) {
			cout << "x: " << setw(5) << xas_temp[i] << ", e: " << setw(10) << xas_y[i] << endl;
 		}
	}
	for (auto &x : xas_y) cout << x << " ";
	cout << endl;	

	return;
}
