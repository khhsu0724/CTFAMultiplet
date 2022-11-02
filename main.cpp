#include <cstdlib> 
#include <ctime> 
#include <stdlib.h>
#include <bitset>
#include <regex>
#include <fstream>
#include "multiplet.hpp"
#include "photon.hpp"

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif

using namespace std;

void tb() {
	// This is pretty out of date
	ofstream myfile;
    myfile.open ("./tb.txt");
    HParam hparam;
	double SC2[5] = {0,0,1*49,0,0.078*441};
	double SC1[3] = {0,0,0};
	vector<double*> SC;
	SC.emplace_back(SC1);
	SC.emplace_back(SC2);
	hparam.SC = SC;
	double CF[5] = {1,1,-2.0/3,-2.0/3,-2.0/3};
	copy(CF,CF+5,hparam.CF);
	for (double m = 0.5; m <= 50; m += 0.5) {
		Hilbert input("./INPUT",hparam,"L",false);
		double CF_TB[5]{0};
		for (int c = 0; c < 5; c++) CF_TB[c] = CF[c] * m;
		calc_coulomb(input,SC);
		calc_CF(input,CF_TB);
		input.hblks[0].diagonalize();
		vector<double> unique_eig = ed::printDistinct(input.hblks[0].eig,
							input.hblks[0].eig[0],input.hblks[0].size,false);
		auto min_eig = *min_element(unique_eig.begin(), unique_eig.end());
		sort(unique_eig.begin(), unique_eig.end());
		int i = 1;
		for (auto& e : unique_eig) {
			myfile << e - min_eig;
			if (++i <= unique_eig.size()) myfile << " ";
		}
		myfile << endl;
	}
}

bool read_num(string line, double* tgt, int pnum) {
	try {
		string p;
		bool skip = false;
		int ncnt = 0;
		for (int s = 0; s < line.size() && !skip; ++s) {
			if (line[s] == '#') skip = true;
			else if (line[s] == '=') p = "";
			else if (line[s] != ',' && line[s] != ' ' && line[s] != '	') p.push_back(line[s]);
			else if (regex_match(p,regex(R"(^([+-]?(?:[[:d:]]+\.?|[[:d:]]*\.[[:d:]]+))(?:[Ee][+-]?[[:d:]]+)?$)"))) {
				if (ncnt < pnum) tgt[ncnt] = stod(p);
				else throw invalid_argument("Too many input parameter");
				ncnt++;
				p = "";
			} else if (p.size() != 0) throw invalid_argument("Invalid number");
			else p = "";
		}
		if (p.size() != 0 && regex_match(p,regex(R"(^([+-]?(?:[[:d:]]+\.?|[[:d:]]*\.[[:d:]]+))(?:[Ee][+-]?[[:d:]]+)?$)"))) {
			if (ncnt < pnum) tgt[ncnt] = stod(p);
			else throw invalid_argument("Too many input parameter");
			ncnt++;
			p = "";
		}
		if (ncnt != pnum) throw invalid_argument("Wrong number of input parameter");
	} catch (const exception &ex) {
		cerr << ex.what() << "\n";
		exit(0);
	}
	return true;
}

bool read_bool(string line, bool& tgt) {
	try {
		string p;
		bool skip = false;
		for (int s = 0; s < line.size() && !skip; ++s) {
			if (line[s] == '#') skip = true;
			else if (line[s] == '=') p = "";
			else if (line[s] != ',' && line[s] != ' ' && line[s] != '	') p.push_back(line[s]);
			else if (regex_match(p,regex(".*[a-zA-Z]+.*"))) {
				transform(p.begin(),p.end(),p.begin(),::toupper);
				if (p == "TRUE") tgt = true;
				else if (p == "FALSE") tgt = false;
				else throw invalid_argument("Invalid input (True/False)");
			} else if (p.size() != 0) throw invalid_argument("Invalid input (True/False)");
			else p = "";
		}
		if (p.size() != 0 && regex_match(p,regex(".*[a-zA-Z]+.*"))) {
			transform(p.begin(),p.end(),p.begin(),::toupper);
			if (p == "TRUE") tgt = true;
			else if (p == "FALSE") tgt = false;
			else throw invalid_argument("Invalid input (True/False)");
		}
	} catch (const exception &ex) {
		cerr << ex.what() << "\n";
		exit(0);
	}
	return true;
}

void read_input(string IDIR, PM& pm, HParam& hparam, bool& overwrite) {
	string line;
	ifstream input(IDIR);
	try {
		if (input.is_open()) {
			bool read_control = false, read_photon = false;
			int atind = 0;
			while (getline(input,line)) {
				if (line[0] == '#') continue;
				if (line[0] == '/') {
					read_control = false;
					read_photon = false;
					continue;
				}
				if (read_control) {
					string p;
					bool skip = false;
					for (int s = 0; s < line.size() && !skip; ++s) {
						if (line[s] != ' ' && line[s] != '	' && line[s] != '=') p.push_back(line[s]);
						else {
							transform(p.begin(),p.end(),p.begin(),::toupper);
							if (p == "SO") skip = read_num(line.substr(s+1,line.size()-1),&hparam.SO,1);
							else if (p == "SC1") skip = read_num(line.substr(s+1,line.size()-1),hparam.SC[0],3);
							else if (p == "SC2") skip = read_num(line.substr(s+1,line.size()-1),hparam.SC[1],5);
							else if (p == "FG") skip = read_num(line.substr(s+1,line.size()-1),hparam.FG,4);
							else if (p == "CF") skip = read_num(line.substr(s+1,line.size()-1),hparam.CF,5);
							else if (p == "HYB") skip = read_num(line.substr(s+1,line.size()-1),&hparam.HYB,1);
							else if (p == "MLCT") skip = read_num(line.substr(s+1,line.size()-1),&hparam.MLdelta,1);
							else if (p == "OVERWRITE") skip = read_bool(line.substr(s+1,line.size()-1),overwrite);
							p = "";
						}
						if (line[s] == '#') skip = true;
					}
				}
				if (read_photon) {
					string p;
					bool skip = false;
					double pvtmp[3]{0};
					for (int s = 0; s < line.size() && !skip; ++s) {
						if (line[s] != ' ' && line[s] != '	' && line[s] != '=') p.push_back(line[s]);
						else {
							transform(p.begin(),p.end(),p.begin(),::toupper);
							if (p == "XAS") skip = read_bool(line.substr(s+1,line.size()-1),pm.XAS);
							else if (p == "RIXS") skip = read_bool(line.substr(s+1,line.size()-1),pm.RIXS);
							else if (p == "PVIN") {
								skip = read_num(line.substr(s+1,line.size()-1),pvtmp,3);
								for (int i = 0; i < 3; ++i) pm.pvin[i] = pvtmp[i];
							}
							else if (p == "PVOUT") {
								skip = read_num(line.substr(s+1,line.size()-1),pvtmp,3);
								for (int i = 0; i < 3; ++i) pm.pvout[i] = pvtmp[i];
							}
							else if (p == "EDGE") {
								size_t eqpos = line.find('=');
								string input = line.substr(eqpos+1,line.find('#')-eqpos-1);
								eqpos = input.find_first_not_of(' ');
								size_t string_size = input.find_last_not_of(' ')-eqpos+1;
								pm.edge = input.substr(eqpos,string_size);
								skip = true;
							}
							else if (p == "ELOSS") {
								bool el;
								skip = read_bool(line.substr(s+1,line.size()-1),el);
								pm.set_eloss(el);
							}
							else if (p == "NEDOS") {
								double nedos;
								skip = read_num(line.substr(s+1,line.size()-1),&nedos,1);
								pm.nedos = (int)nedos;
							}
							else if (p == "SPINFLIP") skip = read_bool(line.substr(s+1,line.size()-1),pm.spin_flip);
							p = "";
						}
						if (line[s] == '#') skip = true;
					}
				}
				if (line == "&CONTROL") {
					read_control = true;
					read_photon = false;
				}
				else if (line == "&PHOTON") {
					read_photon = true;
					read_control = false;
				}
			}
		} else throw invalid_argument("Cannot open " + IDIR + " file");
	} catch (const exception &ex) {
		cerr << ex.what() << "\n";
		exit(0);
	}
	return;
}

bool if_photon_file_exist(const string edge, const string work_dir, bool is_XAS, 
							const vecd& pvin, const vecd& pvout = vecd(3,0)) {
	string pfname = work_dir;
	if (is_XAS) pfname += "/XAS_"+edge+"edge_"+pol_str(pvin)+".txt";
	else pfname += "/RIX_"+edge+"edge_"+pol_str(pvin)+"_"+pol_str(pvout)+".txt";
	ifstream pfile(pfname.c_str());
	if (pfile.good()) {
		cout << "file: " << pfname << " exists, skipping..." << endl;
	}
	return pfile.good();
}

void process_hilbert_space(string input_dir, Hilbert& GS, Hilbert& EX, const PM& pm, const HParam& hparam) {
	// Assembling Hamiltonian and Diagonalizing, returns partition function
	double beta = 0;
	cout << "Assembling Hamiltonian..." << endl;
	auto start = chrono::high_resolution_clock::now();
	calc_ham(GS,hparam);
	// Different SC for excited state
	double SC2[5]{6.054,0,10.7368,0,6.7152};
	std::copy(hparam.SC[1],hparam.SC[1]+5,SC2);
	calc_ham(EX,hparam);
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Diagonalizing Hamiltonian..." << endl;
	auto diag_start = chrono::high_resolution_clock::now();
	// #pragma omp parallel for schedule(static) num_threads(4)
	cout << "Diagonalizing Ground State" << endl;
	for (size_t g = 0; g < GS.hblks.size(); g++) {
		cout << "Number of threads: " << mkl_get_max_threads() << endl;
		start = chrono::high_resolution_clock::now();
		cout << "Diagonalizing block number: " << g << ", matrix size: " << GS.hblks[g].size << endl;
		GS.hblks[g].diagonalize();
		for (size_t i = 0; i < pow(GS.hblks[g].size,2); ++i) 
			if (abs(GS.hblks[g].eigvec[i]) < TOL) GS.hblks[g].eigvec[i] = 0;
		stop = chrono::high_resolution_clock::now();
		duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
		cout << "Run time = " << duration.count() << " ms\n" << endl;
	}
	double gs_en = GS.hblks[0].eig[0];
	vector<bindex> gsi;
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.size; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
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
	cout << "Grounds State Spin Quantum Number (S): " << SDegen << endl;
	cout << "Calculating Occupation (with degeneracy): " << gsi.size() << endl;
	occupation(GS,gsi);
	cout << "Ground State composition: ";
	wvfnc_weight(GS,gsi,3,true);
	cout << endl;
	// effective_delta(GS,3);
	// return;

	// #pragma omp parallel for
	cout << "Diagonalizing Core-Hole State" << endl;
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
	auto diag_stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(diag_stop - diag_start);
	cout << "Total diagonalize run time = " << duration.count() << " ms\n" << endl;
	return;
}

int main(int argc, char** argv){
	char CurrentPath[FILENAME_MAX];
	if (!GetCurrentDir(CurrentPath, sizeof(CurrentPath))) {return errno;}
	string work_dir(CurrentPath);
	cout << "The working directory is: " << work_dir << endl;
	// Reading and Checking input parameters
	auto run_start = chrono::high_resolution_clock::now();
	cout << "Reading files..." << endl;
	auto start = chrono::high_resolution_clock::now();
	PM pm;
	HParam hparam; 
	bool overwrite = false;
	string IDIR = "./INPUT";
	if (argc == 2) IDIR = string(argv[1]);
	read_input(IDIR,pm,hparam,overwrite);
	if (pm.edge == "K") {
		if (hparam.FG[2] != 0 || hparam.FG[3] != 0) 
			throw invalid_argument("invalid FG input");
	} else if (pm.edge != "L") throw invalid_argument("invalid edge input");
	// U = F^0 + 4*F^2 + 36*F^4

	cout << "Input parameters" << endl;
	cout << "SO: " << hparam.SO << "eV" << endl;
	cout << "SC1: ";
	for (int i = 0; i < 3; ++i) cout << hparam.SC[0][i] << ", ";
	cout << endl << "SC2: ";
	for (int i = 0; i < 5; ++i) cout << hparam.SC[1][i] << ", ";
	cout << endl << "FG: ";
	for (int i = 0; i < 4; ++i) cout << hparam.FG[i] << ", ";
	cout << endl << "CF: ";
	for (int i = 0; i < 5; ++i) cout << hparam.CF[i] << ", ";
	cout << endl;
	cout << "HYB: " << hparam.HYB << endl; 
	cout << "delta: " << hparam.MLdelta << endl; 
	if (pm.RIXS) {
		if (pm.eloss) cout << "Using energy loss" << endl;
		else cout << "using omega/omega" << endl;
	}
	cout << "photon edge: " << pm.edge << endl;
	cout << "pvin: ";
	for (int i = 0; i < 3; ++i) cout << pm.pvin[i] << ", ";
	cout << endl << "pvout: ";
	for (int i = 0; i < 3; ++i) cout << pm.pvout[i] << ", ";
	cout << endl;

	
	// Adjust Slater Condor Parameter
	hparam.SC[0][2] *= 25;
	hparam.SC[1][2] *= 49;
	hparam.SC[1][4] *= 441;

	if (!pm.XAS && !pm.RIXS) return 0; // Return if no photon methods
	// Reading and Checking input parameters
	Hilbert GS(IDIR,hparam,pm.edge,false);
	Hilbert EX(IDIR,hparam,pm.edge,true);
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;
	process_hilbert_space(IDIR,GS,EX,pm,hparam);

	// Calculate Cross Sections
	vector<double> pvin = pm.pvin, pvout = pm.pvout;
	if (pm.XAS) {
		for (size_t in = 0; in < 3; ++in) {
			fill(pm.pvin.begin(), pm.pvin.end(), 0);
			if (pvin[in]) pm.pvin[in] = 1;
			else continue;
			if (if_photon_file_exist(pm.edge,work_dir,true,pvin) && !overwrite) continue;
			cout << "Doing photon polarization: ";
			for (auto e : pm.pvin) cout << (int)e << " ";
			cout << endl;
			XAS(GS,EX,pm);
		}
	}
	if (pm.RIXS) {
		for (size_t in = 0; in < 3; ++in) {
			fill(pm.pvin.begin(), pm.pvin.end(), 0);
			if (pvin[in]) pm.pvin[in] = 1;
			else continue;
			for (size_t out = 0; out < 3; ++out) {
				fill(pm.pvout.begin(), pm.pvout.end(), 0);
				if (pvout[out]) pm.pvout[out] = 1;
				else continue;
				if (if_photon_file_exist(pm.edge,work_dir,false,pvin,pvout) && !overwrite) continue;
				cout << "Doing photon polarization (in): ";
				for (auto e : pm.pvin) cout << (int)e << " ";
				cout << "(out): ";
				for (auto e : pm.pvout) cout << (int)e << " ";
				cout << endl;
				RIXS(GS,EX,pm);
			}
		}
	}

	auto run_stop = chrono::high_resolution_clock::now();
	cout << endl << "Total elapsed time: " << endl;
	cout << ed::format_duration(chrono::duration_cast<chrono::milliseconds>(run_stop-run_start)) << endl;
	return 0;
}