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

// STACK TRACE DEBUG
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}
// STACK TRACE DEBUG

using namespace std;

void tb() {
	// Tanabe-Sugano This is pretty out of date
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
		input.hblks[0].diagonalize(2);
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

void stot(string p, double& tgt) {tgt = stod(p);};
void stot(string p, float& tgt) {tgt = stof(p);};
void stot(string p, int& tgt) {tgt = stoi(p);};

template<typename T> 
bool read_num(string line, T* tgt, int pnum) {
	try {
		string p;
		bool skip = false;
		int ncnt = 0;
		for (int s = 0; s < line.size() && !skip; ++s) {
			if (line[s] == '#') skip = true;
			else if (line[s] == '=') p = "";
			else if (line[s] != ',' && line[s] != ' ' && line[s] != '	') p.push_back(line[s]);
			else if (regex_match(p,regex(R"(^([+-]?(?:[[:d:]]+\.?|[[:d:]]*\.[[:d:]]+))(?:[Ee][+-]?[[:d:]]+)?$)"))) {
				if (ncnt < pnum) stot(p,tgt[ncnt]);
				else throw invalid_argument("Too many input parameter");
				ncnt++;
				p = "";
			} else if (p.size() != 0) throw invalid_argument("Invalid number");
			else p = "";
		}
		if (p.size() != 0 && regex_match(p,regex(R"(^([+-]?(?:[[:d:]]+\.?|[[:d:]]*\.[[:d:]]+))(?:[Ee][+-]?[[:d:]]+)?$)"))) {
			if (ncnt < pnum) stot(p,tgt[ncnt]);
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
					bool skip = false, SC2EX_read = false;
					for (int s = 0; s < line.size() && !skip; ++s) {
						if (line[s] != ' ' && line[s] != '	' && line[s] != '=') p.push_back(line[s]);
						else {
							transform(p.begin(),p.end(),p.begin(),::toupper);
							if (p == "SO") skip = read_num(line.substr(s+1,line.size()-1),&hparam.SO,1);
							else if (p == "SC1") skip = read_num(line.substr(s+1,line.size()-1),hparam.SC[0],3);
							else if (p == "SC2") {
								skip = read_num(line.substr(s+1,line.size()-1),hparam.SC[1],5);
								if (!SC2EX_read) read_num(line.substr(s+1,line.size()-1),hparam.SC2EX,5);
							}
							else if (p == "SC2EX") {
								SC2EX_read = true;
								skip = read_num(line.substr(s+1,line.size()-1),hparam.SC2EX,5);
							}
							else if (p == "FG") skip = read_num(line.substr(s+1,line.size()-1),hparam.FG,4);
							else if (p == "CF") skip = read_num(line.substr(s+1,line.size()-1),hparam.CF,5);
							else if (p == "HYB") skip = read_bool(line.substr(s+1,line.size()-1),hparam.HYB);
							else if (p == "TPD") skip = read_num(line.substr(s+1,line.size()-1),&hparam.tpd,1);
							else if (p == "TPP") skip = read_num(line.substr(s+1,line.size()-1),&hparam.tpp,1);
							else if (p == "MLCT") skip = read_num(line.substr(s+1,line.size()-1),&hparam.MLdelta,1);
							else if (p == "EFFDEL") skip = read_bool(line.substr(s+1,line.size()-1),hparam.effective_delta);
							else if (p == "BLOCK") skip = read_bool(line.substr(s+1,line.size()-1),hparam.block_diag);
							else if (p == "DIAG") skip = read_num(line.substr(s+1,line.size()-1),&hparam.diag_option,1);
							else if (p == "SITEOCC") skip = read_bool(line.substr(s+1,line.size()-1),hparam.print_site_occ);
							else if (p == "OVERWRITE") skip = read_bool(line.substr(s+1,line.size()-1),overwrite);
							else if (p == "GSNEV") skip = read_num(line.substr(s+1,line.size()-1),&hparam.gs_nev,1);
							else if (p == "EXNEV") skip = read_num(line.substr(s+1,line.size()-1),&hparam.ex_nev,1);
							p = "";
						}
						if (line[s] == '#') skip = true;
					}
				}
				if (read_photon) {
					string p;
					bool skip = false;
					double pvtmp[3]{0};
					double abrange_temp[2]{0};
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
							else if (p == "NEDOS") skip = read_num(line.substr(s+1,line.size()-1),&pm.nedos,1);
							else if (p == "SPINFLIP") skip = read_bool(line.substr(s+1,line.size()-1),pm.spin_flip);
							else if (p == "AB") {
								skip = read_num(line.substr(s+1,line.size()-1),abrange_temp,2);
								for (int i = 0; i < 2; ++i) pm.ab_range[i] = abrange_temp[i];
							}
							else if (p == "EM") skip = read_num(line.substr(s+1,line.size()-1),&pm.em_energy,1);
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

bool if_photon_file_exist(const string edge, const string work_dir, bool is_XAS, bool overwrite,
							const vecd& pvin, const vecd& pvout = vecd(3,0)) {
	string pfname = work_dir;
	if (is_XAS) pfname += "/XAS_"+edge+"edge_"+pol_str(pvin)+".txt";
	else pfname += "/RIXS_"+edge+"edge_"+pol_str(pvin)+"_"+pol_str(pvout)+".txt";
	ifstream pfile(pfname.c_str());
	if (pfile.good()) {
		cout << "file: " << pfname << " exists";
		if (overwrite) cout << ", overwritting..." << endl;
		else cout << ", skipping..." << endl;
	}
	return pfile.good();
}

double calculate_effective_delta(const string& input_dir, HParam& hparam, const PM& pm) {
	double inp_tpd = hparam.tpd, inp_tpp = hparam.tpp;
	hparam.tpd = 0;
	hparam.tpp = 0;
	Hilbert GS(input_dir,hparam,pm.edge,false);
	if (!GS.cluster->lig_per_site) {
		cout << "No ligands, skipping calculate delta" << endl;
		return 0;
	}
	calc_ham(GS,hparam);
	for (size_t g = 0; g < (GS.hblks.size()+1)/2; g++) {
		GS.hblks[g].diagonalize(abs(hparam.MLdelta)*30);
	}
	hparam.tpd = inp_tpd;
	hparam.tpp = inp_tpp;
	return effective_delta(GS,2);
}

void process_hilbert_space(Hilbert& GS, Hilbert& EX, HParam& hparam, PM& pm) {
	// Assembling Hamiltonian and Diagonalizing, returns partition function
	double beta = 0;
	// Assemble GS Hamiltonian
	cout << "Assembling GS Hamiltonian..." << endl;
	auto start = chrono::high_resolution_clock::now();
	calc_ham(GS,hparam);
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Diagonalizing Hamiltonian..." << endl;
	auto diag_start = chrono::high_resolution_clock::now();
	cout << "Diagonalizing Ground State" << endl;
	if (!hparam.gs_nev) {
		if (pm.RIXS) hparam.gs_nev = 400;
		else hparam.gs_nev = 10 * GS.tot_site_num();
	}
	for (size_t g = 0; g < GS.hblks.size(); g++) {
		start = chrono::high_resolution_clock::now();
		cout << "Diagonalizing block number: " << g << ", matrix size: " << GS.hblks[g].size;
		cout << ", gs nev: " << hparam.gs_nev << endl;
		GS.hblks[g].diagonalize(hparam.gs_nev);
		stop = chrono::high_resolution_clock::now();
		duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
		cout << "Run time = " << duration.count() << " ms\n" << endl;
	}
	auto diag_stop = chrono::high_resolution_clock::now();
	auto diag_duration = chrono::duration_cast<chrono::milliseconds>(diag_stop - diag_start);

	double gs_en = GS.hblks[0].eig[0];
	vector<bindex> gsi;
	for (auto &gsb : GS.hblks) for (size_t i = 0; i < gsb.nev; ++i) if (gsb.eig[i] <= gs_en) gs_en = gsb.eig[i];
	// Calculate Partition function
	double Z = 0, emin = -25.0, emax = 25.0, SDegen = 0;
	for (size_t i = 0; i < GS.hblks.size(); ++i) {
	for (size_t j = 0; j < GS.hblks[i].nev; ++j) {
		if (abs(GS.hblks[i].eig[j]-gs_en) < TOL) {
			Z += exp(-beta*GS.hblks[i].eig[j]);
			if (abs(GS.hblks[i].get_sz()) >= SDegen) SDegen = abs(GS.hblks[i].get_sz());
			gsi.push_back({i,j});
		}
	}}
	if (hparam.block_diag) cout << "Grounds State Spin Quantum Number (S): " << SDegen << endl;
	cout << "Calculating Occupation (with degeneracy): " << gsi.size() << endl;
	occupation(GS,gsi,true,"./gs_occ.txt");
	wvfnc_weight(GS,gsi,2,true);
	cout << endl;

	// Assemble Core Hole Hamiltonian
	cout << "Assembling Core Hole Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	hparam.SC.pop_back();
	hparam.SC.push_back(hparam.SC2EX);
	calc_ham(EX,hparam);
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "size: " << EX.hblks[0].size << endl;
	cout << "Run time = " << duration.count() << " ms\n" << endl;

	cout << "Diagonalizing Core-Hole Hamiltonian" << endl;
	diag_start = chrono::high_resolution_clock::now();
	if (!hparam.ex_nev) {
		if (pm.edge == "K") hparam.ex_nev = 200;
		else if (pm.edge == "L3") hparam.ex_nev = 500;
		else if (pm.edge == "L") hparam.ex_nev = 1500;
	}
	for (size_t e = 0; e < EX.hblks.size(); e++) {
		start = chrono::high_resolution_clock::now();
		cout << "Diagonalizing block number: " << e << ", matrix size: " << EX.hblks[e].size;
		cout << ", ex nev: " << hparam.ex_nev << endl;
		EX.hblks[e].diagonalize(hparam.ex_nev);
		stop = chrono::high_resolution_clock::now();
		duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
		cout << "Run time = " << duration.count() << " ms\n" << endl;
	}
	diag_stop = chrono::high_resolution_clock::now();
	diag_duration += chrono::duration_cast<chrono::milliseconds>(diag_stop - diag_start);
	
	double ex_min_en = EX.hblks[0].eig[0];
	vector<bindex> exmin_ind;
	for (auto &exb : EX.hblks) for (size_t i = 0; i < exb.nev; ++i) if (exb.eig[i] <= ex_min_en) ex_min_en = exb.eig[i];
	SDegen = 0;
	for (size_t i = 0; i < EX.hblks.size(); ++i) {
	for (size_t j = 0; j < EX.hblks[i].nev; ++j) {
		if (abs(EX.hblks[i].eig[j]-ex_min_en) < TOL) {
			if (hparam.block_diag && abs(EX.hblks[i].get_sz()) >= SDegen) 
				SDegen = abs(EX.hblks[i].get_sz());
			exmin_ind.push_back({i,j});
		}
	}}
	// cout << "Band Gap: " << ex_min_en - gs_en << "eV" << endl;
	// auto all_eig = EX.get_all_eigval(true);
	// ed::printDistinct(all_eig,0.0,all_eig.size(),true);

	if (hparam.block_diag && EX.edge == "K") cout << "Grounds State Spin Quantum Number (S): " << SDegen << endl;
	cout << "Calculating Occupation (with degeneracy): " << exmin_ind.size() << endl;
	occupation(EX,exmin_ind);
	wvfnc_weight(EX,exmin_ind,2);
	cout << endl << "Total diagonalize run time: " << ed::format_duration(diag_duration) << endl << endl;
	return;
}

int main(int argc, char** argv){
	signal(SIGSEGV, handler);
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
	} else if (pm.edge != "L" && pm.edge != "L3") throw invalid_argument("invalid edge input");
	// U = F^0 + 4*F^2 + 36*F^4

	cout << "Input parameters" << endl;
	cout << "SO: " << hparam.SO << "eV" << endl;
	cout << "SC1: ";
	for (int i = 0; i < 3; ++i) cout << hparam.SC[0][i] << ", ";
	cout << endl << "SC2 for GS: ";
	for (int i = 0; i < 5; ++i) cout << hparam.SC[1][i] << ", ";
	cout << endl << "SC2 for EX: ";
	for (int i = 0; i < 5; ++i) cout << hparam.SC2EX[i] << ", ";
	cout << endl << "FG: ";
	for (int i = 0; i < 4; ++i) cout << hparam.FG[i] << ", ";
	cout << endl << "CF: ";
	for (int i = 0; i < 5; ++i) cout << hparam.CF[i] << ", ";
	cout << endl;
	if (hparam.HYB) {
		cout << "Hybridization turned on: "; 
		cout << "tpd = " << hparam.tpd; 
		cout << ", tpp = " << hparam.tpp << endl; 
	} else cout << "Hybridization turned off" << endl; 

	// Adjust Slater Condor Parameter
	hparam.SC[0][2] *= 25;
	hparam.SC[1][2] *= 49;
	hparam.SC[1][4] *= 441;
	hparam.SC2EX[2] *= 49;
	hparam.SC2EX[4] *= 441;
	cout << "Ground State Udd: " << (hparam.SC2[0]-hparam.SC2[2]*2/63-hparam.SC2[4]*2/63);
	cout << ", JH: " << (hparam.SC2[2]/14+hparam.SC2[4]/14) << endl;
	cout << "Core-Hole State Udd: " << (hparam.SC2EX[0]-hparam.SC2EX[2]*2/63-hparam.SC2EX[4]*2/63);
	cout << ", JH: " << (hparam.SC2EX[2]/14+hparam.SC2EX[4]/14);
	cout << ", Udp: " << (hparam.FG[0] - hparam.FG[1]/15 - hparam.FG[3]*3/70);
	cout << ", Fdp: " << (hparam.FG[0] + hparam.FG[1]/15 + hparam.FG[3]*3/70) << endl;
	cout << "Diagonalization Option: " << hparam.diag_option << endl;

	// Calculate Delta or Effective Delta
	if (hparam.effective_delta) {
		cout << "effective delta: " << hparam.MLdelta << endl; 
		double del = calculate_effective_delta(IDIR,hparam,pm);
		hparam.MLdelta = 2*hparam.MLdelta - del;
		cout << "calculated delta (used in Hamiltonian): " << hparam.MLdelta << endl; 
	} else {
		cout << "delta (used in Hamiltonian): " << hparam.MLdelta << endl; 
		// double del = calculate_effective_delta(IDIR,hparam,pm);
		// cout << "Calculated effective delta: " << del << endl;
	}

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
	cout << "absorption range: ";
	for (int i = 0; i < 2; ++i) cout << pm.ab_range[i] << ", ";
	cout << endl;
	cout << "emission range: " << pm.em_energy << endl;

	if (!pm.XAS && !pm.RIXS) return 0; // Return if no photon methods
	// Reading and Checking input parameters
	Hilbert GS(IDIR,hparam,pm.edge.substr(0,1),false);
	Hilbert EX(IDIR,hparam,pm.edge.substr(0,1),true);
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Number of Holes: " << GS.num_vh << endl;
	cout << "Run time = " << duration.count() << " ms\n" << endl;
	process_hilbert_space(GS,EX,hparam,pm);

	// Calculate Cross Sections
	vector<double> pvin = pm.pvin, pvout = pm.pvout;
	if (pm.XAS) {
		for (size_t in = 0; in < 3; ++in) {
			fill(pm.pvin.begin(), pm.pvin.end(), 0);
			if (pvin[in]) pm.pvin[in] = 1;
			else continue;
			cout << "Doing photon polarization: ";
			for (auto e : pm.pvin) cout << (int)e << " ";
			cout << endl;
			if (if_photon_file_exist(pm.edge,work_dir,true,overwrite,pm.pvin) && !overwrite) continue;
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
				cout << "Doing photon polarization (in): ";
				for (auto e : pm.pvin) cout << (int)e << " ";
				cout << "(out): ";
				for (auto e : pm.pvout) cout << (int)e << " ";
				cout << endl;
				if (if_photon_file_exist(pm.edge,work_dir,false,overwrite,pm.pvin,pm.pvout) && !overwrite) continue;
				RIXS(GS,EX,pm);
			}
		}
	}

	auto run_stop = chrono::high_resolution_clock::now();
	cout << endl << "Total elapsed time: " << endl;
	cout << ed::format_duration(chrono::duration_cast<chrono::milliseconds>(run_stop-run_start)) << endl;
	return 0;
}