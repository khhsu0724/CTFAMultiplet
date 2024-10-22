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
bool read_num(string line, T* tgt, int pnum, bool fix_pnum = true) {
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
		if (fix_pnum && ncnt != pnum) throw invalid_argument("Wrong number of input parameter");
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
		exit(1);
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
							if (p == "SO") skip = read_num(line.substr(s+1,line.size()-1),hparam.SO,2,false);
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
							else if (p == "HFSCALE") skip = read_num(line.substr(s+1,line.size()-1),&hparam.HFscale,1);
							else if (p == "TPD") skip = read_num(line.substr(s+1,line.size()-1),&hparam.tpd,1);
							else if (p == "TPP") skip = read_num(line.substr(s+1,line.size()-1),&hparam.tpp,1);
							else if (p == "TPDZR") skip = read_num(line.substr(s+1,line.size()-1),&hparam.tpdz_ratio,1);
							else if (p == "TPPSIGMA") skip = read_bool(line.substr(s+1,line.size()-1),hparam.tppsigma_on);
							else if (p == "OCTJT") skip = read_num(line.substr(s+1,line.size()-1),&hparam.octJT,1);
							else if (p == "MLCT") skip = read_num(line.substr(s+1,line.size()-1),&hparam.MLdelta,1);
							else if (p == "SIGPI") skip = read_num(line.substr(s+1,line.size()-1),&hparam.sig_pi,1);
							else if (p == "EFFDEL") skip = read_bool(line.substr(s+1,line.size()-1),hparam.effective_delta);
							else if (p == "BLOCK") skip = read_bool(line.substr(s+1,line.size()-1),hparam.block_diag);
							else if (p == "GSDIAG") skip = read_num(line.substr(s+1,line.size()-1),&hparam.gs_diag_option,1);
							else if (p == "EXDIAG") skip = read_num(line.substr(s+1,line.size()-1),&hparam.ex_diag_option,1);
							else if (p == "SITEOCC") skip = read_bool(line.substr(s+1,line.size()-1),hparam.print_site_occ);
							else if (p == "OVERWRITE") skip = read_bool(line.substr(s+1,line.size()-1),overwrite);
							else if (p == "GSNEV") skip = read_num(line.substr(s+1,line.size()-1),&hparam.gs_nev,1);
							else if (p == "EXNEV") skip = read_num(line.substr(s+1,line.size()-1),&hparam.ex_nev,1);
							else if (p == "DIAG") {
								// Generic Diagonalize Option
								skip = read_num(line.substr(s+1,line.size()-1),&hparam.gs_diag_option,1);
								hparam.ex_diag_option = hparam.gs_diag_option;
							}
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
					double incident_temp[3]{0};
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
							// Spectroscopy solver, 1: Exact Solution, 2: Classic KH, 3: 1+2, 4 (To do): Cont FE, Biconj
							else if (p == "SOLVER") skip = read_num(line.substr(s+1,line.size()-1),&pm.spec_solver,1);
							else if (p == "EPSAB") skip = read_num(line.substr(s+1,line.size()-1),&pm.eps_ab,1);
							else if (p == "EPSLOSS") skip = read_num(line.substr(s+1,line.size()-1),&pm.eps_loss,1);
							else if (p == "NITERCFE") skip = read_num(line.substr(s+1,line.size()-1),&pm.niterCFE,1);
							else if (p == "CGTOL") skip = read_num(line.substr(s+1,line.size()-1),&pm.CG_tol,1);
							else if (p == "NEDOS") skip = read_num(line.substr(s+1,line.size()-1),&pm.nedos,1);
							else if (p == "ABMAX") skip = read_num(line.substr(s+1,line.size()-1),&pm.abmax,1);
							else if (p == "SPINFLIP") skip = read_bool(line.substr(s+1,line.size()-1),pm.spin_flip);
							else if (p == "AB") {
								skip = read_num(line.substr(s+1,line.size()-1),abrange_temp,2);
								for (int i = 0; i < 2; ++i) pm.ab_range[i] = abrange_temp[i];
							}
							else if (p == "INCIDENT") {
								skip = read_num(line.substr(s+1,line.size()-1),incident_temp,3);
								for (int i = 0; i < 3; ++i) pm.incident[i] = incident_temp[i];
								pm.set_incident_points();
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
		exit(1);
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
	// Measures the energy difference between dn/dn+1 if Delta is set to 0
	double inp_tpd = hparam.tpd, inp_tpp = hparam.tpp, mlct = hparam.MLdelta;
	int diag_option = hparam.gs_diag_option;
	hparam.gs_diag_option = 4;
	hparam.tpd = 0;
	hparam.tpp = 0;
	Hilbert GS(input_dir,hparam,pm.edge,false);
	GS.HYB_on = true;
	double U_guess = hparam.SC[1][0]*(ed::choose(GS.num_vh,2)-ed::choose(GS.num_vh-1,2));
	hparam.MLdelta = U_guess;
	if (!GS.cluster->lig_per_site) {
		cout << "No ligands, skipping calculate delta" << endl;
		return 0;
	}
	// Redo delta here, should not need to reset the values
	calc_ham(GS,hparam);
	// for (size_t g = 0; g < (GS.hblks.size()+1)/2; g++) {
	for (size_t g = 0; g < GS.hblks.size(); g++) {
		GS.hblks[g].diagonalize(100);
	}
	hparam.tpd = inp_tpd;
	hparam.tpp = inp_tpp;
	hparam.MLdelta = mlct;
	hparam.gs_diag_option = diag_option;
	return effective_delta(GS,2) - U_guess;
}


void dos_eigenenergy(Hilbert& hilbs, int num_eig, string dos_fname = "", string eig_fname = "") {
	multistream dos_out(false,dos_fname,"w"); // Make sure that it clears the files
	multistream eig_out(false,eig_fname,"w"); 
	dos_out.set_mode("a");
	eig_out.set_mode("a");
	auto all_eig = hilbs.get_all_eigval(true);
	auto unique_eig = ed::printDistinct(all_eig,0.0,all_eig.size(),false);
	double gs_en = unique_eig[0];
	for (int u = 0; u < unique_eig.size(); u++) {
		vector<bindex> eigind;
		for (auto & hblks : hilbs.hblks) {
			for (size_t i = 0; i < hblks.nev; ++i) {
				if (abs(hblks.eig[i]-unique_eig[u]) < TOL) {
					eigind.push_back(bindex(&hblks-&hilbs.hblks[0],i));
				}
			}
		}
		dos_out << unique_eig[u] << " " << eigind.size() << endl;
		if (u < num_eig) {
			eig_out << "Eigen-energy number: " << u+1;
			eig_out << ", energy (rel to gs): " << unique_eig[u] - gs_en << endl; 
			eig_out << "Spin: " << abs(hilbs.hblks[eigind[0].first].get_sz()); 
			eig_out << ", degeneracy: " << eigind.size() << endl;
			eig_out << "------------------------------" << endl; 
			occupation(hilbs,eigind,false,eig_fname,false);
			eig_out << "------------------------------" << endl; 
		}
	}
	return;
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

	bool clear_mat = !(pm.spec_solver == 4);
	for (size_t g = 0; g < GS.hblks.size(); g++) {
		start = chrono::high_resolution_clock::now();
		cout << "Diagonalizing block number: " << g << ", matrix size: " << GS.hblks[g].size << endl;
		GS.hblks[g].diagonalize(hparam.gs_nev,clear_mat);
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
	if (GS.BLOCK_DIAG) cout << "Grounds State Spin Quantum Number (S): " << SDegen << endl;
	cout << "Calculating Occupation (with degeneracy): " << gsi.size() << endl;
	occupation(GS,gsi,true,"");
	occupation(GS,gsi,true,"",true);
	wvfnc_weight(GS,gsi,3,true);
	// auto all_eig = GS.get_all_eigval(true);
	// ed::printDistinct(all_eig,0.0,all_eig.size(),true);
	dos_eigenenergy(GS,20,"dos.txt","eig.txt");
	cout << endl;

	// Assemble Core Hole Hamiltonian
	cout << "Assembling Core Hole Hamiltonian..." << endl;
	start = chrono::high_resolution_clock::now();
	hparam.SC.pop_back();
	hparam.SC.push_back(hparam.SC2EX);
	calc_ham(EX,hparam);
	stop = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
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
		cout << "Diagonalizing block number: " << e << ", matrix size: " << EX.hblks[e].size << endl;
		// Only need the lowest eigenvalue if using Lanczos
		if (hparam.ex_nev > 0) {
			if (pm.spec_solver == 4) EX.hblks[e].diagonalize(5,clear_mat);
			else EX.hblks[e].diagonalize(hparam.ex_nev,clear_mat);
		} else {
			if (pm.spec_solver != 4) {
				cout << "WARNING: EXNEV = 0" << endl;
				exit(1);
			} else cout << "skipping core-hole state diagonalization" << endl;
		}
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
	if (pm.incident[2] == -1) pm.set_incident_points(true,gs_en,ex_min_en);
	dos_eigenenergy(EX,10,"exdos.txt","exeig.txt");
	wvfnc_weight(EX,exmin_ind,3,true);
	// all_eig = EX.get_all_eigval(true);
	// ed::printDistinct(all_eig,0.0,all_eig.size(),true);

	if (EX.BLOCK_DIAG) cout << "Core-Hole Ground State Spin Quantum Number (S): " << SDegen << endl;
	cout << "Calculating Occupation (with degeneracy): " << exmin_ind.size() << endl;
	occupation(EX,exmin_ind,true);
	wvfnc_weight(EX,exmin_ind,3);
	cout << endl << "Total diagonalize run time: " << ed::format_duration(diag_duration) << endl << endl;

	if (pm.abmax != -1) {
		pm.ab_range[0] = ex_min_en - gs_en - 4;
		pm.ab_range[1] = pm.ab_range[0] + pm.abmax;
		cout << "======> NEW AB range: " << pm.ab_range[0] << ", " << pm.ab_range[1] << endl;
	}

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
	} else if (pm.edge != "L" && pm.edge != "L3") throw invalid_argument("invalid edge input: " + pm.edge);
	// U = F^0 + 4*F^2 + 36*F^4
	// Scaling Slater Parameters, without the bare Coulomb term
	for (int i = 1; i < 3; ++i) hparam.SC[0][i] *= hparam.HFscale;
	for (int i = 1; i < 5; ++i) hparam.SC[1][i] *= hparam.HFscale;
	for (int i = 1; i < 5; ++i) hparam.SC2EX[i] *= hparam.HFscale;
	for (int i = 1; i < 4; ++i) hparam.FG[i] *= hparam.HFscale;

	Hilbert GS(IDIR,hparam,pm.edge.substr(0,1),false);
	Hilbert EX(IDIR,hparam,pm.edge.substr(0,1),true);

	cout << "Input parameters" << endl;
	cout << "TM 2p Spin-Orbit Coupling: " << hparam.SO[0] << " eV" << endl;
	cout << "TM 3d Spin-Orbit Coupling: " << hparam.SO[1] << " eV" << endl;
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
		cout << ", tpp = " << hparam.tpp; 
		if (GS.coord == "oct") cout << ", z distortion = " << hparam.octJT << endl; 
		if (GS.coord != "ion") {
			cout << ", sigma pi bond ratio = " << hparam.sig_pi; 
			if (hparam.tppsigma_on) cout << ", tppsigma on (only tppzpi!!!)" << endl; 
			else cout << ", tppsigma off" << endl; 
		}
		else cout << endl;
	} else cout << "Hybridization turned off" << endl; 

	// Adjust Slater Condor Parameter

	cout << "J1 (xz/yz-x2/xy): " << 3*hparam.SC2[2]+20*hparam.SC2[4] << endl;
	cout << "J2 (z2-x2/xy): " << 4*hparam.SC2[2]+15*hparam.SC2[4] << endl;
	cout << "J3 (x2-xy): " << 35*hparam.SC2[4] << endl;
	cout << "J4 (z2-xz/yz): " << hparam.SC2[2]+30*hparam.SC2[4] << endl;

	cout << "Racah Parameter, A: " << (hparam.SC2[0]-49*hparam.SC2[4]) << ", B: " << (hparam.SC2[2]-5*hparam.SC2[4]) << ", C: " << (35*hparam.SC2[4]) << endl;
	hparam.SC[0][2] *= 25;
	hparam.SC[1][2] *= 49;
	hparam.SC[1][4] *= 441;
	hparam.SC2EX[2] *= 49;
	hparam.SC2EX[4] *= 441;
	cout << "Ground State Uavg: " << (hparam.SC2[0]+hparam.SC2[2]*2/63+hparam.SC2[4]*2/63);
	cout << ", JH : " << (hparam.SC2[2]*2.5/49+hparam.SC2[4]*36/441) << " (Kanamori)"; // Adjusted for real space Harmonics
	// cout << "Core-Hole State Udd: " << (hparam.SC2EX[0]-hparam.SC2EX[2]*2/63-hparam.SC2EX[4]*2/63);
	// cout << ", JH: " << (hparam.SC2EX[2]/14+hparam.SC2EX[4]/14*5/7);
	cout << ", Udp: " << (hparam.FG[0] - hparam.FG[1]/15 - hparam.FG[3]*3/70);
	cout << ", Fdp: " << (hparam.FG[0] + hparam.FG[1]/15 + hparam.FG[3]*3/70) << endl;
	cout << "GS Diagonalization Option: " << hparam.gs_diag_option << endl;
	cout << "EX Diagonalization Option: " << hparam.ex_diag_option << endl;

	// Calculate Delta or Effective Delta
	if (hparam.effective_delta) {
		cout << "effective delta: " << hparam.MLdelta << endl; 
		double del = calculate_effective_delta(IDIR,hparam,pm);
		hparam.MLdelta = hparam.MLdelta - del;
		cout << "calculated delta (used in Hamiltonian): " << hparam.MLdelta << endl; 
	} else {
		cout << "delta (used in Hamiltonian): " << hparam.MLdelta << endl; 
		double del = calculate_effective_delta(IDIR,hparam,pm);
		cout << "Calculated effective delta: " << del + hparam.MLdelta << endl;
	}
	cout << "solver options: " <<  pm.spec_solver << endl;
	if (pm.spec_solver == 4) {
		cout << "Lanczos/BiCGS Methods" << endl;
		cout << "absorption broadening: " << pm.eps_ab << endl;
		cout << "loss broadening: " << pm.eps_loss << endl;
		cout << "conjugate gradient tolerance: " << pm.CG_tol << endl;
		cout << "CFE iterations: " << pm.niterCFE << endl;
	}
	else if (pm.RIXS) {
		if (pm.eloss) cout << "	Using energy loss" << endl;
		else cout << "	WARNING: All spectroscopy will be outputted as energy loss" << endl;
		pm.eloss = true;
		if (pm.spec_solver == 1 || pm.spec_solver == 3)
			cout << "	Output exact RIXS peaks" << endl;
		if (pm.spec_solver == 2 || pm.spec_solver == 3)
			cout << "	Evaluate K-H equation, Gamma: " << pm.eps_loss << endl;
	}
	cout << "photon edge: " << pm.edge << endl;
	cout << "pvin: ";
	for (int i = 0; i < 3; ++i) cout << pm.pvin[i] << ", ";
	cout << endl << "pvout: ";
	for (int i = 0; i < 3; ++i) cout << pm.pvout[i] << ", ";
	cout << endl;
	if (pm.abmax != -1) {
		cout << "ABMAX (overriding AB): ";
		cout << pm.abmax << endl;
	} else {
		cout << "absorption range: ";
		for (int i = 0; i < 2; ++i) cout << pm.ab_range[i] << ", ";
		cout << endl;
	}
	cout << "Loss energy (for RIXS): " << pm.em_energy << endl;

	cout << "Incident (For solver = 4): ";
	for (int i = 0; i < 3; ++i) cout << pm.incident[i] << ", ";
	cout << endl;
	// Reading and Checking input parameters
	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Number of Holes: " << GS.num_vh << endl;
	cout << "Run time = " << duration.count() << " ms\n" << endl;
	process_hilbert_space(GS,EX,hparam,pm);
	// Calculate Cross Sections
	vector<double> pvin = pm.pvin, pvout = pm.pvout;

	cout << "RIXS Energy points: ";
	if (pm.inc_e_points.size() != 0)
		for (auto &inc_e : pm.inc_e_points) cout << inc_e << ", ";
	cout << endl;

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
