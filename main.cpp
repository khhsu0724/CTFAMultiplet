#include <cstdlib> 
#include <ctime> 
#include <stdlib.h>
#include <bitset>
#include <regex>
#include "multiplet.hpp"
#include "photon.hpp"

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

void read_input(string IDIR, PM& pm, HParam& hparam) {
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
							else if (p == "NEDOS") skip = read_num(line.substr(s+1,line.size()-1),&hparam.nedos,1);
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
		} else throw invalid_argument("Cannot open INPUT file");
	} catch (const exception &ex) {
		cerr << ex.what() << "\n";
		exit(0);
	}
	return;
}

int main(int argc, char** argv){
	PM pm;
	HParam hparam; 
	string IDIR = "./INPUT";
	if (argc == 2) IDIR = string(argv[1]);
	read_input(IDIR,pm,hparam);
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

	if (pm.XAS) XAS(IDIR,pm,hparam);
	if (pm.RIXS) RIXS(IDIR,pm,hparam);	

	return 0;
}