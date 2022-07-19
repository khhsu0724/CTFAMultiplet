#include <cstdlib> 
#include <ctime> 
#include <stdlib.h>
#include <bitset>
#include <regex>
#include "multiplet.hpp"
#include "photon.hpp"

using namespace std;

struct PM {
	bool XAS = false;
	bool RIXS = false;
	bool PE = false; // Photo Emission
	string edge;
};

void tb() {
	ofstream myfile;
    myfile.open ("./tb.txt");
	double SC2[5] = {0,0,1*49,0,0.078*441};
	double SC1[3] = {0,0,0};
	vector<double*> SC;
	SC.emplace_back(SC1);
	SC.emplace_back(SC2);
	double CF[5] = {1,1,-2.0/3,-2.0/3,-2.0/3};
	double FG[4]{0};
	for (double m = 0.5; m <= 50; m += 0.5) {
		Hilbert input("./INPUT",SC,FG,CF,0,false,false);
		double CF_TB[5]{0};
		for (int c = 0; c < 5; c++) CF_TB[c] = CF[c] * m;
		calc_coulomb(input,SC);
		calc_CF(input,0,CF_TB);
		input.hblks[0].diag_dsyev();
		vector<double> unique_eig = ed::printDistinct(input.hblks[0].eig,input.hblks[0].size,false);
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

void read_input(string IDIR, PM& pm, vector<double*> SC, double* FG, 
				double* CF, double& SO, bool& HYB, double& nedos, 
				vecd& pvin, vecd& pvout) {
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
							if (p == "SO") skip = read_num(line.substr(s+1,line.size()-1),&SO,1);
							else if (p == "SC1") skip = read_num(line.substr(s+1,line.size()-1),SC[0],3);
							else if (p == "SC2") skip = read_num(line.substr(s+1,line.size()-1),SC[1],5);
							else if (p == "FG") skip = read_num(line.substr(s+1,line.size()-1),FG,4);
							else if (p == "CF") skip = read_num(line.substr(s+1,line.size()-1),CF,5);
							else if (p == "HYB") skip = read_bool(line.substr(s+1,line.size()-1),HYB);
							else if (p == "NEDOS") skip = read_num(line.substr(s+1,line.size()-1),&nedos,1);
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
								for (int i = 0; i < 3; ++i) pvin[i] = pvtmp[i];
							}
							else if (p == "PVOUT") {
								skip = read_num(line.substr(s+1,line.size()-1),pvtmp,3);
								for (int i = 0; i < 3; ++i) pvout[i] = pvtmp[i];
							}
							else if (p == "EDGE") {
								cout << p << endl;
								// TODO: if FG size does not match edge then throw error
								// if edge => K, FG size = 2
								// if edge => L, FG size = 4
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
	double SO = 0, nedos = 0;
	// Change SC Parameter here to l = 1 & l = 2
	double SC2[5]{0}, SC1[3]{0}, FG[4]{0}, CF[5]{0};
	vecd HYB_param; // Use vector since the size  is unknown
	bool HYB;
	vector<double*> SC;
	SC.emplace_back(SC1);
	SC.emplace_back(SC2);
	vecd pvin(3,0), pvout(3,0);

	string IDIR = "./INPUT";
	if (argc == 2) IDIR = string(argv[1]);
	read_input(IDIR,pm,SC,FG,CF,SO,HYB,nedos,pvin,pvout);

	// U = F^0 + 4*F^2 + 36*F^4

	cout << "Input parameters" << endl;
	cout << "SO: " << SO << "eV" << endl;
	cout << "SC1: ";
	for (int i = 0; i < 3; ++i) cout << SC[0][i] << ", ";
	cout << endl << "SC2: ";
	for (int i = 0; i < 5; ++i) cout << SC[1][i] << ", ";
	cout << endl << "FG: ";
	for (int i = 0; i < 4; ++i) cout << FG[i] << ", ";
	cout << endl << "CF: ";
	for (int i = 0; i < 5; ++i) cout << CF[i] << ", ";
	cout << endl;
	if (HYB) cout << "yes HYB" << endl; 
	else cout << "no HYB" << endl; 
	cout << "pvin: ";
	for (int i = 0; i < 3; ++i) cout << pvin[i] << ", ";
	cout << endl << "pvout: ";
	for (int i = 0; i < 3; ++i) cout << pvout[i] << ", ";
	cout << endl;
	
	// Adjust Slater Condor Parameter
	SC[0][2] *= 49;
	SC[1][2] *= 49;
	SC[1][4] *= 441;
	if (pm.XAS) XAS(IDIR,SC,FG,CF,SO,HYB,pvin,nedos,pm.edge);
	if (pm.RIXS) RIXS(IDIR,SC,FG,CF,SO,HYB,pvin,pvout,nedos,pm.edge);	

	return 0;
}