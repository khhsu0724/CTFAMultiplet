#include "phonon.hpp"
using namespace std;

void PhononModes::print_phonons_occuptation(const vecd& occ, bool is_print, string fname) {
	multistream mout(is_print,fname);
	// Right now only accept fixed length array
	mout << "***************************" << endl;
	if (occ.size() != get_unique_modes()) throw runtime_error("Wrong phonon occ size");
	for (size_t i = 0; i < get_unique_modes(); i++) {
		mout << "Phonon mode " << i << ": " << occ[i];
	}
	mout << endl;
	return;
};

// should look something like this
/*
&PHONON
    PN = 10 0.3 # 10 modes, omega = 0.3 following line are modes
    h 0.22 5 # h: holestein mode, g_h, orbital index
    b 0.13 13,22 # b: breathing mode, g_b, pairs of 2 orbitals 
    # Takes the first frequency value per mode
    # you can have has any number of coupling in a mode
    PN = 5 0.1 # Starts next phonon by PN

/
*/
void PhononModes::read_phonon_from_input(string file_dir) {
	string line;
	ifstream input(file_dir);
	try {
		if (input.is_open()) {
			bool read_phonon = false;
			Phonon_Param pntemp;
			int num_phonons;
			while (getline(input,line)) {
				if (line[0] == '#') continue;
				if (line[0] == '/') {
					read_phonon = false;
					continue;
				}
				if (read_phonon) {
					string p;
					bool skip = false;
					double readin[5]{0}; // Pretty random length
					for (int s = 0; s < line.size() && !skip; ++s) {
						if (line[s] != ' ' && line[s] != '	' && line[s] != '=') p.push_back(line[s]);
						else {
							transform(p.begin(),p.end(),p.begin(),::toupper);
							if (p == "PN") {
								// Push the previous phonon mode into Phonons
								if (pntemp.is_set) {
									this->addPhonon(num_phonons,pntemp);
									pntemp.clear();
								}
								// Read now
								pntemp.is_set = true;
								skip = ed::read_num(line.substr(s+1,line.size()-1),readin,2,p=p);
								pntemp.omega = readin[1];
								num_phonons = int(readin[0]);
							}
							else if (p == "H") {
								// cout << line.substr(s+1,line.size()-1) << endl;
								auto outline = ed::read_arb(line.substr(s+1,line.size()-1));
								if (!pntemp.g_h) pntemp.g_h = stod(outline[0]);
								else cout << "WARNING: g_h already set to: " << pntemp.g_h <<  "." << endl;
								for (int i = 1; i < outline.size(); ++i) {
									pntemp.h_orb_i.push_back(stoi(outline[i]));
								}
							}
							else if (p == "B") {
								// cout << line.substr(s+1,line.size()-1) << endl;
								auto outline = ed::read_arb(line.substr(s+1,line.size()-1));
								if (!pntemp.g_b) pntemp.g_b  = stod(outline[0]);
								else cout << "WARNING: g_b already set to: " << pntemp.g_b <<  "." << endl;
								for (int i = 1; i < outline.size(); ++i) {
									size_t commaPos = outline[i].find(',');
								    if (commaPos == std::string::npos) {
								    	throw std::invalid_argument("Invalid format for breathing mode: two integers separated by a comma");
								    }
								    int first = std::stoi(outline[i].substr(0, commaPos));
								    int second = std::stoi(outline[i].substr(commaPos + 1));
								    pntemp.b_orb_ij.push_back(std::pair<int,int>(first,second));
								}
							}
							p = "";
						}
						if (line[s] == '#') skip = true;
					}
				}
				if (line == "&PHONON") read_phonon = true;
			}
			// Push back the latest phonon
			if (pntemp.is_set) {
				this->addPhonon(num_phonons,pntemp);
				pntemp.clear();
			}
		} else throw invalid_argument("Cannot open " + file_dir + " file");
	} catch (const exception &ex) {
		cerr << ex.what() << "\n";
		exit(1);
	}
	return;
};

void PhononModes::print_phonon_data() {
	cout << "-----PHONON read from input-----" << endl;
	cout << "Total number of phonon modes: " << this->get_unique_modes() << endl;
	int count = 1;
	for (auto & pn : phonons) {
		cout << "Phonon number: " << count << ", omega: " << pn->pnParam.omega << endl;
		if (pn->is_holestein()) {
			cout << "Holestein g: " << pn->pnParam.g_h;
			cout << ", on orbitals: ";
			for (auto & orb : pn->pnParam.h_orb_i) {
				cout << this->cluster->get_orb_name(orb) << "(" << orb << ")" << ", ";
			}
			cout << endl;
		}
		if (pn->is_breathing()) {
			cout << "Breathing g: " << pn->pnParam.g_b;
			cout << ", on orbitals: ";
			for (auto & orbs : pn->pnParam.b_orb_ij) {
				cout << this->cluster->get_orb_name(orbs.first) << ".";
				cout << this->cluster->get_orb_name(orbs.second);
				cout << "(" << orbs.first << "." << orbs.second << ")" << ", ";
			}
			cout << endl;
		}
		count++;
	}
	cout << "--------------------------------" << endl;
};

bindex PhononModes::norm_Hash(const veci &s) {
	// No phonons conservation
	int ind = 0;
	for (int i = 0; i < this->get_unique_modes(); i++)
		ind += s[i] * this->strides[i];
	return bindex(0,ind);
};

veci PhononModes::norm_Hashback(bindex ind) {
	int index = ind.second;
	veci state(get_unique_modes(),0);
	for (int m = 0; m < get_unique_modes(); m++) {
		state[m] = index/strides[m];
		index = index % strides[m];
	}
	return state;
};

