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

