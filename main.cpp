#include <cstdlib> 
#include <ctime> 
#include <stdlib.h>
#include <bitset>
#include "multiplet.hpp"
#include "photon.hpp"

using namespace std;

void read_input() {
	return;
}

void tb() {
	ofstream myfile;
    myfile.open ("./tb.txt");
	double SC[5] = {0,0,1*49,0,0.078*441};
	double racah_B = (SC[2]/49) - (5*SC[4]/441);
	for (double cf = 0.5; cf <= 50; cf += 0.5) {
		Hilbert input("./INPUT",false);
		calc_coulomb(input,SC);
		calc_CF(input,0,cf,6);
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

int main(int argc, char** argv){
	// Process flags

	bool CF_on = false;
	bool SO_on = false;

	string file_dir = "./INPUT";
	string photon_method;
	Hilbert input(file_dir,true);

	// double SC[5] = {1.0,0,1.0*49,0,1.0*441};
	double SC[5] = {6.0,0,0.13*49,0,0.025*441}; //F_0,2,4 = 6, 0.13, 0.025 eV
	// double FG[4] = {6.0,6.2,6.18,4.63}; // Core-Valence Interaction Parameter
	double FG[4] = {6.0,6.2,6.18,4.63};
	double racah_B = (SC[2]/49) - (5*SC[4]/441);
	int nd = 0;
	for (auto &at:input.atlist) if(at.is_val && at.l == 2) nd = at.num_h - input.is_ex;
	double del = 3+(nd-1)*(SC[0]-SC[2]*2.0/63-SC[4]*2.0/63);
	cout << "nd: " << nd << ", del: " << del << endl;

	vector<double> pvec = {1,1,1};
	calc_coulomb(input,SC);
	// calc_SO(input,10.5); 
	// calc_CF(input,del,1,6); 
	calc_CV(input,FG);
	// ed::write_mat(input.hblks[0].ham,input.hblks[0].size,input.hblks[0].size,"./ham.txt");
	input.hblks[0].diag_dsyev();
	// cout << "eigenvalues: ";
	// ed::printDistinct(input.hblks[0].eig,input.hblks[0].size);
	// cout << endl;

	cout << "sorted eigenvalues: " << endl;
	sort(input.hblks[0].eig,input.hblks[0].eig+input.hblks[0].size);
	double uniq_eig = input.hblks[0].eig[0],  count = 1;
	for (int i = 1; i < input.hblks[0].size; i++) {
		if (abs(uniq_eig-input.hblks[0].eig[i]) < 1e-7) count++;
		else {
			cout << uniq_eig << ": " << count << ", ";
			count = 1;
			uniq_eig = input.hblks[0].eig[i];
		}
	}
	cout << uniq_eig << ": " << count;
	cout << endl;


	XAS(SC,FG,0,10.5,pvec,100);
	
	// delete [] mat;
}