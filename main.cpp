#include <iostream>
#include <fstream>
#include <cstdlib> 
#include <ctime> 
#include <stdlib.h>
#include "diagonalize.h"


using namespace std;
int main(int argc, char** argv){

	// Generate random array
  int n,m;
  double *mat;

  n = 10;
  m = 10;

  mat = new double[n*m];

  srand(time(NULL));
  for (int i = 0; i < n*m; ++i) {
    mat[i] = rand() % 10 + 1;
  }
	cout << "print matrix" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			cout << mat[i*n+j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	diagonalize(mat,n,m);
}