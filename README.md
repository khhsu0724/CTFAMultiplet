# ED
Exact Diagonalization Code

Compile with this command
'''
g++ -pedantic -llapack -std=c++1y -O3 -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp
'''
Remeber to compile with flag -llapack

With OPENMP for macosx
'''
export OMP_NUM_THREADS=NUM_THREADS
g++ -pedantic -llapack -std=c++1y -O3 -Xclang -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp -lomp
'''