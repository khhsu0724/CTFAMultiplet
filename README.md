# ED
Exact Diagonalization Code

Compile with this command
'''
g++ -pedantic -llapack -std=c++1y -O3 -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp
'''
Remeber to compile with flag -llapack