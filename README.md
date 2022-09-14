# ED
Exact Diagonalization Code
Different compilation method provided below, the code is precompiled for MacOS and can be run through INPUT

1. Compile with this command on MAC
```console
g++ -pedantic -llapack -std=c++1y -O3 -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp
```
Remeber to compile with flag -llapack

2. With OPENMP for macosx
```console
export OMP_NUM_THREADS=NUM_THREADS
g++ -pedantic -llapack -std=c++1y -O3 -Xclang -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp -lomp
```

3. Compile it in NERSC (-pedantic will give a lot of warnings which I will fix some day)
```console
g++ -llapack -std=c++1y -O3 -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp
```
