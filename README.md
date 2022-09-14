# ED
Exact Diagonalization Code

Compile with this command on MAC
```console
g++ -pedantic -llapack -std=c++1y -O3 -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp
```
Remeber to compile with flag -llapack

With OPENMP for macosx
```console
export OMP_NUM_THREADS=NUM_THREADS
g++ -pedantic -llapack -std=c++1y -O3 -Xclang -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp -lomp
```

Compile it in NERSC (-pedantic will give a lot of warnings which I will fix some day)
```console
g++ -llapack -std=c++1y -O3 -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp
```