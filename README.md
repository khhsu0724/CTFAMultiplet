# ED
Exact Diagonalization Code
Different compilation method provided below

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

With OPENMP + MKL for macos
```console
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/latest/mac/compiler/lib:$LD_LIBRARY_PATH
g++ -pedantic -L${MKLROOT}/lib -lmkl_scalapack_lp64 -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_mpich_lp64 -L/opt/intel/oneapi/compiler/latest/mac/compiler/lib/ -liomp5 -lpthread -lm -ldl -llapack -std=c++1y -O3 -Xclang -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp -lomp
```

Compile it in NERSC Haswell (-pedantic will give a lot of warnings)
```console
CC -std=c++1y -mkl="parallel" -O3 -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp
```

Compile it in Sherlock (needs gcc 6.3.0 to avoid compile errors)
```console
module load imkl icc boost gcc
icpc -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl -std=c++1y -O3 -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp -lm
```