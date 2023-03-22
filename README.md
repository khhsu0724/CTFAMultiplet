# ED
Exact Diagonalization Code
Different compilation method provided below
To see stack trace error, compile with -rdynamic
Arpack & MKL installation is recommended, arpack will use blas if MKL is not present

1. Compile with this command on MAC
```console
g++ -pedantic -llapack -larpack -lblas -std=c++1y -O3 -ffast-math -march=native -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp
```
Remeber to compile with flag -llapack

2. With OPENMP for macosx
```console
export OMP_NUM_THREADS=NUM_THREADS
g++ -pedantic -llapack -larpack -std=c++1y -O3 -Xclang -fopenmp -ffast-math -march=native -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp -lomp
```

With OPENMP + MKL for macos (remove -pedantic & -Xclang when compiling on Linux machine)
```console
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/latest/mac/compiler/lib:$LD_LIBRARY_PATH
g++ -pedantic -L${MKLROOT}/lib -lmkl_scalapack_lp64 -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_mpich_lp64 -L/opt/intel/oneapi/compiler/latest/mac/compiler/lib/ -liomp5 -lpthread -lm -ldl -llapack -larpack -std=c++1y -O3 -Xclang -fopenmp -ffast-math -march=native -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp -lomp
```

With OPENMP + MKL for Linux Machine
```console
g++ -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -llapack -larpack -std=c++1y -O3 -fopenmp -ffast-math -march=native -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp
```

I have not figured out how to compile arpack on NERSC & Sherlock
Compile it in NERSC Haswell (-pedantic will give a lot of warnings)
```console
CC -std=c++1y -mkl="parallel" -O3 -fopenmp -ffast-math -march=native -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp
```

Compile it in Sherlock
Compile with MKL
```console
module load imkl icc boost gcc/6.3.0
icpc -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl -std=c++1y -O3 -fopenmp -ffast-math -march=native -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp -lm
```
Compile with Arpack
```console
module load arpack icc boost gcc/8.1.0
icpc -llapack -lpthread -ldl -larpack -std=c++1y -m64 -O3 -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp -lm -lgfortran
```