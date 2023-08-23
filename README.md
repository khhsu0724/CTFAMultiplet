# Charge Transfer Full Atomic Multiplet
Charge transfer full atomic multiplet code to calculate XAS & RIXS spectroscopy (Dipole allowed transitions)  
Using Exact Diagonalization or Lanczos, depending on matrix size. 
TODO: Input descriptions  
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

Compile it in NERSC Cori (-pedantic will give a lot of warnings)
```console
CC -std=c++1y -mkl="parallel" -O3 -fopenmp -ffast-math -march=native -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp
```
Using arpack in Cori is a bit complicated, you will need to install your own arpack with ico_c_binding and link it to the compiler, for example:
```console
module load cpu # For Perlmutter 
-mkl="parallel" ==> flag if you want to use mkl
CC -I${INSTALL_DIR}/include -L${INSTALL_DIR}/lib64 -larpack -std=c++1y -O3 -fopenmp -ffast-math -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp
export LD_LIBRARY_PATH=${INSTALL_DIR}/lib64:$LD_LIBRARY_PATH
```

Compile it in Sherlock with MKL
```console
module load imkl icc boost gcc/6.3.0
icpc -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl -std=c++1y -O3 -fopenmp -ffast-math -march=native -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp -lm
```
Compile it in Sherlock with Arpack
```console
module load arpack icc boost gcc/8.1.0
icpc -llapack -lpthread -ldl -larpack -std=c++1y -m64 -O3 -fopenmp -o main main.cpp diagonalize.cpp gaunt.cpp hilbert.cpp multiplet.cpp photon.cpp helper.cpp cluster.cpp -lm -lgfortran
```
