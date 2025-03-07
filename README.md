# Charge Transfer Full Atomic Multiplet

This project provides a **charge transfer full atomic multiplet code** designed to calculate **X-ray Absorption Spectroscopy (XAS)** and **Resonant Inelastic X-ray Scattering (RIXS)** spectra, for both transition metal *L* edge and ligand *K* edge (dipole transitions). A description of the Hamiltonian and results can be found in this [Chemrxiv](https://chemrxiv.org/engage/chemrxiv/article-details/6671eb0e5101a2ffa8e63407).

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
## Package Dependencies
- **C++ Compiler**
- **Boost**
- **OpenMP** (Optional)
- **BLAS/LAPACK**
- **MKL** (Optional)
- **ARPACK** 
- **cmake** (Recommended)

[Arpack](https://github.com/opencollab/arpack-ng) is required for this code, enable iso_c_binding when installing. \
Often, using [homebrew](https://brew.sh/) to install arpack is sufficient. \
Before compiling this code, please make sure all required dependencies are properly installed.

## System and software Requirements
The code is designed for Unix based machine with enough RAM to support the INPUT defined by a user. \
For large matrix problem, it is recommended to use a workstation or computing node with large memory, for example, [NERSC cpu node](https://docs.nersc.gov/systems/perlmutter/architecture/)
The only software version requirements for this program are:
```
c++14
cmake 3.15
```

## Compile and running
If required libraries are not installed, Linux and Mac users can install them via *brew*:
```
$ brew install boost libomp arpack
```

To install the code (with cmake):
```
$ git clone https://github.com/khhsu0724/CTFAMultiplet.git
$ cd CTFAMultiplet/build
$ cmake ..
$ make
$ ../main ./test/INPUT # To execute the code
```

For custom installed Arpack, cmake requires path to the installation directory using the *ARPACK_ROOT* variable:
```
$ cmake -DARPACK_ROOT=/path/to/arpack/ ..
```

Note that for sherlock, icpc compiler is needed instead of c++. \
Once in the build directory, issue these command: 
```
$ module load arpack icc boost gcc/8.1.0
$ cmake -DCMAKE_CXX_COMPILER=icpc ..
$ make
```

There are also Makefile examples for: Linux, macOS, sherlock, perlmutter (NERSC) \
The installation process takes 20 seconds.
Please contact the author if you run into issues compiling,\
For more compiling information see this [wiki](https://github.com/khhsu0724/CTFAMultiplet/wiki/Getting-Started) page


## INPUT cards
The code reads in an *INPUT* file with any name that can be specified in the first argument after the executable:
```
$ ./main INPUT_FILE_NAME
```
The default file name is INPUT if none are specified.\
For detailed input cards, please refer to the **[Wiki](https://github.com/khhsu0724/CTFAMultiplet/wiki/Input-Parameters)**.
