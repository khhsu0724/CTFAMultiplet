# XAS/RIXS Calculation for Fe<sup>3+</sup> Material

## 1. Introduction
This folder provides an example for Fe<sup>3+</sup> (3d<sup>5</sup> configuration) calculation with 3 ligands. The simulated spectra should look similar to measured XAS and RIXS O-_K_ edge and Fe-_L_ edge spectra measured experimentally in Hematite or Li<sub>4</sub>FeSbO<sub>6</sub>. Only Fe-_L_ edge spectra is provided in this example. 

## 2. Running the example
To run the example, compile the program and execute: 
```
$ ../../main
```
To run K edge, simply modify the INPUT file.

## 3. Output files
### 3.1 Miscellaneous files
- `dos.txt`, `exdos.txt`: density of states for ground state and core-hole state Hamiltonian. The output depends of the eigenvalues specified in the INPUT.
- `eig.txt`, `exeig.txt`: Eigenvalues and occupations of the first few eigenvalues for ground state and core-hole state Hamiltonian. Usually the program prints around up to 20 eigenvalues (or less).

### 3.2 XAS Calculation
-  `XAS_Ledge_{pvin}.txt`: simulated XAS. First column is the (unshifted) absorption energy, and the second column is the intensity.

### 3.3 RIXS Calculation
-  `RIXS_Ledge_{pvin}_{pvout}.txt`: simulated RIXS. First column is the (unshifted) absorption energy. Second column is the loss energy and the third column is the intensity.

## 3. References
[chemrxiv preprint](https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/6671eb0e5101a2ffa8e63407/original/a-formal-fe-iii-v-redox-couple-in-an-intercalation-electrode.pdf)
