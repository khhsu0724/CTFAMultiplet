-- Atomic Ni2+ 3d^8 -> 2p^5 3d^9 L-edge XAS, with the PRX d8 parameter set.
-- No crystal field, no ligands: tests the d^8 initial multiplet, the 2p^5 3d^9 final
-- multiplet, the 2p SOC, and the dipole matrix elements.
NF = 16
NB = 0
IndexDn_2p = {0, 2, 4};              IndexUp_2p = {1, 3, 5}
IndexDn_3d = {6, 8, 10, 12, 14};     IndexUp_3d = {7, 9, 11, 13, 15}

-- effective values this code uses: SC2 = 6.0 0 0.12857 0 0.025 -> F2 = 0.12857*49, F4 = 0.025*441
F0dd, F2dd, F4dd = 0.0, 6.2999, 11.0250
-- FG = F0, G1, F2, G3 = 6.0, 4.635, 6.18, 2.63
F0pd, F2pd, G1pd, G3pd = 0.0, 6.18, 4.635, 2.63
zeta_2p = 11.5

Udd     = NewOperator("U", NF, IndexUp_3d, IndexDn_3d, {F0dd, F2dd, F4dd})
Upd     = NewOperator("U", NF, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {F0pd, F2pd}, {G1pd, G3pd})
ldots2p = NewOperator("ldots", NF, IndexUp_2p, IndexDn_2p)
N2p     = NewOperator("Number", NF, IndexUp_2p, IndexDn_2p, {1,1,1}) +
          NewOperator("Number", NF, IndexDn_2p, IndexUp_2p, {1,1,1})

Hi = Udd
Hf = Udd + Upd + zeta_2p * ldots2p

-- dipole 2p -> 3d
t = math.sqrt(1/2)
Tx = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, {{1,-1,t},{1,1,-t}})
Ty = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, {{1,-1,t*I},{1,1,t*I}})
Tz = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, {{1,0,1}})

-- initial state: 2p^6 3d^8
RestrI = {NF, NB, {"111111 0000000000", 6, 6}, {"000000 1111111111", 8, 8}}
Psis = Eigensystem(Hi, RestrI, 3)
if type(Psis) ~= "table" then Psis = {Psis} end
print("initial-state energies:")
for i = 1, #Psis do print(string.format("  %d  %12.6f", i, Complex.Re(Psis[i]*Hi*Psis[i]))) end

XAS = CreateSpectra(Hf, {Tx, Ty, Tz}, Psis, {{"Emin",-15},{"Emax",35},{"NE",2000},{"Gamma",0.2}})
XAS.Print({{"file","q_d8xas.spec"}})
print("wrote q_d8xas.spec")
