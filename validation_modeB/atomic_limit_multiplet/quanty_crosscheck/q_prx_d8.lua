-- Atomic 2p^5 3d^9 multiplet with ONLY the 2p-3d Coulomb interaction.
-- No SOC, no crystal field, no hybridization -- the same limit used in
-- validation_modeB/atomic_limit_multiplet/.
NF = 16
NB = 0
IndexDn_2p = {0, 2, 4}
IndexUp_2p = {1, 3, 5}
IndexDn_3d = {6, 8, 10, 12, 14}
IndexUp_3d = {7, 9, 11, 13, 15}

-- our NiO values after HFscale = 0.8 (unreduced Slater integrals, same convention)
F0 = 6.0
F2 = 6.18
G1 = 4.635
G3 = 2.63

Udp = NewOperator("U", NF, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {F0, F2}, {G1, G3})
H = Udp

-- restrict to exactly 5 electrons in 2p and 9 in 3d  ->  60 states
Restr = {NF, NB, {"111111 0000000000", 5, 5}, {"000000 1111111111", 9, 9}}
Psis = Eigensystem(H, Restr, 60)

print("QUANTY 2p^5 3d^9 eigenvalues (F0=" .. F0 .. " F2=" .. F2 .. " G1=" .. G1 .. " G3=" .. G3 .. ")")
for i = 1, #Psis do
  print(string.format("%3d  %14.6f", i, Complex.Re(Psis[i] * H * Psis[i])))
end
