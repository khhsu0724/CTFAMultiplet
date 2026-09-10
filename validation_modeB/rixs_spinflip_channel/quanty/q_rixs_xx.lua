-- PRX d8 square-planar cluster, HOLE language, mirroring CTFAMultiplet-magnons exactly.
-- modes: 2p core 0-5, 3d 6-15, ligand-x p 16-21, ligand-y p 22-27  (dn even, up odd)
NF = 28; NB = 0
IndexDn_2p = {0,2,4};              IndexUp_2p = {1,3,5}
IndexDn_3d = {6,8,10,12,14};       IndexUp_3d = {7,9,11,13,15}
IndexDn_Lx = {16,18,20};           IndexUp_Lx = {17,19,21}
IndexDn_Ly = {22,24,26};           IndexUp_Ly = {23,25,27}

F0dd,F2dd,F4dd     = 6.0, 6.2999, 11.0250     -- SC2 = 6.0 0 0.12857 0 0.025
F0pp,F2pp          = 1.0, 0.1                 -- SC1 = 1.0 0 0.1  (per ligand, on-site)
F0pd,F2pd,G1pd,G3pd= 6.0, 6.18, 4.635, 2.63   -- FG (F0,G1,F2,G3 reordered)
zeta               = 11.5                     -- SO

Udd  = NewOperator("U", NF, IndexUp_3d, IndexDn_3d, {F0dd,F2dd,F4dd})
Uppx = NewOperator("U", NF, IndexUp_Lx, IndexDn_Lx, {F0pp,F2pp})
Uppy = NewOperator("U", NF, IndexUp_Ly, IndexDn_Ly, {F0pp,F2pp})
Upd  = NewOperator("U", NF, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {F0pd,F2pd}, {G1pd,G3pd})
ld   = NewOperator("ldots", NF, IndexUp_2p, IndexDn_2p)

dofile("h_terms_fix.lua")
Hhop = NewOperator("Number", NF, hop_i, hop_j, hop_v)

Hval = Udd + Uppx + Uppy + Hhop
Hex  = Hval + Upd - zeta*ld

CORE = "111111 0000000000 000000 000000"
VAL  = "000000 1111111111 111111 111111"

t = math.sqrt(1/2)
-- absorption-type operators (hole language: create a 2p hole, destroy a 3d hole)
Tx = NewOperator("CF", NF, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {{1,-1,t},{1,1,-t}})
Ty = NewOperator("CF", NF, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {{1,-1,t*I},{1,1,t*I}})

Psis = Eigensystem(Hval, {NF, NB, {CORE,0,0}, {VAL,2,2}}, 3)
if type(Psis) ~= "table" then Psis = {Psis} end
Egs = Complex.Re(Psis[1]*Hval*Psis[1])
NF_FIN = 231
Chis = Eigensystem(Hval, {NF, NB, {CORE,0,0}, {VAL,2,2}}, NF_FIN)
Phis = Eigensystem(Hex,  {NF, NB, {CORE,1,1}, {VAL,1,1}}, 132)
print(string.format("EGS %.8f", Egs))
for f=1,#Chis do print(string.format("EFIN %d %.8f", f, Complex.Re(Chis[f]*Hval*Chis[f]) - Egs)) end
for n=1,#Phis do print(string.format("EINT %d %.8f", n, Complex.Re(Phis[n]*Hex*Phis[n]) - Egs)) end
-- A[n][i] = <Phi_n| Tx |Psi_i>   (incoming X)
for n=1,#Phis do for i=1,#Psis do
  local a = Phis[n] * Tx * Psis[i]
  print(string.format("A %d %d %.10e %.10e", n, i, Complex.Re(a), Complex.Im(a)))
end end
-- C[n][f] = <Phi_n| Ty |Chi_f>   (outgoing Y);  <Chi_f| Ty^dag |Phi_n> = conj(C[n][f])
for n=1,#Phis do for f=1,#Chis do
  local c = Phis[n] * Tx * Chis[f]
  if Complex.Re(c)^2 + Complex.Im(c)^2 > 1e-24 then
    print(string.format("C %d %d %.10e %.10e", n, f, Complex.Re(c), Complex.Im(c)))
  end
end end
print("DONE")
-- validate the Y-polarization operator against our XAS_Ledge_Y.txt
for n=1,#Phis do
  local w=0
  for i=1,#Psis do local a=Phis[n]*Ty*Psis[i]; w=w+Complex.Re(a)^2+Complex.Im(a)^2 end
  print(string.format("YPOLE %14.8f %18.10e", Complex.Re(Phis[n]*Hex*Phis[n])-Egs, w/#Psis))
end
