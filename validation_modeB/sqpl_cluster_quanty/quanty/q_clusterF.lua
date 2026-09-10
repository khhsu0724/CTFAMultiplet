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

dofile("hF_terms.lua")
Hhop = NewOperator("Number", NF, hop_i, hop_j, hop_v)

Hval = Udd + Uppx + Uppy + Hhop
Hex  = Hval + Upd - zeta*ld

CORE = "111111 0000000000 000000 000000"
VAL  = "000000 1111111111 111111 111111"

function levels(H, ncore, nval, n, tag)
  local R = {NF, NB, {CORE, ncore, ncore}, {VAL, nval, nval}}
  local P = Eigensystem(H, R, n)
  if type(P) ~= "table" then P = {P} end
  local E0 = Complex.Re(P[1]*H*P[1])
  local keys, cnt = {}, {}
  for i=1,#P do
    local k = string.format("%.5f", Complex.Re(P[i]*H*P[i]) - E0)
    if not cnt[k] then cnt[k]=0; keys[#keys+1]=k end
    cnt[k]=cnt[k]+1
  end
  print("### "..tag.."  (E0 = "..string.format("%.6f",E0)..")")
  for _,k in ipairs(keys) do print(string.format("LEVEL %s %12s %d", tag, k, cnt[k])) end
  return P, E0
end


Psis, Egs = levels(Hval, 0, 2, 40, "GSF")
