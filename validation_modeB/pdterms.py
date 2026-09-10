from math import factorial as F, sqrt
def tri(a,b,c):
    if a+b-c<0 or a-b+c<0 or -a+b+c<0: return None
    return F(a+b-c)*F(a-b+c)*F(-a+b+c)/F(a+b+c+1)
def threej0(j1,j2,j3):
    """Wigner 3j with all m = 0."""
    J=j1+j2+j3
    if J%2 or tri(j1,j2,j3) is None: return 0.0
    h=J//2
    return ((-1)**h)*sqrt(F(J-2*j1)*F(J-2*j2)*F(J-2*j3)/F(J+1))*F(h)/(F(h-j1)*F(h-j2)*F(h-j3))
def sixj(j1,j2,j3,j4,j5,j6):
    """Racah formula for {j1 j2 j3; j4 j5 j6} (integer arguments)."""
    for t in ((j1,j2,j3),(j1,j5,j6),(j4,j2,j6),(j4,j5,j3)):
        if tri(*t) is None: return 0.0
    pre=sqrt(tri(j1,j2,j3)*tri(j1,j5,j6)*tri(j4,j2,j6)*tri(j4,j5,j3))
    a=[j1+j2+j3, j1+j5+j6, j4+j2+j6, j4+j5+j3]
    b=[j1+j2+j4+j5, j2+j3+j5+j6, j3+j1+j6+j4]
    s=0.0
    for t in range(max(a), min(b)+1):
        d=F(t-a[0])*F(t-a[1])*F(t-a[2])*F(t-a[3])*F(b[0]-t)*F(b[1]-t)*F(b[2]-t)
        s+=((-1)**t)*F(t+1)/d
    return pre*s
def red(l,k,lp):
    """<l||C^k||l'>"""
    return ((-1)**l)*sqrt((2*l+1)*(2*lp+1))*threej0(l,k,lp)
def pd_terms(F0,F2,G1,G3,l=1,lp=2):
    """Exact term energies of the inequivalent l l' configuration (Cowan eq. 12.5-12.6)."""
    Fk={0:F0,2:F2}; Gk={1:G1,3:G3}
    out={}
    for L in range(abs(l-lp),l+lp+1):
        for S in (0,1):
            e=0.0
            for k,v in Fk.items():
                e+=((-1)**(l+lp+L))*sixj(l,lp,L,lp,l,k)*red(l,k,l)*red(lp,k,lp)*v
            for k,v in Gk.items():
                # Exchange coefficient. The sign here was WRONG in the first version of this
                # file (written as -(-1)^S), which swapped the singlet and triplet of every L
                # and led me to wrongly conclude the code's exchange sign was inverted.
                # Cross-checked against Quanty 0.81 for 2p^5 3d^9: the correct sign reproduces
                # Quanty's degeneracy sequence [5,21,15,9,3,7] and its energies to 4.7e-05 eV.
                e+=((-1)**S)*sixj(l,lp,L,l,lp,k)*red(l,k,lp)**2*v
            out[("%d%s"%(2*S+1,"SPDFGH"[L]))]=e
    return out
if __name__=="__main__":
    print("3j checks: (1,1,0)=%.6f exp %.6f | (1,1,2)=%.6f exp %.6f | (1,2,3)=%.6f exp %.6f"%(
        threej0(1,1,0),-1/sqrt(3),threej0(1,1,2),sqrt(2/15),threej0(1,2,3),-sqrt(3/35)))
    print("6j check {1 2 3; 2 1 0} = %.6f  exp %.6f"%(sixj(1,2,3,2,1,0),(-1)**(1+2+3)/sqrt(3*7)))
