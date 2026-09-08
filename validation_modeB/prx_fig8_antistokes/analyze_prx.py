#!/usr/bin/env python3
"""PRX Fig.7/8 anti-Stokes: our corrected Mode B vs the paper-1 Fortran reference.

Check A -- Fig.8: loss profile integrated over the pre-edge incident window (Fig.7(d) dashed box)
                  at t<t0 / t~t0 / t0<t, XX and XY.
Check B -- ANTISTOKES_FORTRAN_REFERENCE_CHECK.md finding 2: anti-Stokes intensity
                  (max over loss in [-1.0,-0.75]) vs incident energy; is there a local max at
                  ab~-7.3 dominating the main resonance at ab~-4.99?
"""
import numpy as np, os
SP  = os.path.dirname(os.path.abspath(__file__))
B   = os.path.expanduser("~/Documents/GradSchool/Devereaux Group/Herr_Jost/tr_spectroscopy/Ni-L-RIXS")
REF = {"XX": B+"/new_shape_function/XX/16_points_new_shape_function_XX",
       "XY": B+"/new_shape_function/XY"}
BOX = (-8.702, -5.4)          # pre-edge sideband, Fig.7(d) dashed box
AS  = (-1.0, -0.75)           # anti-Stokes window used by the reference-check doc
TM  = [-320,-52,0,32,140,160,200,284,362,446,536,626,710,788,878,968,1052,1136,1214]
DTF = 2*np.pi/320
MYDELAYS = np.linspace(-6.2831853071795865, 22.305305514018902, 15)

def togrid(x, y, z):
    ab, ls = np.unique(x), np.unique(y)
    G = np.zeros((len(ls), len(ab)))
    for a, b_, c in zip(x, y, z):
        G[np.searchsorted(ls, b_), np.searchsorted(ab, a)] = c
    return ab, ls, G

def ref(pol, idx):
    d = np.loadtxt(os.path.join(REF[pol], "trSPEC_%s_%04d.dat" % (pol, 100*idx+10)))
    return togrid(d[:,0], d[:,1], d[:,2])

def mine(tag, di):
    f = os.path.join(SP, tag, "trRIXS2_chan0_t%f.txt" % MYDELAYS[di])
    d = np.loadtxt(f, skiprows=1)
    return togrid(d[:,0], d[:,1], d[:,2])

def losscut(ab, ls, G, box=BOX):
    c = (ab >= box[0]) & (ab <= box[1])
    return ls, G[:, c].sum(axis=1)

def as_scan(ab, ls, G):
    m = (ls >= AS[0]) & (ls <= AS[1])
    return ab, G[m, :].max(axis=0)

def prof(ls, p, lo=-2.0, hi=0.0, step=2):
    m = (ls >= lo) & (ls <= hi)
    x, q = ls[m], p[m]/p[m].max()
    return x[::step], q[::step]

def bumpstats(x, q, c=-0.9, w=0.35):
    """local max nearest c, and its contrast against the minimum between it and the elastic rise"""
    sel = abs(x-c) < w
    if not sel.any(): return None
    k = np.arange(len(x))[sel][np.argmax(q[sel])]
    right = q[(x > x[k]) & (x < x[k]+0.55)]
    left  = q[(x < x[k]) & (x > x[k]-0.65)]
    if len(right) == 0 or len(left) == 0: return x[k], np.nan
    return x[k], q[k]/max(right.min(), left.min())

print("="*78)
print("CHECK A -- Fig.8: energy-loss profile, integrated over ab in [%.2f,%.2f]"%BOX)
print("="*78)
REFIDX = {"t<t0":0, "t~t0":2, "t0<t":17}          # tm = -320, 0, 1136
MYIDX  = {"t<t0":0, "t~t0":3, "t0<t":14}          # -6.283, 0.0(-0.157 nearest is idx3=-0.157? see below)
MYIDX  = {"t<t0":0, "t~t0":3, "t0<t":14}
for pol in ("XX","XY"):
    print("\n--- %s polarization ---"%pol)
    for tag,label in [("REFERENCE (Fortran, coherent convention)","ref"),
                      ("OURS  corrected  sum_l|A_l|^2","prx_%s"%pol),
                      ("OURS  coherent   |sum_l A_l|^2","prx_%s_coh"%pol)]:
        for tname in ("t<t0","t~t0","t0<t"):
            try:
                if label == "ref":
                    ab,ls,G = ref(pol, REFIDX[tname]); dt = TM[REFIDX[tname]]*DTF
                else:
                    ab,ls,G = mine(label, MYIDX[tname]); dt = MYDELAYS[MYIDX[tname]]
            except Exception as e:
                print("   %-40s %-6s  [missing: %s]"%(tag,tname,type(e).__name__)); continue
            ls2,p = losscut(ab,ls,G)
            x,q = prof(ls2,p)
            bs = bumpstats(x,q)
            head = "%-40s"%tag if tname=="t<t0" else " "*40
            print("%s %-6s dt=%+6.2f  "%(head,tname,dt) + " ".join("%5.3f"%v for v in q))
            if bs: print(" "*48 + "-> local max near %+0.3f eV, contrast %.2f"%bs)
    print("      w (eV):"+" "*41+" ".join("%5.2f"%v for v in x))

print()
print("="*78)
print("CHECK B -- anti-Stokes intensity (max over loss in [%.2f,%.2f]) vs incident energy"%AS)
print("="*78)
for pol in ("XX","XY"):
    print("\n--- %s ---"%pol)
    for tag,label,pre,post in [("REFERENCE","ref",0,17),
                               ("OURS corrected","prx_%s"%pol,0,14),
                               ("OURS coherent","prx_%s_coh"%pol,0,14)]:
        try:
            g = ref if label=="ref" else (lambda i: mine(label,i))
            ab0,ls0,G0 = (ref(pol,pre) if label=="ref" else mine(label,pre))
            ab1,ls1,G1 = (ref(pol,post) if label=="ref" else mine(label,post))
        except Exception as e:
            print("  %-16s [missing]"%tag); continue
        a,ppre  = as_scan(ab0,ls0,G0)
        a,ppost = as_scan(ab1,ls1,G1)
        sel = (a>=-8.75)&(a<=-3.0)
        a2,ppre2,ppost2 = a[sel],ppre[sel],ppost[sel]
        k73 = np.argmin(abs(a2+7.3)); k50 = np.argmin(abs(a2+4.99))
        print("  %-16s post-pump peak at ab=%+0.3f | I(-7.3)=%.3e I(-4.99)=%.3e ratio=%.1fx | pump enh. at -7.3 = %.1fx"%(
            tag, a2[np.argmax(ppost2)], ppost2[k73], ppost2[k50],
            ppost2[k73]/max(ppost2[k50],1e-30), ppost2[k73]/max(ppre2[k73],1e-30)))
        s = ppost2/ppost2.max()
        print("     shape vs ab: "+" ".join("%4.2f"%v for v in s[::max(1,len(s)//14)]))
        print("     ab         : "+" ".join("%4.1f"%v for v in a2[::max(1,len(s)//14)]))
