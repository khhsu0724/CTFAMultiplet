#!/usr/bin/env python3
"""NiO equilibrium RIXS maps + XAS overlay, following plot_RIXS_XAS.ipynb."""
import numpy as np, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

NEDOS = 1000
SHIFT = 862.9
XLIM  = [0, 3.5]
YLIM  = [855, 865]

def read_dir_rixs(dir_name, edge="L", pol_in="XYZ", pol_out="XYZ", nedos=NEDOS):
    """Sum the cross-polarisation RIXS files, exactly as the notebook does."""
    z = None
    for pin in pol_in:
        for pout in pol_out:
            if pin == pout: continue
            f = os.path.join(dir_name, "RIXS_%sedge_%s_%s.txt" % (edge, pin, pout))
            r = np.genfromtxt(f, usecols=np.arange(0, 3))
            x, y, zz = r[1:, 0], r[1:, 1], r[1:, 2]
            z = zz if z is None else z + zz
    n = int(x.size / nedos)
    return (np.resize(x, (n, nedos)), np.resize(y, (n, nedos)), np.resize(z, (n, nedos)))

def get_XAS_all(filedir, edge="L"):
    tot = None
    for p in "XYZ":
        r = np.genfromtxt(os.path.join(filedir, "XAS_%sedge_%s.txt" % (edge, p)),
                          usecols=np.arange(0, 2))
        x, y = r[1:, 0], r[1:, 1]
        tot = y if tot is None else tot + y
    return x, tot

try:    plt.rcParams['font.family'] = 'Helvetica'
except Exception: pass
plt.rcParams['font.size'] = 11

PANELS = sys.argv[1:]
labels = {"REF": "Your stored reference",
          "LEG": "New code, upstream exch. sign",
          "FIX": "New code, exchange-sign fix"}
dirs   = dict(zip([p.split("=")[0] for p in PANELS], [p.split("=")[1] for p in PANELS]))
keys   = [p.split("=")[0] for p in PANELS]

fig, axes = plt.subplots(1, len(keys), figsize=(4.7*len(keys), 3.4), squeeze=False)
for k, (key, ax) in enumerate(zip(keys, axes[0])):
    d = dirs[key]
    xr, yr, zr = read_dir_rixs(d)
    xx, yx = get_XAS_all(d)
    energy_loss     = yr
    incident_energy = xr + SHIFT
    cbar_max = np.max(zr) * 0.3
    pcm = ax.pcolormesh(energy_loss, incident_energy, zr, shading='gouraud',
                        cmap='coolwarm', vmin=0, vmax=cbar_max)
    # XAS overlay, scaled into the first 0.8 eV
    xas_n = (yx / np.max(yx)) * 0.8
    ax.plot(xas_n, xx + SHIFT, color='white', linewidth=2)
    ax.text(0.1, YLIM[1]-0.5, 'XAS', color='white', fontsize=13, va='top', ha='left')
    axins = inset_axes(ax, width="25%", height="3%", loc='upper right', borderpad=1.6)
    cb = fig.colorbar(pcm, cax=axins, orientation='horizontal'); cb.set_ticks([])
    axins.text(-0.05, 0.5, 'Min', transform=axins.transAxes, va='center', ha='right', fontsize=8)
    axins.text(1.05, 0.5, 'Max', transform=axins.transAxes, va='center', ha='left', fontsize=8)
    ax.set_xlabel('Energy transfer (eV)')
    if k == 0: ax.set_ylabel('Incident energy (eV)')
    ax.set_xlim(XLIM); ax.set_ylim(YLIM)
    ax.tick_params(axis='both', direction='in', length=5, width=1.2, top=True, right=True)
    ax.tick_params(axis='both', direction='in', length=3, width=1.2, which='minor', top=True, right=True)
    ax.xaxis.set_major_locator(MultipleLocator(1)); ax.xaxis.set_minor_locator(MultipleLocator(0.2))
    ax.yaxis.set_major_locator(MultipleLocator(5)); ax.yaxis.set_minor_locator(MultipleLocator(1))
    for s in ax.spines.values(): s.set_color('white'); s.set_linewidth(1.5)
    ax.set_title(labels.get(key, key), fontsize=10.5, pad=6)
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nio_rixs_maps.png")
fig.savefig(out, dpi=170, facecolor='white'); print("wrote", out)
