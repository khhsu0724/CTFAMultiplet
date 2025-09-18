import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from numpy import genfromtxt
import os
from math import pi

"""
    Author: Sean Hsu
    Date created: 11/20/2024
This file contains sample code to process RIXS files.
There are 3 types of XAS output: 
    1. exact peaks (Solver = 1,3)
    2. Continuous curve in the absorption range (Solver = 2,3)
    3. Continuous curve (Solver = 4)
Using Solver = 4 is strongly encouraged for speed and ease for plotting
"""

            ##########################
            # Option with solver = 4 #
            ##########################

def get_RIXS_iter(filedir,edge="",wipe_loss=False):
    """
    Get rixs that are solved by iterative method (BiCGS)
    """
    with open(filedir) as f:
        res_raman = np.genfromtxt(f,usecols=np.arange(0,3))
    x = res_raman[1:,0]
    y = res_raman[1:,1]
    z = res_raman[1:,2]
    if (wipe_loss):
        for i in range(x.shape[0]):
            if (x[i] < y[i]): z[i] = 0
    print(x.shape,y.shape,z.shape)
    for i in range(x.shape[0]):
        if x[i] != x[0]: break
        if (i == x.shape[0]-1): i = x.shape[0]
    print("ediv: ",i)
    y_dim = int(x.shape[0]/i)
    x = x.reshape((y_dim,i))
    y = y.reshape((y_dim,i))
    z = z.reshape((y_dim,i))
    return x,y,z

def get_RIXS_iter_all(filedir,edge="",pvin="XYZ",pvout="XYZ",cross=True):
    """
    INPUTS:
        filedir: Directory filepath
        pvin/pvout: Incoming and outgoing polarization
        cross: Cross polarization, if True then X->X, Y->Y and Z->Z are not considered
    """
    xsum = None
    ysum = None
    zsum = None
    for pin in pvin:
        for pout in pvout:
            if (cross and pin == pout): continue
            fname = filedir + "/RIXS_"+edge+"edge_"+pin+"_"+pout+".txt"
            x,y,z = get_RIXS_iter(fname,edge=edge)
            xsum = x
            ysum = y
            if (zsum is None): zsum = z
            else: zsum += z
    return xsum,ysum,zsum

### Sample Usage 
x,y,z = get_RIXS_dir(filedir,edge,cross=True)
plt.pcolormesh(x-y,x,z,cmap="terrain") # Plotting in Emission/Absorption


            ##########################
            # Option with solver < 4 #
            ##########################
            # It's a bit archaic

def lorentz(x, amp, cen, wid):
    return (amp)*(wid/((x-cen)**2 + wid**2))

def convert_eloss_to_emit(old_x_ax,old_y_ax,new_x_ax,new_y_ax,z_old):
    # emit_shift shifts emission axis since [-15,2] is not enough range for plottinig
    z_new = np.zeros(z_old.shape)
#     zmax = np.max(z_old)
    # x = eloss, y = incident
    nedos = old_y_ax.shape[0]
    if (new_x_ax.shape[0] != nedos or new_y_ax.shape[0] != nedos): 
        print("Axis not the same size")
        return z_new
    xnewmin = new_x_ax[0]
    xnewmax = new_x_ax[-1]
    ynewmin = new_y_ax[0]
    ynewmax = new_y_ax[-1]
    for xind,xval in enumerate(old_x_ax):      
        for yind,yval in enumerate(old_y_ax):
#             if (z_old[yind,xind] > zmax*0.9): print("eloss",xval, "ab: ", yval, "em", yval-xval)
            y_newind = int((yval-ynewmin)/(ynewmax-ynewmin)*nedos)
            x_newind = int((yval-xval-xnewmin)/(xnewmax-xnewmin)*nedos)
            if (y_newind < 0 or y_newind >= nedos): continue
            if (x_newind < 0 or x_newind >= nedos): continue
            z_new[y_newind,x_newind] += z_old[yind,xind]
    return z_new

def read_RIXS_file_old(flist,xdata,ydata,b_ab,b_em,in_eloss=False,exact=True,wipe_elastic=False,
                weak_elastic=[0,0],lor3D=False,legacy=False):
    # Automatically Convert to eloss, return emission or eloss depending on choice
    # Inputs: eloss => if the file is eloss or not
    nedos = xdata.size
    x,y = np.meshgrid(xdata,ydata)
    z_ab = np.zeros((xdata.size,xdata.size))
    z = np.zeros((xdata.size,xdata.size))
    if (exact):
        for fname in flist:
            with open(fname) as rixs:
                lines = rixs.readlines()[1:]
                for l in lines:
                    ab = float(l.split()[0])
                    el = float(l.split()[1])
                    peak = float(l.split()[2])
                    if (not in_eloss): el = ab-el # Automatically Convert to eloss
                    if (weak_elastic[0] != 0 and abs(el)< weak_elastic[0]): peak = peak*weak_elastic[1]
                    elind = int((el-elmin)/(elmax-elmin)*nedos)
                    abind = int((ab-abmin)/(abmax-abmin)*nedos)
                    if (peak < 1e-6): continue
    #                 print(ab,el,peak)
                    # Broaden in absorption
                    if (not lor3D): z_ab[:,elind] += lorentz(ydata,peak,ab,b_ab)
                    else: z += lorentz3D(x, y, peak, el, ab, b_ab)
    else:
        if (legacy):
            lor3D = False
            for fname in flist:
                with open(fname, 'r') as f:
                    z_ab += np.asarray([[float(str(num)) for num in line.split(' ')] for line in f])
        else: 
            lor3D = False
            for fname in flist:
                with open(fname, 'r') as f:
                    next(f)
                    for line in f:
                        eline = np.asarray([float(str(num)) for num in line.split()])
#                         if (eline[0] > 3): continue; #Wipe out d-d excitation > 3
                        loss_ind = int((eline[0]-xdata[0])/(xdata[-1]-xdata[0])*nedos)
                        z_ab[:,loss_ind] += eline[1:]
        if (weak_elastic[0] != 0):
            loss_ind = int((weak_elastic[0]-xdata[0])/(xdata[-1]-xdata[0])*nedos)
            z_ab[:,:loss_ind] *= weak_elastic[1]
    # Broaden in energy loss
    if (not lor3D):
        bd_list = []
        for yy in range(0,nedos):
            if(np.any(z_ab[:,yy])): 
                bd_list.append(yy)

        for bd_ind in bd_list: 
            for xx in range(0,nedos):
                if (z_ab[xx][bd_ind] != 0):
                    z[xx,:] += lorentz(xdata,z_ab[xx][bd_ind],xdata[bd_ind],b_em)

    # Wipe out anything with energy loss < 0
    if (wipe_elastic):
        if (not exact): zeroind = int((-2)/(15+2)*nedos)
        else: zeroind = int((-elmin)/(elmax-elmin)*nedos)
        for yy in range(0,zeroind):
            z[:,yy].fill(0)
            
    return z

def read_dir_rixs_old(dir_name,xdata,ydata,b_ab=0.3,b_em=0.3,edge = "L",pol_in = "XYZ",pol_out = "XYZ",
                 in_eloss=False,out_eloss=True,exact=True,wipe_elastic=False,weak_elastic=[0,0],lor3D=False,
                 legacy=False,abs_shift=0):
    """
    INPUTS:
        dir_name: Directory filepath
        xdata/ydata: linspace of absorption/loss axis
        b_ab/b_em: Lorentz Broadening in absorption or emission/loss
        edge: photon edge (K/L)
        pvin/pvout: Incoming and outgoing polarization
        in_eloss: input data is in loss format (THIS IS ALWAYS TRUE NOW)
        out_eloss: output data is in loss format
        wipe_elastic: Wipeout elastic line completely
        weak_elastic [a,b]: weaken peaks < (a) eV by a factor of (b)
    """
    zdata = np.zeros((xdata.size,xdata.size))
    flist = []
    for pin in pol_in:
        for pout in pol_out:  
            if (exact): filename = dir_name+"/RIXS_"+edge+"edge_"+pin+"_"+pout+".txt"
            else: filename = dir_name+"/RIXS_"+edge+"edge_"+pin+"_"+pout+"-kh.txt"
            flist.append(filename)
    zdata += read_RIXS_file_old(flist,xdata,ydata,b_ab,b_em,in_eloss,exact,wipe_elastic,weak_elastic,lor3D,legacy)
    if (abs_shift != 0):
        shift_num = int(abs_shift/(ydata[-1]-ydata[0])*ydata.size)
        zdata = np.roll(zdata, shift_num, axis=0)
    if (out_eloss): return zdata
    else: 
        if (not exact): return convert_eloss_to_emit(xdata,ydata,-xdata,ydata,zdata)
        else: return convert_eloss_to_emit(xdata,ydata,xdata,ydata,zdata)
        
def get_absorption_axis(indir,edge):
    fname = indir+"/INPUT-"+edge[0]+"E"
    nedos = 0
    abmin = 0
    abmax = 0
    with open(fname, 'r') as f:
        for line in f:
            if (line.split()[0].upper() == "AB"):
                abmin = float(line.split()[2])
                abmax = float(line.split()[3])
            if (line.split()[0].upper() == "NEDOS"):
                nedos = int(line.split()[2])
    return np.linspace(abmin,abmax,nedos)


## Sample Usage

nedos = 2000
abmin = -50
abmax = 50
elmin = -50
elmax = 50
edge = "L"
eloss = True
x_ax = np.linspace(elmin,elmax,nedos)
y_ax = np.linspace(abmin,abmax,nedos)
x,y = np.meshgrid(x_ax,y_ax)
plt.rcParams["figure.figsize"] = (6,6)
shift = 0#858.2#857.8
z = read_dir_rixs_old(dirname,x_ax,y_ax,0.5,0.5,edge[0],"XYZ","XYZ",True,eloss,exact=True,weak_elastic=[0.2,0])
plt.figure()
if (eloss): 
    c = "gnuplot"
    plt.pcolormesh(x,y+shift, z.reshape(x.shape),shading="auto",cmap=c)
    plt.ylabel("Incident Energy (eV)")
    plt.xlabel("Loss Energy (eV)")
    plt.xlim([-1,10]) 
    if (edge == "L"): plt.ylim([-15+shift,10+shift])
    else: plt.ylim([-15+shift,-5+shift])  
else: 
    c = mymap
    plt.pcolormesh(x+shift, y+shift, z.reshape(x.shape),shading="auto",cmap=c)
    plt.xlabel("Emission Energy (eV)")
    plt.ylabel("Incident Energy (eV)")
    if (edge == "L"):
        plt.ylim([-8+shift,5+shift])  
        plt.xlim([-12+shift,0+shift])
    else:
        plt.ylim([-20+shift,-5+shift])  
        plt.xlim([-25+shift,-8+shift])
    plt.plot(lims,lims,linewidth=2,linestyle='--',color="yellow")