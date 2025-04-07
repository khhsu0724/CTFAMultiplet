import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from math import pi
from numpy import genfromtxt

"""
    Author: Sean Hsu
    Date created: 11/20/2024
This file contains sample code to process XAS files.
There are 2 types of XAS output: 
	1. exact peaks (Solver < 4)
	2. Continouos curve (Solver = 4)
When processing data, select the correct option
"""

def lorentz(x, y0, amp, cen, wid):
    return y0 + (amp)*(wid/((x-cen)**2 + wid**2))

def get_XAS_exact(fname,xdata,b=0.1):
    ydata = np.zeros(xdata.size)
    if not os.path.isfile(fname): return ydata
    with open(fname) as xas:
        lines = xas.readlines()[1:]
        for l in lines:
            p = float(l.split()[0])
            i = float(l.split()[1])
            ydata += lorentz(xdata,0,i,p,b)
    return ydata

def get_XAS_iter(filedir):
    with open(filedir) as f:
        res_raman = np.genfromtxt(f,usecols=np.arange(0,2))
    x = res_raman[1:,0]
    y = res_raman[1:,1]
    return x,y


def read_dir_xas(dir_name,xdata=0,b=0.1,edge = "L",pol = "XYZ",extension=".txt",solver=4):
    """
    read XAS output file in the specified directory. outputs absorption (x), intensity (y)
    INPUTS:
        dir_name: Directory filepath
        xdata: for exact peaks, xdata can be np.linspace(-25,25,1000) etc
        b: broadening
        edge: absorption edge
        pol: string of "XYZ" for example
    """
    if (solver < 4): 
        ydata = np.zeros(xdata.size)
    for p in pol:
        filename = dir_name+"/XAS_"+edge+"edge_"+p+extension
            assert not (xdata is 0), print("xdata cant be 0 for exact XAS file")
            ydata += get_XAS_exact(filename,xdata,b)
        else: 
            ydata = 0
            x,y = get_XAS_iter(filename)
            if (ydata == 0):
                ydata = y
            else: ydata = ydata + y
            xdata = x
    return xdata,ydata

### Sample Usage 
edge = "L"
fdir = "./"
x,y = read_dir_xas(fdir,edge)
plt.plot(x,y,label="sim.",ls="-")