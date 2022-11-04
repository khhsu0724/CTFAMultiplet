#!/usr/bin/env python
import os
import numpy as np
import time
# Quick job submission, using python3
# Variables: number of holes, tenDQ, del_eff, t, K/L edge
def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.makedirs(dir)

class Input_Param:
    def __init__(self,nh, tenDQ, del_eff, delta, t):
        self.nh = nh
        self.tenDQ = tenDQ
        self.del_eff = del_eff
        self.delta = delta
        self.t = t
        
base_dir = os.getcwd()
data_dir = os.path.join(base_dir, "DATA_DIR")
input_arr = []
nh5 = True
nh4 = True
if (nh5):
    for d in np.arange(0,-2.6,-2.5):
        # 5h, 10DQ = 1.1, delta_eff = -2/-4.5
        for t in np.arange(1.0, 1.6, 0.2):
            input_arr.append(Input_Param(5,1.1,-2+d,18.49+d,t))
        # 5h, 10DQ = 1.3, delta_eff = -2/-4.5
        for t in np.arange(0.8, 1.4, 0.2):
            input_arr.append(Input_Param(5,1.3,-2+d,18.69+d,t))
        # 5h, 10DQ = 1.5, delta_eff = -2/-4.5
        for t in np.arange(0.6, 1.2, 0.2):
            input_arr.append(Input_Param(5,1.5,-2+d,18.89+d,t))
        # 5h, 10DQ = 1.7, delta_eff = -2/-4.5
        for t in np.arange(0.5, 1.1, 0.2):
            input_arr.append(Input_Param(5,1.7,-2+d,19.09+d,t))
        # 5h, 10DQ = 1.9, delta_eff = -2/-4.5
        for t in np.arange(0.4, 1, 0.2):
            input_arr.append(Input_Param(5,1.9,-2+d,19.29+d,t))
if (nh4):
    for d in np.arange(0,3.1,3):
        # 4h, 10DQ = 1.1, delta_eff = 1/4
        for t in np.arange(1.0, 1.6, 0.2):
            input_arr.append(Input_Param(4,1.1,1+d,17.1+d,t))
        # 4h, 10DQ = 1.3, delta_eff = 1/4
        for t in np.arange(0.8, 1.4, 0.2):
            input_arr.append(Input_Param(4,1.3,1+d,17.283+d,t))
        # 4h, 10DQ = 1.5, delta_eff = 1/4
        for t in np.arange(0.6, 1.2, 0.2):
            input_arr.append(Input_Param(4,1.5,1+d,17.41+d,t))
        # 4h, 10DQ = 1.7, delta_eff = 1/4
        for t in np.arange(0.5, 1.1, 0.2):
            input_arr.append(Input_Param(4,1.7,1+d,17.654+d,t))
        # 4h, 10DQ = 1.9, delta_eff = 1/4
        for t in np.arange(0.4, 1, 0.2):
            input_arr.append(Input_Param(4,1.9,1+d,17.85+d,t))

kedge = True
ledge = False
for inp in input_arr:
    input_dir = os.path.join(data_dir,"nh="+str(inp.nh),"tenDQ="+f"{inp.tenDQ:.1f}",
    					"del_eff="+f'{inp.del_eff:.1f}',"t="+f'{inp.t:.1f}')
    mkdir_p(input_dir)
    
    for edge in ["K","L"]:
        if (not kedge): 
            if(edge == "K"): continue
        if (not ledge): 
            if(edge == "L"): continue
        job_file = os.path.join(input_dir,"run%sE.slurm" % edge)
        jobname = "nh"+str(inp.nh)+"-tenDQ"+f'{inp.tenDQ:.1f}'+ \
                    "-del_eff"+f'{inp.del_eff:.1f}'+"-t"+f'{inp.t:.1f}'+"-"+edge+"E"
        input_file = os.path.join(input_dir,"INPUT-"+edge+"E")
        tenDQ = f'{inp.tenDQ:.1f}'
            
        with open(input_file,'w') as fh:
            fh.writelines("&CONTROL\n")
            fh.writelines("\tSO = 10.5 # Spin orbit coupling (eV)\n")
            fh.writelines("\tSC1 = 1.0 0 0.1\n")
            fh.writelines("\tSC2 = 6.0215 0 0.206 0 0.01436\n")
            fh.writelines("\tSC2EX = 6.0327 0 0.2191 0 0.0137\n")
            if (edge == "K"): fh.writelines("\tFG = 1.0 0.25 0 0\n")
            else: fh.writelines("\tFG = 7.4332 4.7576 6.3192 2.7072\n")
            fh.writelines("\tCF = 0 0 "+tenDQ+" "+tenDQ+" "+tenDQ+"\n")
            fh.writelines("\ttpd = "+f'{inp.t:.1f}'+"\n")
            fh.writelines("\ttpp = 0.25\n")
            fh.writelines("\tMLCT = "+str(inp.delta)+"\n")
            fh.writelines("\tOVERWRITE = True\n")
            fh.writelines("/\n")
            fh.writelines("&CELL\n")
            fh.writelines("\tCoordination = \"Square Planar\"\n")
            fh.writelines("\tSites = 1\n")
            fh.writelines("\tHoles = "+str(inp.nh)+"\n")
            fh.writelines("/\n")
            fh.writelines("&PHOTON\n")
            fh.writelines("\tXAS = true\n")
            fh.writelines("\tRIXS = true\n")
            fh.writelines("\tpvin = 1 1 1\n")
            fh.writelines("\tpvout = 1 1 1\n")
            fh.writelines("\tEdge = "+edge+"\n")
            fh.writelines("\tEloss = false\n")
            fh.writelines("\tNEDOS = 1000\n")
            fh.writelines("/\n")
        with open(job_file,'w') as fh:
            # if L edge, request much more time
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%s.job\n" % jobname)
            fh.writelines("#SBATCH --mail-user=khhsu@stanford.edu\n")
            fh.writelines("#SBATCH --mail-type=ALL\n")
            fh.writelines("#SBATCH --cpus-per-task=16\n")
            fh.writelines("#SBATCH --ntasks=1\n")
            fh.writelines("#SBATCH --nodes=1\n")
            if (inp.nh == 5): fh.writelines("#SBATCH -t 48:00:00\n")
            else:
                if (edge == "K"): fh.writelines("#SBATCH -t 24:00:00\n")
                else: fh.writelines("#SBATCH -t 36:00:00\n")
            fh.writelines("#SBATCH --mem=32GB\n")
            fh.writelines("#SBATCH -p owners,simes\n")
            fh.writelines("#SBATCH --qos=normal\n")
            fh.writelines("\n")
            fh.writelines("module load imkl icc boost\n")
            fh.writelines("module load gcc\n")
            fh.writelines("export OMP_NUM_THREADS=16\n")
            fh.writelines("srun -c 16 " + os.getcwd()+"/main "+input_dir+"/INPUT-"+edge+"E 2>&1 | tee -a "+input_dir+"/%s.out\n" %jobname)

        os.chdir(input_dir)
        short_job_file = "./run%sE.slurm" % edge
        jobid = os.popen("sbatch --parsable %s" %short_job_file).read()
        os.system('echo \"submitted job: %s\"' %jobname)
        os.system('echo \"job id: %s\"'% jobid)
        os.chdir(base_dir)
        time.sleep(1)