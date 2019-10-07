#!/usr/bin/env python

"""
Script to see how changing the fractional delta of the derivative
calculations affects signal to noise ratio calcs.

Fixes number of at 5 and varies delta
"""

import numpy as np
import subprocess
import os
import shutil
from w0_sn_tuner import gen_tomo_data

def main():
    # first make sure the pdf is set up
    probe = "lensing"
    bin_scheme = "equal_num"
    nbins = 1
    print("Making tomographic bins")
    gen_tomo_data(probe, bin_scheme, nbins)

    # what range of deltas to consider
    delta_start = 0.01 #5e-3
    delta_end = 0.3 #3e-1
    delta_num = 20
    deltas = np.logspace(np.log10(delta_start), np.log10(delta_end), 
                         num=delta_num, endpoint=True)

    print("Setting up output directories")
    out_dir = "out_sn_vs_delta/%s/%d_bins/"%(bin_scheme, nbins)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, mode=0755)

    out_file = os.path.join(out_dir, "sn_vs_delta.dat")
    
    if os.path.exists(out_file):
        ans = np.loadtxt(out_file)
    else:
        ans = np.zeros((len(deltas), 2))
        ans[:, 0] = deltas

    for i in range(len(deltas))[:]:
        if i >= 0:
            delta = deltas[i]
        
            # print out message at each iteration
            msg = "*"*20 + "\ndelta=%.18e, run=%d/%d\n"%(delta, i+1, delta_num) + "*"*20 + "\n"
            subprocess.call(["echo", msg])
            
            # get FOM data
            subprocess.call(["./w0_sn_tuner.py", str(nbins), str(delta), probe, bin_scheme])

            # load FOM data into FOM array
            sn_ratio = np.loadtxt("out_w0_sn_all/%s_bins_%d_delta_%.8e/w0_sn_ratio.dat"%(probe, nbins, delta))
            ans[i, 1] = sn_ratio

            # save data
            np.savetxt(X=ans, fname=out_file)

            # delete data files and make space for new data files
            subprocess.call(["rm", "-r", "out_w0_sn_ratio"])

    # move all idnividual data folders to correct place 
    out_dir = "out_w0_sn_all/%s/%d_bins/"%(bin_scheme, nbins)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, mode=0755)

    for folder in os.listdir("out_w0_sn_all/"):
        if folder.startswith("%s"%probe):
            src = os.path.join("out_w0_sn_all", folder)
            dst = os.path.join(out_dir, folder)
            shutil.move(src, dst)
    
    for file in ["Cl_fid.dat", "ordering_fid.dat", "num_dens_lensing.dat", "tomo_lensing.dat"]:
        shutil.copy(file, os.path.join(out_dir, file))

 
if __name__ == "__main__":
    main()         
 
