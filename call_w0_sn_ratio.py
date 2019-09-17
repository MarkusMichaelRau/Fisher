"""
Script to see how changing the fractional delta of the derivative
calculations affects signal to noise ratio calcs.

Fixes number of at 5 and varies delta
"""

import numpy as np
import subprocess
import os

def main():
    probe = "lensing"
    # first make sure the pdf is set up
    nbins = 5

    # what range of deltas to consider
    delta_start = 1e-3 #5e-3
    delta_end = 3e-1 #3e-1
    delta_num = 100
    deltas = np.logspace(np.log10(delta_start), np.log10(delta_end), 
                         num=delta_num, endpoint=True)

    if os.path.exists("out_sn_vs_delta/sn_vs_delta.dat"):
        ans = np.loadtxt("out_sn_vs_delta/sn_vs_delta.dat")
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
            subprocess.call(["./w0_sn_tuner.py", str(nbins), str(delta), probe])

            # load FOM data into FOM array
            #om_m_A_s_FOM = np.loadtxt("out_FOM/sbins_%d_delta_%.8e/FOM_om_m_A_s.dat"%(nbins, delta))
            #w0_wa_FOM = np.loadtxt("out_FOM/sbins_%d_delta_%.8e/FOM_w0_wa.dat"%(nbins, delta))
            sn_ratio = np.loadtxt("out_w0_sn_all/%s_bins_%d_delta_%.8e/w0_sn_ratio.dat"%(probe, nbins, delta))
            ans[i, 1] = sn_ratio

            #om_m_A_s_FOMs[i] = om_m_A_s_FOM[1]
            #w0_wa_FOMs[i] = w0_wa_FOM[1]

            #ans[:, 1] = om_m_A_s_FOMs
            #ans[:, 2] = w0_wa_FOMs 
            
            # save data
            np.savetxt(X=ans, fname="out_sn_vs_delta/sn_vs_delta.dat")

            # delete data files and make space for new data files
            if (i < len(deltas)-1):
                subprocess.call(["rm", "-r", "out_w0_sn_ratio"])

 
if __name__ == "__main__":
    main()         
 
