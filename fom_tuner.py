"""
Script to see how changing the fractional delta of the derivative
calculations affects figure of merit calcs.

Fixes number of at 5 and varies delta
"""

import numpy as np
import subprocess
import os
import sys
from calc_fom import get_paras_fom

def main(para, use_h2=False):
    # set up init values
    nbins = 5
    if use_h2:
        paras = ["ommh2", "w0", "h0", "A_s", "ombh2", "n_s", "wa"]
        script = "./fisher_cluster_h2_deriv.sh"
        fom_para_1 = "ommh2"
    else:
        paras = ["om_m", "w0", "h0", "A_s", "om_b", "n_s", "wa"]
        script = "./fisher_cluster_deriv.sh"
        fom_para_1 = "om_m"
    para_inds = {}
    for i, para_i in enumerate(paras):
        para_inds[para_i] = i + 1
    para_inds["wa"] += 1
    para_ind = para_inds[para]


    # what range of deltas to consider
    delta_start = 1e-5
    delta_end = 1e-1 
    delta_num = 50
    deltas = np.logspace(np.log10(delta_start), np.log10(delta_end),
                         num=delta_num, endpoint=True)

    # init array of FOMs
    # zeroth column is deltas
    # first column is om_m A_s FoM
    # second column is w0 wa FoM
    om_m_A_s_FOMs = np.zeros(deltas.shape)
    w0_wa_FOMs = np.zeros(deltas.shape)
    foms = np.zeros((len(deltas), 3))
    foms[:, 0] = deltas

    """
    if os.path.exists("out_FOM_vs_delta/FOM_vs_delta.dat"):
        foms = np.loadtxt("out_FOM_vs_delta/FOM_vs_delta_%s.dat"%para)
    else:
        foms = np.zeros((len(deltas), 3))
        foms[:, 0] = deltas
    """

    for i in range(len(deltas))[:]:
        if i >= 0:
            delta = deltas[i]

            # print out message at each iteration
            msg = "*"*20 + "\ndelta=%.18e, run=%d/%d\n"%(delta, i+1, delta_num) + "*"*20 + "\n"
            subprocess.call(["echo", msg])

            # get FOM data
            sp_args = [script, "tomo_lensing.dat", "%d"%para_ind, "%f"%delta]
            subprocess.call(sp_args)

            # load FOM data into FOM array
            fisher = np.loadtxt("out_fisher_cluster/fisher_out.dat")
            foms[i, 1] = get_paras_fom(fisher, paras, fom_para_1, "A_s")
            foms[i, 2] = get_paras_fom(fisher, paras, "w0", "wa")

            # save data
            np.savetxt(X=foms, fname="out_FOM_vs_delta/FoM_vs_delta_%s.dat"%para)

if __name__ == "__main__":
    args = sys.argv[1:]
    para = args[0]
    if len(args) > 1:
        use_h2 = bool(int(args[1]))
    else:
        use_h2 = False
    main(para, use_h2)
