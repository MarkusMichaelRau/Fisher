""" Get the fisher matrix

Parameters:
-----------
lmin:  minimum l
lmax:  maximum l
args*: the parameter set in the respective order
       of the matrix

"""

import numpy as np
import sys

def get_autocorr_inds(orderings):
    autocorr_inds = []
    for i, order in enumerate(orderings):
        if order[0] == order[1]:
            autocorr_inds.append(i)
    return autocorr_inds

if (__name__ == '__main__'):

    # read in args
    args = sys.argv[1:]
    name = args[0]

    # load in the data
    w0_deriv = np.loadtxt("deriv_w0.dat")[:, 1:]
    c_ells = np.loadtxt("Cl_fid.dat")[:, 1:]
    num_dens = np.loadtxt("num_dens_lensing.dat")
    orderings = np.loadtxt("ordering_fid.dat")
    nbins = num_dens.shape[0]

    # we only want autocorrelation info
    autocorr_inds = get_autocorr_inds(orderings)

    # calc the w0 signal to noise ratio
    w0_sn_ratio =  np.sum((w0_deriv[:, autocorr_inds]**2 / \
                          (c_ells[:, autocorr_inds] + num_dens)**2))

    # save the output
    fname = "w0_sn_ratio_%s.dat"%(name)
    np.savetxt(fname, w0_sn_ratio)
