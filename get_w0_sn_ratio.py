""" Get the w0 sn ratio

Parameters:
-----------
name: string
    name of probe
num_dens_path: string
    path to the number density file
auto_corr_only: bool, optional
    whether or not to use only autocorrelations or also cross correlations
maindir: string, optional
    path to where data will be loaded from and saved

"""

import numpy as np
import sys
import os

def get_autocorr_inds(orderings):
    if orderings.shape == (2,) and orderings[0] == 1 and orderings[1] == 1:
        return [0]
    autocorr_inds = []
    for i, order in enumerate(orderings):
        if order[0] == order[1]:
            autocorr_inds.append(i)
    return autocorr_inds

def get_cross_corr_num_dens(orderings, num_dens):
    autocorr_inds = get_autocorr_inds(orderings)
    if orderings.shape == (2,):
        num_dens_with_cross_corr = np.zeros((1,))
    else:
        num_dens_with_cross_corr = np.zeros(len(orderings))
    num_dens_with_cross_corr[autocorr_inds] = num_dens
    return num_dens_with_cross_corr

def calc_sn_ratio(w0_deriv, c_ells, num_dens, orderings, auto_corr_only=True):
    if auto_corr_only:
        inds = get_autocorr_inds(orderings)
    else:
        if orderings.shape == (2,) and orderings[0] == 1 and orderings[1] == 1:
            inds = [0]
        else:
            inds = range(len(orderings))
            autocorr_inds = get_autocorr_inds(orderings)
        num_dens = get_cross_corr_num_dens(orderings, num_dens)
    print c_ells[0]
    print num_dens
    sn_ratio = np.sum((w0_deriv[:, inds]**2 / \
                      (c_ells[:, inds] + num_dens)**2))
    if not auto_corr_only:
        crosscorr_inds = np.delete(range(len(orderings)), autocorr_inds)
        sn_ratio += np.sum((w0_deriv[:, crosscorr_inds]**2 / \
                      (c_ells[:, crosscorr_inds])**2))
    return sn_ratio

def save_sn_ratio(sn_ratio, maindir, tunedir, auto_corr_only=True):
    if auto_corr_only:
        fname = "w0_sn_ratio.dat"
    else:
        fname = "w0_sn_ratio_with_cross_corr.dat"
    fname = os.path.join(maindir, tunedir, fname)
    np.savetxt(fname, [w0_sn_ratio])

if (__name__ == '__main__'):

    # read in args
    args = sys.argv[1:]
    name = args[0]
    num_dens_path = args[1]
    try:
        if args[2] == "False":
            auto_corr_only = False
        else:
            auto_corr_only = True
    except:
        auto_corr_only = True
    try:
        maindir = args[3]
    except:
        maindir = "."
    try:
        tunedir = args[4]
    except:
        tunedir = "."

    print(name, num_dens_path, auto_corr_only, maindir, tunedir)

    # load in the data
    w0_deriv = np.loadtxt(os.path.join(maindir, tunedir, "deriv_w0.dat"))[:, 1:]
    c_ells = np.loadtxt(os.path.join(maindir, "Cl_fid.dat"))[:, 1:]
    num_dens = np.loadtxt(num_dens_path)
    orderings = np.loadtxt(os.path.join(maindir, "ordering_fid.dat"))
    if num_dens.shape == tuple():
        nbins = 1
    else:
        nbins = num_dens.shape[0]

    # calc and save the w0 signal to noise ratio
    w0_sn_ratio =  calc_sn_ratio(w0_deriv, c_ells, num_dens, orderings, auto_corr_only=auto_corr_only)
    save_sn_ratio(w0_sn_ratio, maindir, tunedir, auto_corr_only=auto_corr_only)
