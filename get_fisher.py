""" Get the fisher matrix

Parameters:
-----------
lmin:  minimum l
lmax:  maximum l
args*: the parameter set in the respective order
       of the matrix

"""

import numpy as np
from numpy.linalg import inv
import sys


def get_fisher_comp(para_1, para_2, ells, derivs):
    """ Return the fisher matrix entry
    Parameters:
    -----------
    para1 / para2 : entries for the fisher matrix
    llmin / llmax: min l to maxl
    derivs: list of np arrays of derivs

    Returns:
    --------
    out_comp: fisher component
    """

    #first read in the derivativ_files
    # deriv_omega_m.dat
    deriv_para_1 = derivs[para_1] #np.loadtxt('deriv_' + para_1 + '.dat')
    deriv_para_2 = derivs[para_2] #np.loadtxt('deriv_' + para_2 + '.dat')
    lvals = ells
    if deriv_para_1.shape[1] == 2:
        print("get fisher matrix for single bin case")
        return single_bin_fisher(deriv_para_1, deriv_para_2, lvals)
    else:
        return multi_bin_fisher(deriv_para_1, deriv_para_2, lvals)

def single_bin_fisher(deriv_para_1, deriv_para_2, lvals):
    cov_matrix = np.loadtxt('output_covmat/onebin.mat')
    cov_matrix = inv(cov_matrix) #inverse covariance matrix
    return deriv_para_1[:, 1].dot(cov_matrix).dot(deriv_para_2[:, 1])

def multi_bin_fisher(deriv_para_1, deriv_para_2, lvals):
    fisher = 0.0
    for idx, el in enumerate(lvals):
        #print(idx)
        #first read in the datafiles
        #covariance matrix:
        left_vec = deriv_para_1[idx, 1:]
        right_vec = deriv_para_2[idx, 1:]

        cov_mat = np.loadtxt('output_covmat_binned/' + str(el) + '.mat')  # for cosmosis --> '.0.mat'
        #invert cov_mat
        cov_mat = inv(cov_mat)
        fisher = fisher + left_vec.dot(cov_mat.dot(right_vec))
    return fisher

if (__name__ == '__main__'):

    #python cov_mat fskyvalue .pasted_bins.dat .input_covmat.dat output_covmat/
    #get system args
    #python get_fisher.py $lmin $lmax $om_m $w0 $h0 $A_s $om_b $n_s
    #python get_fisher.py $lmin $lmax $om_m $w0 $h0 $A_s $om_b $n_s
    args = sys.argv[1:]

    lmin = float(args[0])
    lmax = float(args[1])
    ells = np.loadtxt("Cl_fid.dat")[:, 0]
    para_strings = []
    derivs = {}
    for i in range(2, len(args)):
        para_strings.append(args[i])
        derivs[args[i]] = np.loadtxt('deriv_' + args[i] + '.dat')

    #allocate the respective array for the fisher matrix
    fisher = np.zeros((len(para_strings), len(para_strings)))
    for i in range(len(para_strings)):
        for j in range(len(para_strings)):
            fisher[i, j] = get_fisher_comp(para_strings[i], para_strings[j], ells, derivs)

    #write out the fisher matrix
    np.savetxt(X=fisher, fname='fisher_out.dat')
