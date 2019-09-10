"""This module calculates the tranformation matrix for
the Fisher analysis.

Command line parameters:
------------------------
python get_bias_trafo.py bias_z_pdf fiducial_z_pdf

"""

import numpy as np
from numpy.linalg import inv
import sys


def get_bias_z_pdf_comp(para, diff_cl, llmin, llmax):
    """ Return the fisher matrix entry
    Parameters:
    -----------
    para: parameter bias to  be calculated
    cl_biased: biased cl's
    cl_unbiased: unbiased cl's
    llmin / llmax: min l to maxl

    Returns:
    --------
    out_comp: fisher component
    """

    #first read in the derivativ_files
    # deriv_omega_m.dat
    deriv_para = np.loadtxt('out_fisher_cluster/deriv_' + para + '.dat')
    lvals = np.arange(llmin, llmax + 1)
    if deriv_para.shape[1] == 2:
        print('Single bin case')
        bias = get_bias_single_bin(diff_cl, deriv_para)
    else:
        bias = 0.0
        for idx, el in enumerate(lvals):
            #print idx
            #first read in the datafiles
            #covariance matrix:
            left_vec = diff_cl[idx, 1:]
            right_vec = deriv_para[idx, 1:]
            cov_mat = np.loadtxt('out_fisher_cluster/output_covmat/' + str(el) + '.0.mat')
            #invert cov_mat
            cov_mat = inv(cov_mat)
            bias = bias + left_vec.dot(cov_mat.dot(right_vec))
    return bias

def get_bias_single_bin(diff_cl, deriv_para):
    cov_mat = np.loadtxt('out_fisher_cluster/output_covmat/onebin.mat')
    cov_mat = inv(cov_mat)
    #print np.sum(diff_cl[:, 1])
    return diff_cl[:, 1].dot(cov_mat).dot(deriv_para[:, 1])

if (__name__ == '__main__'):

    #get system args
    args = sys.argv[1:]

    lmin = int(float(args[0]))
    lmax = int(float(args[1]))
    para_strings = []
    header_string = ''
    for i in range(2, len(args)):
        para_strings.append(args[i])
        header_string = header_string + ' ' + args[i]

    #get the cl files

    cl_biased = np.loadtxt('.Cl_biased.dat')
    cl_unbiased = np.loadtxt('Cl_fid.dat')

    #print cl_biased
    #print cl_unbiased
    diff_cl = np.column_stack((cl_biased[:, 0], cl_biased[:, 1:] - cl_unbiased[:, 1:]))
    #in output_covmat are the covariance matrices

    #allocate the respective array for the fisher matrix
    bias_vec = []

    for i in range(len(para_strings)):
        bias_vec.append(get_bias_z_pdf_comp(para_strings[i], diff_cl, lmin, lmax))

    bias_vec = np.array(bias_vec)
    #now calculate the transformed fisher matrix:

    fisher_fiducial = np.loadtxt('out_fisher_cluster/fisher_out.dat')
    inv_fisher_fid = inv(fisher_fiducial)

    para_bias = inv_fisher_fid.dot(bias_vec)
    print para_bias
    np.savetxt(X=para_bias, fname='parameter_bias.dat', header=header_string)
