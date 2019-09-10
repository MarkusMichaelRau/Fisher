"""Divide the C_l's by l(l + 1)/(2*pi)

Important for the calculation of cov-matrix
and derivatives

Command Line arguments: inputfilename outputfilename

"""

import numpy as np
import sys

if (__name__ == '__main__'):

    args = sys.argv[1:]

    #python interpolate_cl.py .Cl_out.dat .Cl_out_interp.dat $2 $3
    #TODO: Add different versions for the theoretical covariance matrices
    if (len(args) != 3):
        print "Wrong number of command line arguments"
        sys.exit(1)

    data = np.loadtxt(args[0])
    out_fname = args[1]

    #generate equal spaced value ranges

    #allocate new numpy array
    divident = (data[:, 0]*(data[:, 0] + 1.0))/(2.0 * np.pi)

    data[:, 1:] = np.divide(data[:, 1:].T, divident).T

    if args[2] == 'lensing':
        #request a lensing power spectrum --> class outputs lensing potential
        #power spectrum --> transform it into kappa power spectrum
        # cl_kappa = (l^4)/4 cl_phi (flat sky)(see https://arxiv.org/pdf/1704.01054.pdf, 2.14)
        multiplier = (data[:, 0]**4)/4.
        data[:, 1:] = np.multiply(data[:, 1:].T, multiplier).T
    elif args[2] == 'nolensing':
        pass
    else:
        raise ValueError('Invalid argument to divide_prefac. Either lensing or nolensing')

    np.savetxt(X=data, fname=out_fname)
