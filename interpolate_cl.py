"""Interpolate the '.pasted_bins.dat' files
such that they contain values at each integer l

Important for the calculation of cov-matrix
and derivatives

Command Line arguments: inputfilename outputfilename

"""

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import sys

def interp_cl(data, input_min_l, input_max_l):
    #generate equal spaced value ranges
    l_vec = np.arange(input_min_l, input_max_l + 1, step=1)
    
    #allocate new numpy array
    new_vec = np.zeros((len(l_vec), data.shape[1]))
    new_vec[:, 0] = l_vec

    for i in range(1, data.shape[1]):
        #first generate a spline interpolation object
        interpolation = interp(data[:, 0], data[:, i])
        new_vec[:, i] = interpolation(l_vec)
    
    return new_vec

if (__name__ == '__main__'):

    args = sys.argv[1:]

    #python interpolate_cl.py .Cl_out.dat .Cl_out_interp.dat $2 $3
    #TODO: Add different versions for the theoretical covariance matrices
    if (len(args) != 4):
        print("Wrong number of command line arguments")
        sys.exit(1)

    data = np.loadtxt(args[0])
    out_fname = args[1]
    input_min_l = int(float(args[2]))
    input_max_l = int(float(args[3]))

    new_vec = interp_cl(data, input_min_l, input_max_l)

    np.savetxt(X=new_vec, fname=out_fname)
