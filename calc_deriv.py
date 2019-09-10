
""" Author: MMRAU
Script to calculate the derivative
in the fisher forecast pipeline

Call:
python calc_deriv.py in_up in_low parameter_tag diff

in_up ascii file format l Cl: calculated for the upper limit
in_low ascii file format l Cl: calculated for the lower limit
parameter_tag: tag of the parameter that was varied
diff: value wrt. that the parameter was varied

output: save a ascii file with the derivatives
        Format: l dCldpara
"""

import numpy as np
import sys


def get_deriv(b_high, b_low, diff):
    """ Calculate the derivative
    see doc-string of the module
    """

    prefac = 1.0/(2.0 * diff)
    deriv = prefac*(b_high[:, 1:] - b_low[:, 1:])
    output = np.column_stack((b_high[:, 0], deriv))

    return output


if (__name__ == '__main__'):

    #get system args
    args = sys.argv[1:]

    #TODO: Add different versions for the theoretical covariance matrices
    if (len(args) != 4):
        print "Wrong number of command line arguments"
        sys.exit(1)

    bins_high = np.loadtxt(args[0])
    bins_low = np.loadtxt(args[1])
    tag = str(args[2])
    diff = float(args[3])
    print diff
    print bins_high[0, 1]
    print bins_low[0, 1]
    deriv = get_deriv(bins_high, bins_low, diff)
    print deriv[0, 1]
    np.savetxt(X=deriv, fname="deriv_"+tag+".dat")
