import sys
import numpy as np

# This script takes the covariance matrix in the form which is output by Markus and my Fisher code 
# and massages it into the form we need to run an emcee chain with cosmosis.

np.set_printoptions(linewidth=240)

# Set values which we define
lmin = 76
lmax = 999

# Manually define the number of ell's per bin - this is calculated in ellbins_for_fitscov.ipynb assuming ell range is 76-999
Nell_bin = [29, 41, 57, 80, 111, 154, 215, 237]
binned_ells = np.loadtxt("Cl_fid.dat")[:, 0]

# Number of redshift bins of galaxy positions and galaxy shapes
nbins = 5
ellcov_loc ='output_covmat/'

# Get simple derived quantities 
nell = lmax - lmin + 1

# Load the covariance matrices in the native form to mine and Markus' code.
print "Loading ell matrices"
ell_mats_unbinned = [np.loadtxt(ellcov_loc + str(i)+'.0.mat') for i in range(lmin,lmax+1)]

ell_mats = [1./Nell_bin[0]**2 * \
            sum(ell_mats_unbinned[0:Nell_bin[0]])] + \
           [1./Nell_bin[i]**2 * \
            sum(ell_mats_unbinned[sum(Nell_bin[0:i]):sum(Nell_bin[0:i+1])]) \
            for i in range(1, len(Nell_bin))]
print("Finished binning the covariance matrices.")

# save the binned covariances 
for i, ell in enumerate(binned_ells):
    print("output_covmat_binned/" + str(ell)+'.mat')
    np.savetxt(fname="output_covmat_binned/" + str(ell)+'.mat', X=ell_mats[i])
