import sys
import numpy as np

Cl_file = sys.argv[1]

np.set_printoptions(linewidth=240)

# Manually define the number of ell's per bin - 
# this is calculated in ellbins_for_fitscov.ipynb assuming ell range is 76-999
Nell_bin = [29, 41, 57, 80, 111, 154, 215, 237]

# For simplicity, we are just going to assume we include all bins all 
# observables for now, because this matches what we've been doing till now. 

Cls = np.loadtxt(Cl_file)
try:
    spectrum = np.loadtxt('./ordering_fid.dat')
except:
    spectrum = np.loadtxt('.log_ordering.dat')

# Create bins in ell.
ell_binned = [1./Nell_bin[0] * sum(Cls[0:Nell_bin[0],0])] + \
             [1./Nell_bin[i] * \
             sum(Cls[sum(Nell_bin[0:i]):sum(Nell_bin[0:i+1]),0]) for i in range(1, len(Nell_bin))]

Cl_binned = [[1./Nell_bin[0] * sum(Cls[0:Nell_bin[0],j])] + \
             [1./Nell_bin[i] * \
             sum(Cls[sum(Nell_bin[0:i]):sum(Nell_bin[0:i+1]),j]) \
                         for i in range(1, len(Nell_bin))] \
                         for j in range(1, len(spectrum)+1)]

Cl_struct = np.vstack((ell_binned, np.asarray(Cl_binned))).T

np.savetxt(Cl_file, Cl_struct)

