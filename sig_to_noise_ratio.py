import os
import numpy as np
import sys

def load_c_ells(c_ell_path):
    # opens Cl_fid.dat file and loads C_ells
    # each row corresponds to a unique ell
    # first column has the value of ell in question
    # subsequent columns contain the C_ell^{bin_i x bin_j}
    #
    # Note the ells start at 2 and end at 999
    # but we are intersted only in ells >= 76

    c_ells = np.loadtxt(c_ell_path)
    return c_ells

def load_cov(path, ell):
    #path = "out_FOM/sbins_5_delta_1.56268273e-01/output_covmat/%.1f.mat"%(ell)
    path = path + "%.1f.mat"%(ell)
    cov = np.loadtxt(path)
    return cov

def main():
    c_ell_path = sys.argv[1]
    cov_path = sys.argv[2]
    sig_to_noise = 0
    c_ells = load_c_ells(c_ell_path)
    for all_data in c_ells:
        ell = all_data[0]
        data  = all_data[1:] 
        #if ell >= 76:
        #    print ell, 
        cov = load_cov(cov_path, ell)
        cov_inv = np.linalg.inv(cov)
        contrib = np.linalg.multi_dot([data, cov_inv, data])
        #print contrib
        sig_to_noise += contrib
    sig_to_noise =  np.sqrt(sig_to_noise)
    #np.savetxt(X=np.array([sig_to_noise]), fname="sn_ratio.dat")
    print(sig_to_noise)

if __name__ == "__main__":
   main()
 
