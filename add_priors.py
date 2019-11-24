import sys
import numpy as np

def main():
    fisher_path = sys.argv[1]
    fisher_out_path = sys.argv[2]

    fisher = np.loadtxt(fisher_path)

    # load priors for params
    sigma_om_m = 0           #don't set prior
    sigma_w0 = 0             #don't set prior
    sigma_h0 = 0.0057        #set prior
    sigma_A_s = 0            #don't set prior
    sigma_om_b = 0.00066     #set prior
    sigma_n_s = 0            #don't set prior
    sigma_galbias = 0        #don't set prior
    sigma_wa = 0             #don't set prior

    priors = np.diag([sigma_om_m, sigma_w0, sigma_h0, 
                      sigma_A_s, sigma_om_b, sigma_n_s, 
                      sigma_wa])

    # get corresponding fisher matrix by squaring to get a variance
    # matrix, and then inverting non-zero values
    priors[np.nonzero(priors)] = priors[np.nonzero(priors)]**(-2.)
    fisher_with_priors = fisher + priors

    np.savetxt(X=fisher_with_priors, fname=fisher_out_path)

if __name__ == "__main__":
    main()
