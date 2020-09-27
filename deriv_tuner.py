import subprocess
import numpy as np
from Fisher_Forecaster import *
import subprocess
from subprocess import Popen, PIPE

def gen_derivs_in_delta_range(L=100, delta_start=0.0000000003, 
                              delta_end=0.3, delta_num=150, F=F):
    """
    Generates a list of various derivatives (om_m and w0 
    default) given a particular value for l and a range of 
    deltas.

    The deltas are fractional differences of a given 
    parameter used in the finite difference calculation
    of the derivatives.

    Creates a file where first column is the deltas,
    the other columns are the derivatives in the same
    order as the derivative files they are pulled from.

    Parameters
    ----------
    L : int
        The value of L used to pick out which C_L to use
        for the derivatives
    delta_start : float
        The starting value of the range of deltas used
        in the derivative calculation
    delta_end : float
        The end value of the range of deltas used in the
        derivative calculation
    delta_num : int
        The number of deltas to be considered between the 
        start and end (inclusive)
    """
        

    # choose what value of L to pick out
    L_START = 76
    L_INDEX = L - L_START

    # what range of deltas to consider
    deltas = np.logspace(np.log10(delta_start), np.log10(delta_end), 
                         num=delta_num, endpoint=True)

    #deriv names
    om_m    = "deriv_om_m"
    w0      = "deriv_w0"
    h0      = "deriv_h0"
    A_s     = "deriv_A_s"
    om_b    = "deriv_om_b"
    n_s     = "deriv_n_s"
    galbias = "deriv_galbias"
    wa      = "deriv_wa"

    #put derivs you want into paras
    paras = [om_m, w0, h0, A_s, om_b, n_s, wa]
    tune_para = "om_m"
    for para in paras:
        


    # init output arrays
    # assumes 5 bins and sets 11 cols
    # one for delta values
    # rest for various derivs
    COLS = 16
    outputs = [np.zeros((delta_num, COLS)) for i in range(len(paras))]
    fom_w0_wa = np.zeros(deltas.shape)
    fom_om_m_A_s = np.zeros(deltas.shape)

    # put deltas into first column
    for output in outputs:
        output[:, 0] = deltas

    for i in range(len(deltas)):
        delta = deltas[i]

        # print out message at each iteration
        msg = "*"*20 + "\ndelta=%.18e\n"%delta + "*"*20 + "\n"
        subprocess.call(["echo", msg])

        #create deriv files
        #subprocess.call(["python", "call_fisher_cluster.py", str(delta)])
        subprocess.call(["/home/nbhandar/Fisher/deriv_tuner_wrapper.sh", "%d"%para_ind, "%f"%delta])

        for j in range(len(paras)):
            # load derivative file
            path = "out_fisher_cluster/"+paras[j]+".dat"
            derivs = np.loadtxt(path)

            # pick out derivative we want
            # columns correspond to columns in deriv files 
            for col in range(1, COLS):
                outputs[j][i, col] = derivs[L_INDEX, col]
        
            # save output on the go to deal with time outs
            np.savetxt(paras[j] + "_vs_delta_L" + str(L) + \
                       "_num" + str(delta_num) + ".dat", outputs[j])

        # delete data files and make space for new data files
        if (i < len(deltas)-1):
            subprocess.call(["rm", "-r", "out_fisher_cluster"])

if __name__ == '__main__':
    probe = "lensing"
    bin_type = "equal_size"
    nbins = 5
    deriv_order = 2
    derivs_to_calc = "all"
    use_binned = True
    F = Fisher_Forecaster(probe, bin_type, nbins, deriv_order, derivs_to_calc, use_binned=use_binned)
    para_pairs_list = [("om_m", "A_s"), ("w0", "wa")]
    gen_derivs_in_delta_range(delta_start=0.003, delta_end=0.3,
                              delta_num=30, F=F)

