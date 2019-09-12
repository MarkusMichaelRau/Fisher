#!/usr/bin/env python 

"""
Script to calculate derivatives of the C_ells in the 
fisher forecast pipeline.

Uses an n-point stencil where n is either 3,5,7, or 9,
with error orders of 2,4,6, or 8 respectively.

Command Line Parameters:
------------------------
    para_id : int
        index (starting at 1) of parameter in the Fisher matrix
        in the follwoing order:
        "om_m", "w0", "h0", "A_s", "om_b", "n_s", "galbias", "wa"        
    step_size : float
        fractional difference to be used in the approximation
        with respect to the above parameter
    dNdz_path : string
        path to binned dNdz file
    order : int, optional
        order of the error in the approximation used,
        as described above, defaults to 2    
    cosmo_run_sh_path : str, optional
        path to your cosmo_run.sh script,
        defaults to "cosmo_run_cosmosis.sh"
        change in the __main__ block if you need to

EXAMPLE: if you want to calculate the derivative of C_ell wrt 
the parameter h0 with 5 point stencil (error order of 4),
with a fractional step size of 0.02,
and a dNdz file located at 'tomo_source.dat', call this script as
./get_stencil_deriv_cluster.py 3 0.02 tomo_source.dat 4 

WARNING: if the cosmo_run.sh script does not have a bash shebang, this
script fails
"""

import subprocess
import sys
import os
import numpy as np
# import numdifftools as nd

def get_difference_coeffs(order):
    """ see https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference

    Gets the finite difference coefficients for 
    central first derivative upto accuracy O(h^order)    

    Parameters:
    -----------
    order: int
        the order of the error; must be even and between 2 and 8
        defaults to 2 otherwise

    Returns:
    --------
    coeffs: dict 
        keys are the prefactors of the step for each function call, 
        and the values are the associated coefficients 
    """

    if order == 2: #3 point stencil
        return {-1: -1./2, 0: 0, 1: 1./2}
    elif order == 4: #5 point stencil
        return {-2: 1./12, -1: -2./3, 0: 0, 1: 2./3, 2: -1./12}
    elif order == 6: #7 point stencil
        return {-3: -1./60, -2: 3./20, -1: -3./4, 0: 0, 1: 3./4, 2: -3./20, 3: 1./60}
    elif order == 8: # 9 point stencil
        return {-4: 1./280, -3: -4./105, -2: 1./5, -1: -4./5, 0: 0, 1: 4./5, 2: -1./5, 3: 4./105, 4: -1./280}
    else:  # default to 3 point stencil
        return {-1: -1./2, 0: 0, 1: 1./2}

def load_fid_vals(path):
    """Loads values of fiducial parameters from a given file

    Parameters:
    -----------
    path: string
        path to file containing fiducial values of parameters

    Returns:
    --------
    fids: dict
        dictionary mapping paramaters to their fiducial values  
    
    """

    fids = dict()
    with open(path, 'r') as f:
        for line in f:
            if line[0] != '#' and line[0] != '\n':
                items = line[:-1].split('=', 1)
                fids[items[0]] = float(items[1])
    return fids

def get_file_name(name, step_mult):
    """
    Makes a name for outputted C_ells files or ordering log file
    based on what term in the finite difference approximation
    is being calculated. This is characterized by the coefficient
    multiplying the step size in each function call

    Parameters:
    -----------
    name: string
        name of file to change
    step_mult: int, float
        coefficient of the step size for a given function call
    """
    
    if step_mult > 0:
        # forward difference
        type_str = "upper"
    elif step_mult < 0:
        # backward difference
        type_str = "lower"
    # for a central derivative, 
    # the term corresponding to no step has 
    # a 0 coefficient, so we ignore it
    num = abs(step_mult)
    if name == ".Cl_out.dat":
        name = ".Cl"
    elif name == ".log_ordering.dat":
        name = ".log_ordering"
    else:
        name = name[:-4]
    return "%s_%s_%d_step.dat"%(name, type_str, num)

def step_para_val(fid_vals, para, step_mult, step):
    """
    in our approximation, we have terms like
    f(para1, ..., para + step_mult * step, ...)
    this function returns all the parameters after varying 
    a given parameter by 'para -> para + step_mult * step'
    
    Parameters:
    -----------
    fid_vals: dict
        maps parameters to fiducial values    
    para: string
        name of parameter to change
    step_mult: int, float
        coefficient of the step size for a given function call
    step: float
        value of step size in finite difference

    Returns:
    --------
    fid_valse: dict 
        maps parameters to fiducial values after varying a parameter 
    """

    fid_vals[para] = fid_vals[para] + step_mult * step
    return fid_vals

def calc_deriv_terms(para, order, dNdz_path, cosmo_run_sh_path, step_size=0.02):
    """
    Calculates the C_ells in all terms the finite difference approximation
    with a specified accuracy

    Paramaters
    ----------
    para : string
        name of the parameter
    order : int
        order of the error in the derivative
    dNdz_path : string
        path of dNdz file, first column contains z values
        the rest have the binned dNdz values
    cosmo_run_sh_path : string
        path to bash script that calls cosmosis
    step_size : float, optional
        fractional step size to be used in the finite difference calculation

    Returns:
    --------
    coeffs : dict    
        keys are the prefactors of the step for each function call,
        and the values are the associated coefficients
    step : float
        absolute step size to be used in finite difference calculation
    term_paths : dict
        maps the coefficient of the step to a string that contains a
        path to the outputted C_ells for all the terms in the approximation 
    """

    # initial parameters
    fid_vals_orig = load_fid_vals("fid_values.sh")
    paras = ["om_m", "w0", "h0", "A_s", "om_b", "n_s", "galbias", "wa"]
    lmin = fid_vals_orig['lmin']
    lmax = fid_vals_orig['lmax']
    # convert relative difference to absolute difference
    if para == "wa":
        # since wa has a fiducial value of zero, we use w0 as a reference value
        step = abs(fid_vals_orig["w0"]) * step_size
    else:
        step = abs(fid_vals_orig[para]) * step_size
    coeffs = get_difference_coeffs(order)
    term_paths = dict()

    # for all possible finite difference coefficients,
    # we calculate a term in the approximation
    for step_mult in coeffs:
        coeff = coeffs[step_mult]
        # if the coefficient is 0 we can ignore the term
        if np.allclose(coeff, 0):
            continue
   
        # print step to the terminal 
        msg = "Changing paramater %s from %f by %d * %f"%(para, fid_vals_orig[para], step_mult, step) 
        subprocess.call(["echo", msg])

        # load the fiducial values
        fid_vals = load_fid_vals("fid_values.sh")    
        #step the parameter we want to vary by the given step size
        fid_vals = step_para_val(fid_vals, para, step_mult, step)
        # now we call cosmosis and get the C_ells
        subprocess.call([cosmo_run_sh_path, dNdz_path, str(lmin), str(lmax)] + [str(fid_vals[p]) for p in paras])
        # organize the outputs so that we can load them later
        out_Cl_path = get_file_name(".Cl_out.dat", step_mult)
        term_paths[step_mult] = out_Cl_path
        subprocess.call(["mv", ".Cl_out.dat", out_Cl_path])
        subprocess.call(["mv", ".log_ordering.dat", get_file_name(".log_ordering.dat", step_mult)])

    return coeffs, step, term_paths

"""
def get_calced_deriv_terms(para, order, step_size, maindir):
    maindir = "out_FOM/run_w03"
    folders = [folder for folder in os.listdir(maindir) if folder[:5] == "sbins"]
    fid_vals_orig = load_fid_vals("fid_values.sh")
    
    if para == "wa":
        # since wa has a fiducial value of zero, we use w0 as a reference value 67                                                                                 
        step = abs(fid_vals_orig["w0"]) * step_size
    else:
        step = abs(fid_vals_orig[para]) * step_size
    coeffs = get_difference_coeffs(order)
    term_paths = dict()
    delta_dir = "/sbins_5_delta_%.8e/"%(step_size)
    for step_mult in coeffs:
        coeff = coeffs[step_mult]
        term_path_orig = os.path.join(maindir, delta_diri, "Cl")
        subprocess.call(["cp", ])        

    return coeffs, step, term_paths
"""


def calc_deriv(para, order, coeffs, step, term_paths):
    """
    Uses C_ells at various values of the paramater to approximate the derivative
    using an (order+1) point stencil

    Parameters:
    -----------
    para : string
        name of the parameter
    order : int
        order of the error in the derivative
    coeffs : dict    
        keys are the prefactors of the step for each function call,
        and the values are the associated coefficients
    step : float
        absolute step size to be used in finite difference calculation
    term_paths : dict
        maps the coefficient of the step to a string that contains a 
        path to the outputted C_ells for all the terms in the approximation

    Returns:
    --------
    deriv : array
        columnwise derivatives of the C_ells,
        the first column contains the ell values 
    """
    deriv = None
    i = 0
    for step_mult in coeffs:
        if step_mult == 0:
            continue # ignore case where coefficient is 0
        coeff = coeffs[step_mult]
        term = np.loadtxt(term_paths[step_mult])
        if i == 0:
            deriv = np.zeros(term.shape)
            deriv[:, 0] = term[:, 0] # set the ell values
        deriv[:, 1:] += coeff * term[:, 1:]
        i += 1
    deriv[:, 1:] /= step
    return deriv

def save_deriv(deriv, para):
    """
    Saves the derivative array into a file
    
    Parameters:
    -----------
    deriv : array
        columnwise derivatives of the C_ells,
        the first column contains the ell values
    para : string
        name of the parameter
    """
    deriv_path = "deriv_%s.dat"%(para)
    np.savetxt(deriv_path, deriv)

def convert_para_id_to_str(para_id):
    """
    Gets parameter name from id

    Parameters:
    -----------
    para_id : int
        integer corresponding to which parameter we want to take the derivative of

    Returns:
    --------
    para : string
        name of parameter we want
    """
    paras = ["om_m", "w0", "h0", "A_s", "om_b", "n_s", "galbias", "wa"]
    return paras[para_id - 1]

if (__name__ == "__main__"):

    # get system args
    args = sys.argv[1:]
    if (len(args) not in [3, 4, 5]):
        subprocess.call(["echo", "Wrong number of command line arguments"])
        sys.exit(1)

    # initialize variables
    para_id = int(args[0])
    para = convert_para_id_to_str(para_id)
    step_size = float(args[1])
    dNdz_path = args[2]

    # if fourth argument is not supplied, default to 
    # a error order of 2
    if (len(args) == 4):
        order = int(args[3])
    else:
        order = 2

    # if fifth argument is not supplied, default to 
    # some cosmo_run.sh script
    if (len(args) == 5):
        cosmo_run_sh_path = args[4]
    else:
        cosmo_run_sh_path = "./cosmo_run_cluster.sh" #change this to your cosmo_run.sh script

    coeffs, step, term_paths = calc_deriv_terms(para, order, dNdz_path, cosmo_run_sh_path, step_size)
    deriv = calc_deriv(para, order, coeffs, step, term_paths)
    save_deriv(deriv, para)

