#!/usr/bin/env python

"""
Script to see how changing number of bins/delta for derivative
affects w0 signal to noise ratio
"""

import os
import sys
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

OUTDIR = "out_w0_sn_all"

########################################################
# First, generate tomo and num_dens data 
########################################################

def gen_tomo_data(probe, nbins=5):
    print("Making %s pdf"%probe)
    script = "./tomo_and_num_dens.py"
    subprocess.call([script, str(nbins), probe])
    print("Made %s pdf"%probe)


########################################################
# Make w0 sn ratio
######################################################## 

def get_subdir(nbins, probe, delta=0.15):
    return "/%s_bins_%d_delta_%.8e/"%(probe, nbins, delta) 

def make_dir(newdir):
    # make a new output folder if it doesn't exist
    if not os.path.exists(newdir):
        os.makedirs(newdir)

def gen_w0_sn_ratio(probe, delta=0.1):
    subprocess.call(["echo", "Generating w0 signal to noise ratio"])
    dndz_path = "tomo_%s.dat"%(probe)
    num_dens_path = "num_dens_%s.dat"%(probe)
    script = "./w0_sn_ratio_cluster.sh"
    args = [script, dndz_path, num_dens_path, probe, str(delta)]
    val = subprocess.call(args)
    subprocess.call(["echo", "Generate w0 signal to noise ratio"])

def save_w0_sn_data(nbins, delta, probe):
    # moves data from out_wo_sn_ratio to outdir
    # assumes fisher matrix (and other data) has been made
    outdir = OUTDIR
    make_dir(outdir)

    old_out_dir = "out_w0_sn_ratio/"
    if os.path.exists(old_out_dir):
        out_subdir = outdir + get_subdir(nbins, probe, delta)
        # make_dir(out_subdir)
        subprocess.call(["cp", "-R", old_out_dir, out_subdir])
        return out_subdir + " w0_sn_ratio.dat"
    else:
        print "Error: Can't find w0 sn ratio matrix"
        sys.exit(1)  

########################################################

def main():

    # Takes upto three arguments: for nbins, 
    # fractional difference for derivative/FOM tunings
    # probet type
    if len(sys.argv) != 2 and len(sys.argv) != 3 and len(sys.argv) != 4:
        print("Incorrect number of arguments."),
        print("Provide the number of bins")
        print("Provide fractional delta for derivative")
        print("Provide probe type: lensing or clustering")
        sys.exit(1)

    ##########################################
    # read args
    ##########################################
    nbins = int(sys.argv[1])
    
    if len(sys.argv) == 3 or len(sys.argv) == 4:
        delta = float(sys.argv[2])
    elif len(sys.argv) == 2:
        delta = 0.15
        probe = "lensing"

    if len(sys.argv) == 4:
        probe = sys.argv[3] #use fisher_cluster_with_repeat.sh to reduce calcs
    else:
        probe = "lensing"

    ##########################################
    # make data
    ##########################################

    # make galaxy-redshift pdfs
    # gen_tomo_data(probe, nbins)

    # get signal-to-noise ratio
    gen_w0_sn_ratio(probe, delta)    

    ##########################################
    # save data
    ##########################################
    
    w0_sn_path = save_w0_sn_data(nbins, delta, probe)
    
if __name__ == '__main__':
    main()
