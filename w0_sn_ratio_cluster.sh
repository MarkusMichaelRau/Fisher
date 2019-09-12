#!/bin/bash

#this script performs the final fisher analysis
#Parameters:
#$1:dndz filepath
#$2: numdens


#properties of the survey
fsky=0.12


# #define Fiducial values
source fid_values.sh

#run cosmosis on the fiducial values

./cosmo_run_cluster.sh $1 $lmin $lmax $om_m $w0 $h0 $A_s $om_b $n_s $galbias $wa

mv .Cl_out.dat Cl_fid.dat
mv .log_ordering.dat ordering_fid.dat


# #calculate the covariance matrix from the fiducial values
# mkdir output_covmat
# python cov_mat.py $fsky Cl_fid.dat ordering_fid.dat .num_dens.dat output_covmat/

echo 'w0 derivative'
./get_stencil_deriv_cluster.py 2 0.15 $1 2

# python get_fisher.py $lmin $lmax om_m A_s
python get_w0_sn_ratio.py lensing

# # #now the resulting fisher matrix is in fisher_out.dat
# mkdir out_fisher_cluster
mkdir out_wo_sn_ratio

# cp Cl_fid.dat out_fisher_cluster/
rm .Cl_lower*.dat
rm .Cl_upper*.dat
mv deriv_*.dat out_wo_sn_ratio/

mv w0_sn_ratio_lensing.dat
rm .log_ordering_lower*.dat
rm .log_ordering_upper*.dat
# mv ordering_fid.dat out_fisher_cluster/

# mv output_covmat out_fisher_cluster
