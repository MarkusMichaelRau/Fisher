#!/bin/bash

#this script performs the final fisher analysis
#Parameters:
#$1:dndz filepath


#properties of the survey
fsky=0.48

# #output the number density values for
# #each BIN!!!!

#assume 10 arcmin**(-2) per bin
#cat <<EOT >> .num_dens.dat
#118200000.0
#118200000.0
#EOT
#
# cat <<EOT >> .num_dens.dat
# $2
# EOT


# #define Fiducial values
source fid_values.sh

#run cosmosis on the fiducial values
echo $lmin $lmax 

./cosmo_run_cluster.sh $1 $lmin $lmax $om_m $w0 $h0 $A_s $om_b $n_s $galbias $wa

mv .Cl_out.dat Cl_fid.dat
mv .log_ordering.dat ordering_fid.dat


#calculate the covariance matrix from the fiducial values
mkdir output_covmat
python cov_mat.py $fsky .Cl_fid.dat ordering_fid.dat num_dens_lensing.dat output_covmat/
mv .Cl_fid.dat Cl_fid_unbinned.dat
mkdir output_covmat_binned
python construct_Cov_fits_object_listcomp.py # bin the covmats

#rm .num_dens.dat

echo 'Omega matter derivative'
./get_stencil_deriv_cluster.py 1 0.15 $1 2

echo 'w0 derivative'
./get_stencil_deriv_cluster.py 2 0.15 $1 2

echo 'H0 derivative'
./get_stencil_deriv_cluster.py 3 0.15 $1 2

echo 'A_s derivative'
./get_stencil_deriv_cluster.py 4 0.15 $1 2

echo 'Omega_b derivative'
./get_stencil_deriv_cluster.py 5 0.15 $1 2

echo 'n_s derivative'
./get_stencil_deriv_cluster.py 6 0.15 $1 2

# echo 'galbias derivative'
# ./get_deriv_cluster.sh 7 0.15 $1

echo 'wa derivative'
./get_stencil_deriv_cluster.py 8 0.1 $1 2

python get_fisher.py $lmin $lmax om_m w0 h0 A_s om_b n_s wa

# #now the resulting fisher matrix is in fisher_out.dat
mkdir out_fisher_cluster

cp Cl_fid.dat out_fisher_cluster/
rm .Cl_lower.dat
rm .Cl_upper.dat
mv deriv_*.dat out_fisher_cluster/

mv fisher_out.dat out_fisher_cluster/
rm .log_ordering_lower.dat
rm .log_ordering_upper.dat
mv ordering_fid.dat out_fisher_cluster/

mv output_covmat out_fisher_cluster
mv output_covmat_binned out_fisher_cluster
