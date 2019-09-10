#!/bin/bash
#Parameters:
#-----------
#$1: bias_z_pdf
#$2: fiducial_z_pdf
#$3: number_dens
#
#first calculate the fisher matrix

./fisher_cluster.sh $2 $3

#Calculate the Cl's for the biased z_pdf
#Fiducial values


#define Fiducial values
source fid_values.sh

#run cosmosis on the fiducial values

./cosmo_run.sh $1 $lmin $lmax $om_m $w0 $h0 $A_s $om_b $n_s $galbias $wa
#run cosmosis on the fiducial values
mv .Cl_out.dat .Cl_biased.dat
mv .log_ordering.dat .log_ordering_biased.dat

#Cl's of the fiducial values are in .Cl_fid.dat
#cp out_fisher_cluster/fisher_out.dat ./
#cp out_fisher_cluster/deriv*.dat ./
python get_bias_trafo.py $lmin $lmax om_m A_s
mv parameter_bias.dat out_fisher_cluster/
mv .Cl_biased.dat out_fisher_cluster/Cl_biased.dat
rm Cl_fid.dat
mv .log_ordering_biased.dat out_fisher_cluster/ordering_biased.dat
