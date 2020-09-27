#!/bin/bash

#this script performs the final fisher analysis
#Parameters:
# $1: dndz filepath
# $2: deriv para number
# $3: deriv step size


#properties of the survey
fsky=0.48


# #define Fiducial values
source fid_values.sh

ommh2=$(python -c "print $om_m*($h0**2)")
ombh2=$(python -c "print $om_b*($h0**2)")

echo ${2} derivative
./get_stencil_deriv_cluster.py $2 $3 $1 4 1 

python get_fisher.py $lmin $lmax ommh2 w0 h0 A_s ombh2 n_s wa

#rm .Cl_lower.dat
#rm .Cl_upper.dat
#mv deriv_*.dat out_fisher_cluster/

mv fisher_out.dat out_fisher_cluster/
#rm .log_ordering_lower.dat
#rm .log_ordering_upper.dat

