#!/bin/bash

#This script prints the value file (cosmology)
#for cosmosis. We assume a 7 parameter family
#
#Parameters:
#$1: index
#$2: delta difference
#$3: dndz file
#index encoding:
# 1: om_m
# 2: w0
# 3: h0
# 4: A_s
# 5: om_b
# 6: n_s
# 7: galbias
# 8: wa

#Fiducial values
source fid_values.sh

#which parameter gets modified?
if [ $1 == 1 ]; then
	om_m="$(python -c "print $om_m + abs($om_m)*$2")"
	echo $om_m
fi
if [ $1 == 2 ]; then
 	w0="$(python -c "print $w0 + abs($w0)*$2")"
fi
if [ $1 == 3 ]; then
 	h0="$(python -c "print $h0 + abs($h0)*$2")"
fi
if [ $1 == 4 ]; then
 	A_s="$(python -c "print $A_s + abs($A_s)*$2")"
fi
if [ $1 == 5 ]; then
 	om_b="$(python -c "print $om_b + abs($om_b)*$2")"
fi
if [ $1 == 6 ]; then
 	n_s="$(python -c "print $n_s + abs($n_s)*$2")"
fi
if [ $1 == 7 ]; then
	galbias="$(python -c "print $galbias + abs($galbias)*$2")"
fi
if [ $1 == 8 ]; then
	wa="$(python -c "print abs($w0)*$2")"
fi

echo 'upper'
./cosmo_run.sh $3 $lmin $lmax $om_m $w0 $h0 $A_s $om_b $n_s $galbias $wa


mv .Cl_out.dat .Cl_upper.dat
#cp .Cl_upper.dat store_upper.dat
mv .log_ordering.dat .log_ordering_upper.dat

source fid_values.sh

#which parameter gets modified?
if [ $1 == 1 ]; then
	om_m="$(python -c "print $om_m - abs($om_m)*$2")"
fi
if [ $1 == 2 ]; then
 	w0="$(python -c "print $w0 - abs($w0)*$2")"
fi
if [ $1 == 3 ]; then
 	h0="$(python -c "print $h0 - abs($h0)*$2")"
fi
if [ $1 == 4 ]; then
 	A_s="$(python -c "print $A_s - abs($A_s)*$2")"
fi
if [ $1 == 5 ]; then
 	om_b="$(python -c "print $om_b - abs($om_b)*$2")"
fi
if [ $1 == 6 ]; then
 	n_s="$(python -c "print $n_s - abs($n_s)*$2")"
fi
if [ $1 == 7 ]; then
	galbias="$(python -c "print $galbias - abs($galbias)*$2")"
fi
if [ $1 == 8 ]; then
	wa="$(python -c "print -abs($w0)*$2")"
fi


echo 'lower'
./cosmo_run.sh $3 $lmin $lmax $om_m $w0 $h0 $A_s $om_b $n_s $galbias $wa

mv .Cl_out.dat .Cl_lower.dat
#cp .Cl_lower.dat store_lower.dat
mv .log_ordering.dat .log_ordering_lower.dat

source fid_values.sh

#which parameter gets modified?
if [ $1 == 1 ]; then
	python calc_deriv.py .Cl_upper.dat .Cl_lower.dat om_m "$(python -c "print abs($om_m)*$2")"
fi
if [ $1 == 2 ]; then
	python calc_deriv.py .Cl_upper.dat .Cl_lower.dat w0 "$(python -c "print abs($w0)*$2")"
fi
if [ $1 == 3 ]; then
	python calc_deriv.py .Cl_upper.dat .Cl_lower.dat h0 "$(python -c "print abs($h0)*$2")"
fi
if [ $1 == 4 ]; then
	python calc_deriv.py .Cl_upper.dat .Cl_lower.dat A_s "$(python -c "print abs($A_s)*$2")"
fi
if [ $1 == 5 ]; then
	python calc_deriv.py .Cl_upper.dat .Cl_lower.dat om_b "$(python -c "print abs($om_b)*$2")"
fi
if [ $1 == 6 ]; then
	python calc_deriv.py .Cl_upper.dat .Cl_lower.dat n_s "$(python -c "print abs($n_s)*$2")"
fi
if [ $1 == 7 ]; then
	python calc_deriv.py .Cl_upper.dat .Cl_lower.dat galbias "$(python -c "print abs($galbias)*$2")"
fi
if [ $1 == 8 ]; then
	python calc_deriv.py .Cl_upper.dat .Cl_lower.dat wa "$(python -c "print abs($w0)*$2")"
fi

# echo 'delta w_a'
# echo "$(python -c "print abs($w0)*$2")"
