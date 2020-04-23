#!/usr/bin/env bash

#This script is meant to be run in sequencially
#CAREFUL: If ran in parallel things will get messed
#up
#It runs cosmosis and writes the C_l's in the file
#.Cl_out.dat
#The ordering of the bins is logged in
#.log_ordering.dat

#Parameters:
#
#$1: dndz_filename
#$2: lmin
#$3: lmax
#$4: om_m
#$5: w
#$6: h
#$7: A_s
#$8: om_b
#$9: n_s
#${10}: galbias
#${11}: wa

# put back factor of 1-e9 into A_s
AS=$(python -c "print $7*1e-9")

cat <<EOT >> .values.dat
[cosmological_parameters]
omega_c = $4
h0 = $6
omega_b = $8
tau = 0.08
n_s = $9
A_s = ${7}
w = $5
wa = ${11}
[galaxy_bias]
b=${10}
EOT

cat <<EOT >> .run_cosmosis.ini
[runtime]
sampler = test
root = /home/nbhandar/cosmosis

[test]
save_dir=.parameters
fatal_errors=T

[pipeline]
; the main pipeline. It's a sequence of modules to run.
modules = consistency camb halofit extrapolate_power load_nz pk_to_cl
;modules = consistency camb halofit load_nz pk_to_cl
;modules = consistency camb linear_pk extrapolate_power load_nz pk_to_cl

; the steps are:
; 1) consistency: calculate the simply derived cosmological parameters (e.g. omega_c = omega_m-omega_b)
; 2) camb: run the Boltzmann code to get the matter power spectrum
; 3) halofit: get the nonlinear matter power spectrum
; 4) extrapolate_power: extend the power spectra to high k
; 5) load_nz: get the photometric n(z) for LSST 
; 6) pk_to_cl: convert the 3D spectra into 2D tomographic C_ell with the Limber approximation

; initial parameter values and their ranges and priors
values = .values.dat

; If you want to combine with additional likelihoods such as Planck;
; then you will need to add them here, e.g.  likelihoods = xipm planck euclid lsst
likelihoods =

; extra (derived) parameter to save
extra_output =

;Control of extra info printed out
quiet=T
timing=T
debug=F


;***********************************
;Theory
;***********************************


[consistency]
file = cosmosis-standard-library/utility/consistency/consistency_interface.py

[camb]
file = cosmosis-standard-library/boltzmann/camb/camb.so
mode=all
lmax=2500
feedback=0
kmax=300.0
high_accuracy_default=T

[halofit]
file = cosmosis-standard-library/boltzmann/halofit_takahashi/halofit_interface.so
nk=500

[linear_pk]
file= /pylon5/as5fp8p/nbhandar/lensing/cosmosis_fisherforecast/halofit_linear.py

[extrapolate_power]
file=cosmosis-standard-library/boltzmann/extrapolate/extrapolate_power.py
kmax=500.0


;***********************************
; Choice of photometric redshift
;***********************************


[load_nz]
file = cosmosis-standard-library/number_density/load_nz/load_nz.py
output_section = "nz_source"
filepath = $1


;***********************************
; Calculate shear C_ells
;***********************************


[pk_to_cl]
file = cosmosis-standard-library/structure/projection/project_2d.py
ell_min = $(python -c "print float($2)")
ell_max = $(python -c "print float($3)")
n_ell=700
shear-shear = source-source
verbose = F
get_kernel_peaks=F
EOT

#cleanup the parameter and
#input files for cosmosis
#even if the run fails
#source  ~/cosmosis/config/setup-cosmosis

{ # this is my bash try block
    cosmosis .run_cosmosis.ini &&
    rm .values.dat
    rm .run_cosmosis.ini
} || { # this is catch block
    rm .values.dat
    rm .run_cosmosis.ini
}

#This part takes all the calculated cl's from
#cosmosis and pastes them into a convinient
#file. It also outputs a log of the bin ordering
#imposed

#### lensing 

for file in $(ls .parameters/shear_cl/bin_*.txt)
do
    sed '/#/d' $file > $file.new_file.txt
    awk 'NR==1{print}' $file | sed "s/#//" >> .log_ordering.dat
done
sed "s/bin_//;s/_/ /" .log_ordering.dat | awk '{print $2, $1}' > .log_ordering_new.dat
rm .log_ordering.dat
mv .log_ordering_new.dat .log_ordering.dat
sed '/#/d' .parameters/shear_cl/ell.txt > .parameters/shear_cl/ell.new_file.txt
paste .parameters/shear_cl/ell.new_file.txt .parameters/shear_cl/bin_*.new_file.txt > .Cl_out.dat
rm -r .parameters

#### clustering 
# for file in $(ls .parameters/matter_cl/bin_*.txt)
# do
#     sed '/#/d' $file > $file.new_file.txt
#     awk 'NR==1{print}' $file | sed "s/#//" >> .log_ordering.dat
# done
# sed "s/bin_//;s/_/ /" .log_ordering.dat | awk '{print $2, $1}' > .log_ordering_new.dat
# rm .log_ordering.dat
# mv .log_ordering_new.dat .log_ordering.dat
# sed '/#/d' .parameters/matter_cl/ell.txt > .parameters/matter_cl/ell.new_file.txt
# paste .parameters/matter_cl/ell.new_file.txt .parameters/matter_cl/bin_*.new_file.txt > .Cl_out.dat
# rm -r .parameters

#now interpolate the Cl's accordingly

python interpolate_cl.py .Cl_out.dat .Cl_out_interp.dat $2 $3

rm .Cl_out.dat
mv .Cl_out_interp.dat .Cl_out.dat

cp .Cl_out.dat .Cl_fid.dat

# bin the c_ells
python construct_Cell_fits_object.py .Cl_out.dat

