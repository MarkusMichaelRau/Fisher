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

cat <<EOT >> .values.dat
[cosmological_parameters]
omega_m = $4
h0 = $6
omega_b = $8
tau = 0.08
n_s = $9
A_s = $7
w = $5
wa = ${11}
[bias_parameters]
c=1.0
b0=${10}
alpha=1.0
EOT

cat <<EOT >> .run_cosmosis.ini
[runtime]
sampler = test
root = /home/markus/Arbeitsfläche/final_analysis/forecast/cosmosis
[test]
save_dir=.parameters
fatal_errors=T
[pipeline]
modules = consistency camb halofit extrapolate_power growth_factor clerkin load_nz shear_shear
values = .values.dat
likelihoods =
extra_output =
quiet=F
timing=T
debug=F
[consistency]
file = cosmosis-standard-library/utility/consistency/consistency_interface.py
[camb]
file = cosmosis-standard-library/boltzmann/camb/camb.so
mode=all
lmax=2500
feedback=0
[halofit]
file = cosmosis-standard-library/boltzmann/halofit/halofit_module.so
[extrapolate_power]
file=cosmosis-standard-library/boltzmann/extrapolate/extrapolate_power.py
kmax=500.0
[growth_factor]
file=cosmosis-standard-library/structure/growth_factor/interface.so
[clerkin]
file=cosmosis-standard-library/bias/clerkin/clerkin_interface.py
mode='both'
model='gtd'
[load_nz]
file = cosmosis-standard-library/number_density/load_nz/load_nz.py
filepath = $1
[shear_shear]
file = cosmosis-standard-library/shear/spectra/interface.so
ell_min = $(python -c "print float($2)")
ell_max = $(python -c "print float($3)")
n_ell = 800
galaxy_bias=1.0
matter_spectra=T
ggl_spectra=T
EOT

#cleanup the parameter and
#input files for cosmosis
#even if the run fails
source ./cosmosis/setup-my-cosmosis
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

for file in $(ls .parameters/matter_cl/bin_*.txt)
do
    sed '/#/d' $file > $file.new_file.txt
    awk 'NR==1{print}' $file | sed "s/#//" >> .log_ordering.dat
done
sed "s/bin_//;s/_/ /" .log_ordering.dat | awk '{print $2, $1}' > .log_ordering_new.dat
rm .log_ordering.dat
mv .log_ordering_new.dat .log_ordering.dat
sed '/#/d' .parameters/matter_cl/ell.txt > .parameters/matter_cl/ell.new_file.txt
paste .parameters/matter_cl/ell.new_file.txt .parameters/matter_cl/bin_*.new_file.txt > .Cl_out.dat
rm -r .parameters

#now interpolate the Cl's accordingly

python interpolate_cl.py .Cl_out.dat .Cl_out_interp.dat $2 $3

rm .Cl_out.dat
mv .Cl_out_interp.dat .Cl_out.dat
