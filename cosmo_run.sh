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



#find out how many redshift bins dndz_filename has
num_cols=$(awk '{print NF}' $1 | sort -nu | tail -n 1)

echo $num_cols

#pad with zeros in front and back of the dndz file
python pad_dndz.py $1

cat <<EOT >> .explanatory.ini
h = $6
T_cmb = 2.726
Omega_b = $8
N_ur = 3.046
Omega_cdm = $4
Omega_dcdmdr = 0.0
Gamma_dcdm = 0.0
N_ncdm = 0
Omega_k = 0.
Omega_Lambda = 0.70
w0_fld = $5
wa_fld = 0.
cs2_fld = 1
YHe = BBN
recombination = RECFAST
reio_parametrization = reio_camb
z_reio = 10.
reionization_exponent = 1.5
reionization_width = 1.5
helium_fullreio_redshift = 3.5
helium_fullreio_width = 0.5
annihilation = 0.
annihilation_variation = 0.
annihilation_z = 1000
annihilation_zmax = 2500
annihilation_zmin = 30
annihilation_f_halo = 20
annihilation_z_halo = 8
on the spot = yes
decay = 0.
output = sCl
number count contributions =
non linear = HALOFIT
modes = s
ic = ad
gauge = synchronous
P_k_ini type = analytic_Pk
k_pivot = 0.05
A_s = $7
n_s = $9
alpha_s = 0.
l_max_lss = $3
selection=yourfile
selection_mean = $(python -c "print ','.join([str(0.1 + 0.01*m) for m in range($num_cols-1)])")
selection_width = $(python -c "print ','.join([str(0.2 + 0.02*m) for m in range($num_cols-1)])")
non_diagonal=$(python -c "print $num_cols - 2")
selection_bias =0
selection_magnification_bias = 0
l_switch_limber_for_nc_local_over_z = 1
l_switch_limber_for_nc_los_over_z = 1
selection_filename=$1
selection_biask_arr=0.0
selection_bias0_arr= ${10}
selection_bias1_arr=0.
root = .output
headers = yes
format = class
write background = no
write thermodynamics = no
write primordial = no
write parameters = yeap
input_verbose = 1
background_verbose = 1
thermodynamics_verbose = 1
perturbations_verbose = 1
transfer_verbose = 1
primordial_verbose = 1
spectra_verbose = 1
nonlinear_verbose = 1
lensing_verbose = 1
output_verbose = 1
EOT


#cleanup the parameter and
#input files for cosmosis
#even if the run fails

{ # this is my bash try block
	./class .explanatory.ini &&
	rm .explanatory.ini
    rm .outputparameters.ini
    rm .outputunused_parameters
} || { # this is catch block
	rm .explanatory.ini
    rm .outputparameters.ini
    rm .outputunused_parameters
}

#This part takes all the calculated cl's from
#cosmosis and pastes them into a convinient
#file. It also outputs a log of the bin ordering
#imposed
awk 'NR==7{print}' .outputcl.dat | grep -o '[0-9]\{1,5\}]' | tr -d ] > .prelim_log_ordering.dat
awk 'NR%2==1{print}' .prelim_log_ordering.dat > .prelim_first_column.dat
awk 'NR%2==0{print}' .prelim_log_ordering.dat > .prelim_second_column.dat
paste .prelim_first_column.dat .prelim_second_column.dat > .log_ordering.dat
rm .prelim_first_column.dat
rm .prelim_second_column.dat
rm .prelim_log_ordering.dat
mv .outputcl.dat .Cl_out.dat
#now interpolate the Cl's accordingly

#divide the prefactor l(l + 1)/(2*pi)
#if third argument is lensing then transform lensing potential powerspectrum
#to kappa power spectrum, if nolensing then don't transform (e.g. if galaxy clustering)
# python divide_prefact.py .outputcl.dat .pre_outputcl.dat lensing
# rm .outputcl.dat
# mv .pre_outputcl.dat .outputcl.dat
#
#
# python interpolate_cl.py .outputcl.dat .Cl_out_interp.dat $2 $3
#
# rm .outputcl.dat
# mv .Cl_out_interp.dat .Cl_out.dat
