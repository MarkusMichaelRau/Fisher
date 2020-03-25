import os
import sys
import shutil
import numpy as np
import pyccl as ccl
import get_stencil_deriv_cluster as sd
import get_fisher as gf 
from tomo_and_num_dens import PhotoZ_Binner
from cov_mat import multi_bin_cov, one_bin_cov


class Fisher_Forecaster:
    def __init__(self, probe, bin_type, nbins, deriv_order, derivs_to_calc, \
                 use_binned=True, outdir="out_fisher_cluster/", use_h2=False):
        self.fsky = 0.48
        self.probe = probe
        self.bin_type = bin_type
        self.nbins = nbins
        self.use_binned = use_binned
        self.Nell_bin = [29, 41, 57, 80, 111, 154, 215, 237]
        self.deriv_order = deriv_order
        self.paras = ["om_m", "w0", "h0", "A_s", "om_b", "n_s", "wa"]
        self.use_h2 = use_h2
        if self.use_h2:
            self.paras[0] = "om_m_h2"
            self.paras[4] = "om_b_h2"
        self.outdir = outdir
        if derivs_to_calc == "all":
            self.derivs_to_calc = self.paras
        else:
            self.derivs_to_calc = derivs_to_calc
        self.deltas = [0.15]*len(self.paras)
        self.deltas[-1] = 0.1 # set delta_wa = 0.1 otherwise wa will go below -1
        self.cosmo_params = self.get_cosmo_params()
        self.cosmo = self.get_cosmology(self.cosmo_params)
        self.get_orderings()
        self.get_ells()

    def make_out_dir(self):
        os.makedirs(self.outdir, exist_ok=True)

    def get_tomo_data(self, to_save=True):
        """
        Calls PhotoZ_binner class from tomo_and_num_dens.py
        to bin the photometric redhifts into equal sized or
        equal number bins
        """
        self.photoZ_binner = PhotoZ_Binner('zmid', self.probe, self.nbins, self.bin_type)
        self.photoZ_binner.calc_bins()
        if to_save:
            self.photoZ_binner.save_bins()
        self.photoZ_binner.calc_num_dens()
        if to_save:
            self.photoZ_binner.save_num_dens()
        self.z_dNdz = self.photoZ_binner.bins
        self.num_dens = self.photoZ_binner.num_dens.flatten()
        
    def get_cosmology(self, cosmo_params):
        """
        simple method to make a ccl cosmology object using 
        a dict of cosmological parameters
        """
        if self.use_h2:
            om_c = cosmo_params['om_m_h2']/(cosmo_params['h0']**2)
            om_b = cosmo_params['om_b_h2']/(cosmo_params['h0']**2)
        else:
            om_c = cosmo_params['om_m']
            om_b = cosmo_params['om_b']                     

        return ccl.Cosmology(Omega_c = om_c,
                             Omega_b = om_b, 
                             h = cosmo_params['h0'], 
                             A_s = cosmo_params['A_s'],# * 1e-9, 
                             n_s = cosmo_params['n_s'], 
                             Neff = cosmo_params['N_eff'], 
                             w0= cosmo_params['w0'], 
                             wa= cosmo_params['wa'])

    def get_cosmo_params(self, param_file="fid_values.sh"):
        """
        Load cosmologocial parameters from fid_values.sh file
        """
        cosmo_params = dict()
        with open(param_file, 'r') as f:
            for line in f:
                if line[0] != '#' and line[0] != '\n':
                    words = line.split("=")
                    cosmo_params[words[0]] = float(words[1][:])
        cosmo_params["N_eff"] = 3.046
        if self.use_h2:
            cosmo_params['om_m_h2'] = cosmo_params['om_m']*(cosmo_params['h0']**2)
            cosmo_params['om_b_h2'] = cosmo_params['om_b']*(cosmo_params['h0']**2)
            del(cosmo_params['om_m'])
            del(cosmo_params['om_b'])
        return cosmo_params
    
    def get_orderings(self):
        """
        generate cosmosis-like orderings for auto and cross correlations
        """
        self.orderings = np.array([(i, j) for j in range(1, self.nbins+1) for i in range(1, j+1)])
        
    def get_ells(self):
        """
        make a list of ell values
        """
        self.ells = np.arange(self.cosmo_params['lmin'], self.cosmo_params['lmax']+1)

    def return_pk(self, k, a, cosmo=None, lin=False):
        """
        return matter power spectrum of cosmology object
        defaults to cosmology object with fiducial values
        """
        if cosmo is None:
            cosmo = self.cosmo

        if lin:
            return ccl.linear_matter_power(cosmo, k, a)
        else:
            return ccl.nonlinear_matter_power(cosmo, k, a)


    def get_tracers(self, cosmo):
        """
        return a list of ccl weak lensing tracers to be used to calulate C_ells
        uses the tomographic bins and the cosmology object
        """
        tracers = []
        z = self.z_dNdz[:, 0]
        dNdz = self.z_dNdz[:, 1:]
        for i in range(self.nbins):
            # We need (z, dNdz) as a tuple for CCL:
            z_dNdz = (z, dNdz[:, i])
            
            # We mostly will use the defaults here. By using the defaults, we are assuming 
            # that there are no intrinsic alignments affecting this tracer (not physically realistic
            # but the simplest case for debugging.)
            tracers.append(ccl.WeakLensingTracer(cosmo, dndz=z_dNdz))
        return tracers

    def get_c_ells(self, tracers, cosmo):
        """
        takes a set of tracers and calculates auto and cross correlations using them
        """
        c_ells = []
        for ordering in self.orderings:
            index1 = ordering[0] - 1
            index2 = ordering[1] - 1
            c_ells.append(ccl.angular_cl(cosmo, tracers[index1], tracers[index2], self.ells))
        return np.array([self.ells] + c_ells).T

    def save_c_ells(self, binned=False):
        if binned:
            fname = os.path.join(self.outdir, "Cl_fid_binned.dat")
            np.savetxt(X=self.binned_c_ells, fname=fname)
        else:
            fname = os.path.join(self.outdir, "Cl_fid.dat")
            np.savetxt(X=self.c_ells, fname=fname)
    
    def save_orderings(self):
        fname = os.path.join(self.outdir, "ordering_fid.dat")
        np.savetxt(X=self.orderings, fname=fname, fmt='%i')

    def get_cov_mat(self):
        """
        Calls cov_mat.py to calculate covariance matrices or reuse them if already in place
        """
        if self.nbins == 1:
            if os.path.exists("output_covmat/onebin.mat"):
                print("Reusing old covmats", flush=True)
                self.cov_mats = np.loadtxt("output_covmat/onebin.mat")
            else:
                print("Calculating new covariance matrix", flush=True)
                os.makedirs("output_covmat/", exist_ok=True)
                self.cov_mats = one_bin_cov(self.fsky, self.c_ells, self.num_dens, "output_covmat/")
        else:
            if os.path.exists(os.path.join(self.outdir, "output_covmat/76.0.mat")):
                print("Reusing old covmats", flush=True)
                self.cov_mats = []
                for ell in self.ells:
                    self.cov_mats.append(np.loadtxt(os.path.join(self.outdir, "output_covmat/%s.mat"%str(ell))))
            elif os.path.exists("output_covmat/76.0.mat"):
                print("Reusing old covmats", flush=True)
                self.cov_mats = []
                for ell in self.ells:
                    self.cov_mats.append(np.loadtxt("output_covmat/%s.mat"%str(ell)))
            else:
                print("Calculating new covariance matrix", flush=True)
                cov_out_dir = os.path.join(self.outdir, "output_covmat/")
                os.makedirs(cov_out_dir, exist_ok=True)
                self.cov_mats = multi_bin_cov(self.fsky, self.c_ells, self.orderings, self.num_dens, cov_out_dir)

    def get_binned_ells(self, ells):
        Nell_bin = self.Nell_bin
        binned_ells = [1./Nell_bin[0] * sum(self.ells[0:Nell_bin[0]])] + \
                      [1./Nell_bin[i] * \
                      sum(self.ells[sum(Nell_bin[0:i]):sum(Nell_bin[0:i+1])]) for i in range(1, len(Nell_bin))]
        return np.asarray(binned_ells)

    def get_binned_c_ells(self, c_ells):
        """
        based on Danielle's code to bin the 
        fiducial c_ells 
        """
        # Manually define the number of ell's per bin -
        # this is calculated in ellbins_for_fitscov.ipynb 
        # assuming ell range is 76-999
        Nell_bin = self.Nell_bin
        
        # Create the bins in the C_ells
        binned_c_ells= [[1./Nell_bin[0] * sum(c_ells[0:Nell_bin[0],j])] + \
                        [1./Nell_bin[i] * \
                        sum(c_ells[sum(Nell_bin[0:i]):sum(Nell_bin[0:i+1]),j]) \
                                        for i in range(1, len(Nell_bin))] \
                                        for j in range(1, len(self.orderings)+1)]   
        return np.vstack((self.binned_ells, np.asarray(binned_c_ells))).T
    
    def bin_cov_mats(self):
        """
        based on Danielle's code to bin the 
        theoretical covariance matrices 
        """
        print("Binning the theoretical covariance matrices")
        Nell_bin = self.Nell_bin
        nell = len(self.ells)

        ell_mats_unbinned = self.cov_mats
        ell_mats = [1./Nell_bin[0]**2 * \
                    sum(ell_mats_unbinned[0:Nell_bin[0]])] + \
                   [1./Nell_bin[i]**2 * \
                    sum(ell_mats_unbinned[sum(Nell_bin[0:i]):sum(Nell_bin[0:i+1])]) \
                    for i in range(1, len(Nell_bin))]
        self.binned_cov_mats = ell_mats
        print("Finished binning the covariance matrices.\n")
        cov_dir_binned = os.path.join(self.outdir, "output_covmat_binned/")
        os.makedirs(cov_dir_binned, exist_ok=True)
        for i, x in enumerate(self.binned_cov_mats):
            np.savetxt(fname=os.path.join(cov_dir_binned, str(self.binned_ells[i])+".mat"), X=x)        

    def calc_para_deriv(self, para, order, step_size):
        fid_vals_orig = self.get_cosmo_params()
        lmin = fid_vals_orig['lmin']
        lmax = fid_vals_orig['lmax']
        if para == "wa":
            # since wa has a fiducial value of zero, we use w0 as a reference value
            step = abs(fid_vals_orig["w0"]) * step_size
        else:
            step = abs(fid_vals_orig[para]) * step_size
        coeffs = sd.get_difference_coeffs(order)
        # term_paths = dict()
        step_c_ells = dict()

        # for all possible finite difference coefficients,
        # we calculate a term in the approximation
        for step_mult in coeffs:
            coeff = coeffs[step_mult]
            # if the coefficient is 0 we can ignore the term
            if np.allclose(coeff, 0):
                continue
    
            # print step to the terminal 
            msg = "Changing paramater %s from %e by %d * %e"%(para, fid_vals_orig[para], step_mult, step) 
            # subprocess.call(["echo", msg])
            print(msg, flush=True)

            # load the fiducial values
            fid_vals = self.get_cosmo_params() 
            #step the parameter we want to vary by the given step size
            fid_vals = sd.step_para_val(fid_vals, para, step_mult, step)
            # now we call ccl and get the C_ells
            cosmo = self.get_cosmology(fid_vals)
            tracers = self.get_tracers(cosmo)
            step_c_ells[step_mult] = self.get_c_ells(tracers, cosmo)

            # # organize the outputs so that we can load them later
            # out_Cl_path = get_file_name(".Cl_out.dat", step_mult)
            # term_paths[step_mult] = out_Cl_path
            # subprocess.call(["mv", ".Cl_out.dat", out_Cl_path])
            # subprocess.call(["mv", ".log_ordering.dat", get_file_name(".log_ordering.dat", step_mult)])

        deriv = None
        i = 0
        for step_mult in coeffs:
            if step_mult == 0:
                continue # ignore case where coefficient is 0
            coeff = coeffs[step_mult]
            term = step_c_ells[step_mult]
            if i == 0:
                deriv = np.zeros(term.shape)
                deriv[:, 0] = term[:, 0] # set the ell values
            deriv[:, 1:] += coeff * term[:, 1:]
            i += 1
        deriv[:, 1:] /= step
        print()
        return deriv, step_c_ells

    def get_derivs(self):
        """
        based of the get_stencil_deriv_cluster.py file,
        calculate the n-point stencil derivative
        for all parameters whose derivative we want to 
        calculate (self.derivs_to_calc)
        """
        print("Calculating derivatives", flush=True)

        self.derivs = dict()
        for i in range(len(self.paras)):
            para = self.paras[i]
            delta = self.deltas[i]
            print("Calculating %s derivative"%para, flush=True)
            if para in self.derivs_to_calc:
                self.derivs[para], _ = self.calc_para_deriv(para, self.deriv_order, delta)
            else:
                try:
                    self.derivs[para] = np.loadtxt("deriv_%s_%.5e.dat"%(para, delta))
                except:
                    print("If you don't want to calculate the %s derivative, " + 
                          "then make sure that the derivative exists in a file" + 
                          'with format "deriv_%%s_%%.5e.dat"%(para, delta)', flush=True)
                    sys.exit(1)

    def bin_derivs(self):
        """
        Bins the derivatives after calculating them
        """
        self.binned_derivs = {para:self.get_binned_c_ells(self.derivs[para]) for para in self.derivs_to_calc}

    def save_derivs(self):
        if self.use_binned:
            derivs = self.binned_derivs
        else:
            derivs = self.derivs
        for i in range(len(self.paras)):
            para = self.paras[i]
            delta = self.deltas[i]
            deriv = derivs[para]
            fname = os.path.join(self.outdir, "deriv_%s_%.5e.dat"%(para, delta))
            np.savetxt(X=deriv, fname=fname)

    def get_fisher_mat(self):
        """
        based on the get_fisher.py file, calculates the fisher matrix
        entry by entry
        """
        def get_fisher_comp(para1, para2, ells, derivs, cov_mats):
            if self.nbins == 1:
                return single_bin_fisher(para1, para2, derivs, cov_mats)
            else:
                return multi_bin_fisher(para1, para2, ells, derivs, cov_mats)

        def single_bin_fisher(para1, para2, derivs, cov_mats):
            return derivs[para1].dot(np.linalg.inv(cov_mats)).dot(derivs[para2])

        def multi_bin_fisher(para1, para2, ells, derivs, cov_mats):
            fisher_comp = 0.0
            for idx, _ in enumerate(ells):
                left_vec = derivs[para1][idx, 1:]
                right_vec = derivs[para2][idx, 1:]
                fisher_comp += left_vec.dot(np.linalg.inv(cov_mats[idx]).dot(right_vec))
            return fisher_comp

        # initialize empty array
        n = len(self.paras)
        fisher = np.zeros((n, n))

        # choose binned vs non_binned quantities
        if self.use_binned:
            ells = self.binned_ells
            c_ells = self.binned_c_ells
            derivs = self.binned_derivs
            cov_mats = self.binned_cov_mats
        else:
            ells = self.ells
            c_ells = self.c_ells
            derivs = self.derivs
            cov_mats = self.cov_mats
        # calc fisher components
        for i in range(n):
            for j in range(n):
                fisher[i, j] = get_fisher_comp(self.paras[i], self.paras[j], ells, derivs, cov_mats)
        self.fisher = fisher

    def add_priors(self):
        """
        Adds priors based on Planck data to fisher matrix
        """
        try:
            # load priors for params
            sigma_om_m = 0           #don't set prior
            sigma_w0 = 0             #don't set prior
            sigma_h0 = 0.0057        #set prior
            sigma_A_s = 0            #don't set prior
            sigma_om_b = 0.00066     #set prior
            sigma_n_s = 0            #don't set prior
            sigma_galbias = 0        #don't set prior
            sigma_wa = 0             #don't set prior

            # make a matrix with error on diagonal
            if self.probe == "lensing":
                priors = np.diag([sigma_om_m, sigma_w0, sigma_h0, sigma_A_s,
                                  sigma_om_b, sigma_n_s, sigma_wa])
            elif self.probe == "clustering":
                priors = np.diag([sigma_om_m, sigma_w0, sigma_h0, sigma_A_s,
                                  sigma_om_b, sigma_n_s, sigma_galbias, sigma_wa])
            
            # get corresponding fisher matrix by squaring to get a variance
            # matrix, and then inverting non-zero values
            priors[np.nonzero(priors)] = priors[np.nonzero(priors)]**(-2.)

            self.fisher_with_priors = self.fisher + priors
        except:
            print("Couldn't add priors to Fisher matrix. Make sure that"+
                  " the fisher matrix was calculated properly")

    def save_fisher(self, priors=False):
        if priors:
            fname = os.path.join(self.outdir, "fisher_out_with_priors.dat")
            np.savetxt(X=self.fisher_with_priors, fname=fname)
        else:
            fname = os.path.join(self.outdir, "fisher_out.dat")
            np.savetxt(X=self.fisher, fname=fname)

    def get_fom(self, priors=True, para_pairs_list=None):
        """
        calulates the DETF FoM and the om_m A_s
        Figure of Merit
        """
        def calc_fom(cov):
            # takes marg cov matrix
            # inverts it to get reduced fisher matrix
            # returns sqrt of det of fisher matrix
            fisher = np.linalg.inv(cov)
            FOM = np.sqrt(np.linalg.det(fisher))
            return FOM

        def marg_cov(cov, para1, para2):
            ind1 = self.paras.index(para1)
            ind2 = self.paras.index(para2)
            inds = np.asarray((ind1, ind2))
            # first get all indeces to remove
            all_cov_inds = range(len(cov[0]))
            del_inds = np.delete(all_cov_inds, inds)
            # then delete rows and cols
            marged_cov = np.delete(np.delete(cov, del_inds, 0), del_inds, 1)
            return marged_cov

        if priors:
            fisher = self.fisher_with_priors
            self.foms_with_priors = dict()
        else:
            fisher = self.fisher
            self.foms = dict()
        
        if para_pairs_list is None:
            para_pairs_list = [["om_m", "A_s"], ["w0", "wa"]]
            if self.use_h2:
                para_pairs_list[0][0] = "om_m_h2"

        for para1, para2 in para_pairs_list:
            cov = np.linalg.inv(fisher)
            marged_cov = marg_cov(cov, para1, para2)
            fom = calc_fom(marged_cov)
            if priors:
                self.foms_with_priors[(para1, para2)] = fom
            else:
                self.foms[(para1, para2)] = fom
            print(para1, para2, fom, flush=True)

    def save_fom(self, priors=True):
        if priors:
            foms = self.foms_with_priors
            suffix = "_priors.dat"
        else: 
            foms = self.foms
            suffix = ".dat"
        for para1, para2 in foms:
            fom = np.array([foms[(para1, para2)]])
            fname = os.path.join(self.outdir, "fom_%s_%s"%(para1, para2)+suffix)
            np.savetxt(X=fom, fname=fname)

    #TODO make cleanup function
    def cleanup_cov_mats(self):
        os.makedirs(self.outdir, exist_ok=True)
        covmat_dir = "output_covmat/"
        if os.path.exists(covmat_dir):
            shutil.move(covmat_dir, os.path.join(self.outdir, covmat_dir))

    #TODO tuning derivatives 
    def get_derivative_tunings(self, tune_para):
        # initialize tuning vars
        delta_start = 1e-3
        delta_end = 3e-1
        delta_num = 50
        deltas = np.geom_space(delta_start, delta_end, delta_num)

        # initialize fom storage
        fom_wo_wa = np.zeros(deltas.shape)
        fom_wo_wa_priors = np.zeros(deltas.shape)
        fom_om_m_A_s = np.zeros(deltas.shape)
        fom_om_m_A_s_priors = np.zeros(deltas.shape)

        for i in range(len(deltas)):
            self.derivs[para], _ = self.calc_para_deriv(tune_para, self.deriv_order, deltas[i])
            self.get_fisher_mat()
            self.add_priors()
            self.get_fom(priors=False)
            self.get_fom(priors=True)
            fom_w0_wa[i] = self.foms[('w0', 'wa')]
            fom_om_m_A_s = self.foms[('om_m', 'A_s')]
            # TODO


    def get_fiducial_data(self):
        # calc fiducial stuff
        self.make_out_dir()
        self.get_tomo_data()
        print("Getting fiducial C_ells", flush=True)
        self.tracers = self.get_tracers(self.cosmo)
        self.c_ells = self.get_c_ells(self.tracers, self.cosmo)
        self.save_c_ells()
        self.save_orderings()
        print("Fiducial C_ells and orderings saved in Cl_fid.dat and ordering_fid.dat\n", flush=True)
        if self.use_binned:
            print("Binning the ells and the C_ells...")
            self.binned_ells = self.get_binned_ells(self.ells)
            self.binned_c_ells = self.get_binned_c_ells(self.c_ells)
            self.save_c_ells(binned=self.use_binned)
            print("Finished binning the ells and C_ells")
        print()

    def get_cov_and_deriv_data(self):
        # calc covariance matrix and derivs
        self.get_cov_mat()
        self.get_derivs()
        if self.use_binned:
            self.bin_cov_mats()
            self.bin_derivs()
        self.save_derivs()
        print("Done calculating derivatives\n", flush=True)

    def get_fisher_and_fom_data(self):
        # calc fisher matrix
        print("Caclulating Fisher Matrix", flush=True)
        self.get_fisher_mat()
        self.add_priors()
        print(self.fisher_with_priors)
        self.save_fisher()
        self.save_fisher(priors=True)
        print("Done. Fisher matrix has been saved to fisher_out.dat and fisher_out_with_priors.dat\n", flush=True)
        print("Calculating DETF Figure of Merit", flush=True)
        self.get_fom(priors=False)
        self.save_fom(priors=False)
        self.get_fom(priors=True)
        self.save_fom(priors=True)
        print("Done calculating Figure of Merit", flush=True)

    def run_fisher_pipeline(self):
        self.get_fiducial_data()        
        self.get_cov_and_deriv_data()
        self.get_fisher_and_fom_data()

if __name__=="__main__":
    probe = "lensing"
    bin_type = "equal_size"
    nbins = 5
    deriv_order = 2
    derivs_to_calc = "all"
    use_binned = True
    F = Fisher_Forecaster(probe, bin_type, nbins, deriv_order, derivs_to_calc, use_binned=use_binned)
    F.run_fisher_pipeline()
    print("Finished running successfully")
