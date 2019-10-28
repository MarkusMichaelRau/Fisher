import os
import sys
import numpy as np
import pyccl as ccl
import get_stencil_deriv_cluster as sd
import get_fisher as gf 
from tomo_and_num_dens import PhotoZ_Binner
from cov_mat import multi_bin_cov, one_bin_cov


class Fisher_Forecaster:
    def __init__(self, probe, bin_type, nbins, deriv_order, derivs_to_calc):
        self.fsky = 0.48
        self.probe = probe
        self.bin_type = bin_type
        self.nbins = nbins
        self.deriv_order = deriv_order
        self.paras = ["om_m", "w0", "h0", "A_s", "om_b", "n_s", "wa"]
        if derivs_to_calc == "all":
            self.derivs_to_calc = self.paras
        else:
            self.derivs_to_calc = derivs_to_calc
        self.deltas = [0.1]*len(self.paras)
        self.cosmo_params = self.get_cosmo_params()
        self.cosmo = self.get_cosmology(self.cosmo_params)
        self.get_orderings()
        self.get_ells()

    def get_tomo_data(self):
        """
        Calls PhotoZ_binner class from tomo_and_num_dens.py
        to bin the photometric redhifts into equal sized or
        equal number bins
        """
        self.photoZ_binner = PhotoZ_Binner('zmid', self.probe, self.nbins, self.bin_type)
        self.photoZ_binner.calc_bins()
        self.photoZ_binner.save_bins()
        self.photoZ_binner.calc_num_dens()
        self.photoZ_binner.save_num_dens()
        self.z_dNdz = self.photoZ_binner.bins
        self.num_dens = self.photoZ_binner.num_dens.flatten()
        
    def get_cosmology(self, cosmo_params):
        """
        simple method to make a ccl cosmology object using 
        a dict of cosmological parameters
        """
        return ccl.Cosmology(Omega_c = cosmo_params['om_m'],
                             Omega_b = cosmo_params['om_b'], 
                             h = cosmo_params['h0'], 
                             A_s = cosmo_params['A_s'], 
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
        return cosmo_params
    
    def get_orderings(self):
        """
        generate cosmosis-like orderings for auto and cross correlations
        """
        self.orderings = np.array([(i, j) for j in range(1, nbins+1) for i in range(1, j+1)])
        
    def get_ells(self):
        """
        make a list of ell values
        """
        self.ells = np.arange(self.cosmo_params['lmin'], self.cosmo_params['lmax']+1)

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

    def save_c_ells(self):
        np.savetxt(X=self.c_ells, fname="Cl_fid.dat")
    
    def save_orderings(self):
        np.savetxt(X=self.orderings, fname="ordering_fid.dat", fmt='%i')

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
            if os.path.exists("output_covmat/76.0.mat"):
                print("Reusing old covmats", flush=True)
                self.cov_mats = []
                for ell in self.ells:
                    self.cov_mats.append(np.loadtxt("output_covmat/%s.mat"%str(ell)))
            else:
                print("Calculating new covariance matrix", flush=True)
                os.makedirs("output_covmat/", exist_ok=True)
                self.cov_mats = multi_bin_cov(self.fsky, self.c_ells, self.orderings, self.num_dens, "output_covmat/")

    #TODO bin c_ells and cov mats

    def get_derivs(self):
        """
        based of the get_stencil_deriv_cluster.py file,
        calculate the n-point stencil derivative
        for all parameters whose derivative we want to 
        calculate (self.derivs_to_calc)
        """
        print("Calculating derivatives", flush=True)

        def calc_para_deriv(para, order, step_size):
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
            return deriv

        self.derivs = dict()
        for i in range(len(self.paras)):
            para = self.paras[i]
            delta = self.deltas[i]
            print("Calculating %s derivative"%para, flush=True)
            if para in self.derivs_to_calc:
                self.derivs[para] = calc_para_deriv(para, self.deriv_order, delta)
            else:
                try:
                    self.derivs[para] = np.loadtxt("deriv_%s_%.5e.dat"%(para, delta))
                except:
                    print("If you don't want to calculate the %s derivative, " + 
                          "then make sure that the derivative exists in a file" + 
                          'with format "deriv_%%s_%%.5e.dat"%(para, delta)', flush=True)
                    sys.exit(1)

    def save_derivs(self):
        for i in range(len(self.paras)):
            para = self.paras[i]
            delta = self.deltas[i]
            deriv = self.derivs[para]
            np.savetxt(X=deriv, fname="deriv_%s_%.5e.dat"%(para, delta))

    def get_fisher_mat(self):
        """
        based on the get_fisher.py file, calculates the fisher matrix
        """
        def get_fisher_comp(para1, para2):
            if nbins == 1:
                return single_bin_fisher(para1, para2)
            else:
                return multi_bin_fisher(para1, para2)

        def single_bin_fisher(para1, para2):
            return self.derivs[para1].dot(np.linalg.inv(self.cov_mats)).dot(self.derivs[para2])

        def multi_bin_fisher(para1, para2):
            fisher = 0.0
            for idx, ell in enumerate(self.ells):
                left_vec = self.derivs[para1][idx, 1:]
                right_vec = self.derivs[para2][idx, 1:]
                fisher += left_vec.dot(np.linalg.inv(self.cov_mats[idx]).dot(right_vec))
            return fisher

        n = len(self.paras)
        fisher = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                fisher[i, j] = get_fisher_comp(self.paras[i], self.paras[j])
        self.fisher = fisher

    def save_fisher(self):
        np.savetxt(X=self.fisher, fname="fisher_out.dat")

    def get_fom(self):
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

        self.foms = dict()
        for para1, para2 in [("om_m", "A_s"), ("w0", "wa")]:
            cov = np.linalg.inv(self.fisher)
            marged_cov = marg_cov(cov, para1, para2)
            fom = calc_fom(marged_cov)
            self.foms[(para1, para2)] = fom
            print(para1, para2, fom, flush=True)

    def save_fom(self):
        for para1, para2 in self.foms:
            fom = np.array([self.foms[(para1, para2)]])
            np.savetxt(X=fom, fname="fom_%s_%s.dat"%(para1, para2))

    def run_fisher_pipeline(self):
        # calc fiducial stuff
        self.get_tomo_data()
        print("Getting fiducial C_ells", flush=True)
        self.tracers = self.get_tracers(self.cosmo)
        self.c_ells = self.get_c_ells(self.tracers, self.cosmo)
        self.save_c_ells()
        self.save_orderings()
        print("Fiducial C_ells and orderings saved in Cl_fid.dat and ordering_fid.dat\n", flush=True)

        # calc covariance matrix and derivs
        self.get_cov_mat()
        self.get_derivs()
        self.save_derivs()
        print("Done calculating derivatives\n", flush=True)

        # calc fisher matrix
        print("Caclulating Fisher Matrix", flush=True)
        self.get_fisher_mat()
        print(self.fisher)
        self.save_fisher()
        print("Done. Fisher matrix has been saved to fisher_out.dat\n", flush=True)
        print("Calculating DETF Figure of Merit", flush=True)
        self.get_fom()
        self.save_fom()
        print("Done calculating Figure of Merit", flush=True)

if __name__=="__main__":
    probe = "lensing"
    bin_type = "equal_size"
    nbins = 5
    deriv_order = 2
    derivs_to_calc = "all"
    F = Fisher_Forecaster(probe, bin_type, nbins, deriv_order, derivs_to_calc)
    F.run_fisher_pipeline()
