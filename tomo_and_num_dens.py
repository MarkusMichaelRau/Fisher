#!/usr/bin/env python

"""
script to take photoz distributions in the format

zmin | zmid | zmax | dndz

and bins the distributions using either 
equal spacing or equal numbers of galaxies

We will use zmid as the support for the distribution
"""

import sys
import numpy as np
from scipy.integrate import simps, trapz, quad

class PhotoZ_Binner(object):
    def __init__(self, z_photoZ, ztype, name, nbins):
        self.load_z_support(z_photoZ, ztype)
        self.load_photoZ(z_photoZ)
        self.norm_dist()
        self.name = name
        self.nbins = nbins

    def load_z_support(self, z_photoZ, ztype):
        ztypes_inds = {'zmin': 0, 'zmid': 1, 'zmax': 2}
        self.z_support = z_photoZ[:, ztypes_inds[ztype]]

    def load_photoZ(self, z_photoZ):
        self.photoZ_dist = z_photoZ[:, -1]

    def norm_dist(self):
        area = simps(self.photoZ_dist, self.z_support)
        self.photoZ_dist = self.photoZ_dist/area

    def calc_equal_size_bins(self):
        # initialize np.array of bins
        cols = self.nbins + 1
        rows = len(self.z_support)
        bins = np.zeros((rows, cols))

        # first column = z_support values
        bins[:, 0] = self.z_support

        # rest = binned pdf values padded with zeros
        
        # get width of bins and number of vals leftover from 
        # integer dividing to ensure we don't forget them
        num_per_bin_nominal = len(self.z_support)//self.nbins
        leftover = len(self.z_support)%self.nbins

        start_in = 0
        end_in = num_per_bin_nominal
        
        # sometimes the number of bins does not cleanly
        # divide the number of data points in the pdf
        # In this case, there are leftover data points
        # which need to be accounted for
        
        # loop through columns and store
        # num_per_bin_nominal (+1 if leftover >0)
        # number of data points in each of them
        for col in range(1, cols):
            # check if there are any leftover data points
            if leftover > 0:
                end_in   += 1
                leftover -= 1
            # save the data in the current bin
            bins[start_in:end_in, col] = self.photoZ_dist[start_in:end_in]
            # increment the start and end indeces of the bins
            start_in  = end_in
            end_in    = start_in + num_per_bin_nominal
        self.bins = bins

    def calc_equal_num_bins(self, int_f=simps):
        total_area = int_f(self.photoZ_dist, self.z_support)
        target_area = total_area/self.nbins
        error_tol = 5e-3
        
        cols = nbins + 1
        rows = len(self.z_support)
        bins = np.zeros((rows, cols))
        
        bins[:, 0] = self.z_support
        
        start_in = 0
        bin_spacing_equal = int(np.floor(len(self.photoZ_dist)/self.nbins))
        end_in = bin_spacing_equal
        last_end_in = end_in
        last_diff = 0
        inds = [0]
        
        for cur_bin in range(nbins-1):
            cur_area = int_f(self.photoZ_dist[start_in: end_in], self.z_support[start_in:end_in])
            while abs(cur_area - target_area) > error_tol:
                if cur_area > target_area:
                    last_end_in = end_in
                    if last_diff >= 0:
                        end_in = end_in - (end_in-start_in)//2
                    elif last_diff < 0:
                        end_in = end_in - (last_end_in-end_in)//2
                    last_diff = 1
                else:
                    last_end_in = end_in
                    if last_diff <= 0:
                        end_in = end_in + (end_in-start_in)//2
                    elif last_diff > 0:
                        end_in = end_in + (end_in - last_end_in)//2
                    last_diff = -1
                cur_area = int_f(self.photoZ_dist[start_in:end_in], self.z_support[start_in: end_in])
            inds.append(end_in)
            start_in = end_in -1
            end_in = start_in + bin_spacing_equal
            
        # now that we have the start indeces of the bins, we
        # make the bins
        for col,_ in enumerate(inds):
            print(col, len(inds))
            if col == 0:
                bins[: inds[col+1], col+1] = self.photoZ_dist[: inds[col+1]]
            elif col < cols -2:
                bins[inds[col]-1:inds[col+1], col+1] = self.photoZ_dist[inds[col]-1:inds[col+1]]
            else:
                bins[inds[col]-1:, col+1] = self.photoZ_dist[inds[col]-1:]

        self.bins = bins
        
    def calc_bins(self, bin_name):
        if bin_name == "eq_size":
            self.calc_equal_size_bins()
        elif bin_name == "eq_num":
            self.calc_equal_num_bins()
        else:
            print("Enter a valid binning scheme")
    
    def get_bins(self):
        return self.bins
    
    def save_bins(self, path=""):
        fname = path + 'tomo_' + self.name + '.dat'
        np.savetxt(fname, self.bins)
        
    def calc_num_dens(self):
        def get_bin_nonzeros(bins, ind):
            # the tomographic bins are organized as a column with
            # zeros at unimportant z-vals and with actual values at
            # important z-vals. 
            # this function picks out the nonzero vals of a given bin
            # 
            # if it is the first bin, it also includes the first value
            # in the bin which may not be included if it is zero 
            if ind == 1:
                if bins[:, ind][0] == 0:
                    inds = np.hstack((np.asarray([0]), (np.nonzero(bins[:, 1])[0])))
                else:
                    inds = np.nonzero(bins[:, ind])[0]
                return bins[:, ind][inds], bins[:, 0][inds]
            else:
                inds = np.nonzero(bins[:, ind])
                return bins[:, ind][inds], bins[:, 0][inds]

        # initialize number of galaxies per steradian
        if self.name == "lensing":
            gal_per_ster = 3.1e8
        elif self.name == "clustering":
            gal_per_ster = 6.5e9

        # Calculate the number density
        # store in a dat file,
        # formatted as comma seperated values

        cols = self.nbins + 1
        num_dens = np.zeros((1, self.nbins))
        for i in range(1, cols):
            y, x = get_bin_nonzeros(self.bins, i)
            area = simps(y, x)
            num_dens[0, i-1] = area * gal_per_ster
        
        # for lensing we need shape noise also
        if self.name == "lensing":
            sig_ep_2 = 0.23
            print("sig_ep^2 = ", sig_ep_2)
            print("n_g = ", num_dens)
            num_dens = 2 * sig_ep_2 / num_dens 
            print("num dens = ", num_dens)

        self.num_dens = num_dens

    def save_num_dens(self, path=""):
        fname = path + 'num_dens_' + self.name + '.dat'
        np.savetxt(fname, self.num_dens)

    def update_nbins(self, nbins):
        self.nbins = nbins
        
if __name__ == "__main__":
    # load args    
    if len(sys.argv) == 1:
        nbins = 2
        name = 'lensing'
        bin_name = 'eq_size'
    elif len(sys.argv) == 4:
        nbins = int(sys.argv[1])
        name = sys.argv[2] 
        bin_name = sys.argv[3]
        if name not in ['lensing', 'clustering']:
            print("name must be lensing or clustering")
            sys.exit(1)
        if bin_name not in ['eq_size', 'eq_num']:
            print("bin_name must be eq_size or eq_num")
            sys.exit(1)
    else:
        print("Incorrect number of arguments")
        sys.exit(1)

    if name == "lensing":
        z_photoZ = np.loadtxt('zdistri_model_z0=1.100000e-01_beta=6.800000e-01_Y10_source')
    elif name == "clustering":
        z_photoZ = np.loadtxt('zdistri_model_z0=2.800000e-01_beta=9.000000e-01_Y10_lens')

    photoZ_binner = PhotoZ_Binner(z_photoZ, 'zmid', name, nbins)
    photoZ_binner.calc_bins(bin_name)
    photoZ_binner.save_bins()
    photoZ_binner.calc_num_dens()
    photoZ_binner.save_num_dens()
