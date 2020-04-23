import sys
import numpy as np

def calc_fom(cov):
    # takes marg cov matrix
    # inverts it to get reduced fisher matrix
    # returns sqrt of det of fisher matrix
    fisher = np.linalg.inv(cov)
    FOM = np.sqrt(np.linalg.det(fisher))
    return FOM

def marg_cov(cov, paras, para1, para2):
    ind1 = paras.index(para1)
    ind2 = paras.index(para2)
    inds = np.asarray((ind1, ind2))
    # first get all indeces to remove
    all_cov_inds = range(len(cov[0]))
    del_inds = np.delete(all_cov_inds, inds)
    # then delete rows and cols
    marged_cov = np.delete(np.delete(cov, del_inds, 0), del_inds, 1)
    return marged_cov

def get_paras_fom(fisher, paras, para1, para2): 
    cov = np.linalg.inv(fisher)
    marged_cov = marg_cov(cov, paras, para1, para2)
    fom = calc_fom(marged_cov)
    return fom

def main():
    fisher_path = sys.argv[1]
    fisher = np.loadtxt(fisher_path)
    paras = ["om_m", "w0", "h0", "A_s", "om_b", "n_s", "wa"]
    for para1, para2 in [("om_m", "A_s"), ("w0", "wa")]:
        fom = get_paras_fom(fishers, paras, para1, para2)
        print(para1, para2, fom)

if __name__ == "__main__":
    main()
