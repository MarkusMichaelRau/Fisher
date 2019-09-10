from subprocess import call
import numpy as np
from numpy.linalg import inv
from scipy.stats import norm
import os


class FisherForecast(object):
    def __init__(self):
        abspath = os.path.abspath(__file__)
        dname = os.path.dirname(abspath)
        os.chdir(dname)
        self.base_fname = './'

    def get_forecast(self, biased_z_distribution, unbiased_z_distribution,
                         number_dens=118200000.0, shape_noise=0.23,
                         cleanup=True):
        ###This is a fix to make the fisher forecast consistent despite a
        ###bug in the CLASS
        biased_z_distribution, unbiased_z_distribution = _order_distributions(biased_z_distribution,
                                                                                   unbiased_z_distribution)

        fname_biased_z_distr = self.base_fname+'biased_z_distribution.dat'
        fname_unbiased_z_distr = self.base_fname+'unbiased_z_distribution.dat'
        np.savetxt(X=biased_z_distribution, fname=fname_biased_z_distr)
        np.savetxt(X=unbiased_z_distribution, fname=fname_unbiased_z_distr)
        if shape_noise != None:
            number_dens = (number_dens * 2.)/shape_noise**2
        np.savetxt(X=number_dens, fname='.num_dens.dat')
        call([self.base_fname+"get_bias_cluster.sh", fname_biased_z_distr, fname_unbiased_z_distr, str(number_dens)])
        fisher = np.loadtxt(self.base_fname+'out_fisher_cluster/fisher_out.dat')
        para_bias = np.loadtxt(self.base_fname+'out_fisher_cluster/parameter_bias.dat')
        call(['mv', fname_biased_z_distr, 'out_fisher_cluster/'+fname_biased_z_distr])
        call(['mv', fname_unbiased_z_distr, 'out_fisher_cluster/'+fname_unbiased_z_distr])
        if cleanup is True:
            self.cleanup()
        return para_bias, fisher


    def cleanup(self):
        call(['rm', '-r', self.base_fname+'out_fisher_cluster'])
        os.chdir('../')

def _order_distributions(biased_z_distribution, unbiased_z_distribution):
    means = np.zeros((unbiased_z_distribution.shape[1] - 1,))
    for i in range(1, unbiased_z_distribution.shape[1]):
        means[i - 1] = _get_mean(unbiased_z_distribution[:, 0],
                                      unbiased_z_distribution[:, i])
    idx_sorted = np.argsort(means) + 1 #because z is first column
    idx_sorted = np.insert(idx_sorted, 0, 0)
    return biased_z_distribution[:, idx_sorted], unbiased_z_distribution[:, idx_sorted]

def _get_mean(z, pdf):
    pdf = pdf/np.trapz(pdf, z)
    return np.trapz(pdf*z, z)


def test_ordering():
    from scipy.stats import norm
    rv = norm()
    z = np.linspace(0.0, 5.0, num=1000)
    pdf1_unbiased = np.column_stack((z, rv.pdf(z - 1), rv.pdf(z - 2)))
    pdf1_biased = np.column_stack((z, rv.pdf(z - 1.2), rv.pdf(z - 0.8)))

    pdf2_unbiased = np.column_stack((z, rv.pdf(z - 2), rv.pdf(z - 1)))
    pdf2_biased = np.column_stack((z, rv.pdf(z - 1.2), rv.pdf(z - 0.8)))

    test1 = _order_distributions(pdf1_biased, pdf1_unbiased)
    print test1[-1]
    print test1[0]
    print test1[1]

    test2 = _order_distributions(pdf2_biased, pdf2_unbiased)
    print test2[-1]
    print test2[0]
    print test2[1]

def test_1bin():
    x = np.linspace(0.0, 3.0, num=1000)
    norm_pdf = norm()
    unbiased_z_1 = np.column_stack((x, norm_pdf.pdf(x - 3.), norm_pdf.pdf(x - 3.4)))
    biased_z_1 = np.column_stack((x, norm_pdf.pdf(x - 3.01), norm_pdf.pdf(x - 3.41)))

    # unbiased_z_2 = np.column_stack((x, norm_pdf.pdf(x - 3.4), norm_pdf.pdf(x - 3.)))
    # biased_z_2 = np.column_stack((x, norm_pdf.pdf(x - 3.41), norm_pdf.pdf(x - 3.01)))

    model = FisherForecast()
    result = model.get_forecast(unbiased_z_1, biased_z_1,
    number_dens=np.array([118200000.0, 18200000.0]), cleanup=False)
    #
    # model2 = FisherForecast()
    # result_new = model2.get_forecast(unbiased_z_2, biased_z_2,
    # number_dens=np.array([18200000.0, 118200000.0]), cleanup=False)

    # print 'numerics error bias'
    # print result[0]
    # print result_new[0]
    # print result[0] - result_new[0]
    #
    # print "Fisher error"
    # print result[1]
    # print result_new[1]
    # print result[1] - result_new[1]
    # # #error_parameters = np.sqrt(np.diag(inv(result[1])))
    # # return result[0], error_parameters


if (__name__ == '__main__'):
    test_ordering()
    # #print('single bin case')
    # test_1bin()
    # #
    # # print('multibin case')
    # # result_singlebin = test_1bin()
    # # print result_singlebin
    # # call(['pwd'])
    # # print('multibin case')
    # # result_multibin = test_2bin()
    # # print result_multibin
    # # if np.isclose(np.sum(result_singlebin[0]), 0.) and np.isclose(np.sum(result_multibin[0]), 0.) \
    # # and np.isclose(result_multibin[1][0], 1.02117711e-02) and np.isclose(result_multibin[1][1], 1.48107175e-10) :
    # #     print 'test passed'
