# fisher
## Naren's branch

## To Do:

1. Reproduce Kitching Paper Results
2. Extend to DETF FoM
3. Consider non gaussian covariances

## Reproducing Kitching Paper

The goal here is to show that using the same FoM as Kitching,
we get the similar results, albeit using LSST redshift 
distributions as opposed to Euclid like they do.

To do so, we will use the Y10 LSST redshift distributions from 
the LSST DESC Science Requirements Document v1 (LSST DESC 2018)
([paper](https://arxiv.org/pdf/1809.01669.pdf)).

Since the Kitching paper is lensing only, our work will also
be lensing only - so we will use the source redshift distributions.

In addition, they use a figure of merit defined by the
cosmic shear signal to noise ratio for the dark energy parameter
$*w_0*$ given by 

$$
F = - \sum_{\ell, i, j}\left[\frac{\partial C_\ell (\eta_i, \eta_j)}{\partial \w_0} \right]^2 
\frac{1}{[C_\ell(\eta_i, \eta_j) + N_\ell(\eta_i, \eta_j)]^2}
$$

We will start by only considering auto correlations since the 
cross correlations for lensing are small compared to the autocorrelations.