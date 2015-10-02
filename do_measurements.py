import numpy as np
from scipy.optimize import leastsq
import pickle
import sys


def fit_mstar_re_fixed_z(lmstar_sample, lreff_sample, aimf_sample, guess=(0.3, 0.3, 0.)):

    def modelfunc(p):
        return p[0] + p[1]*(lmstar_sample - 11.5) + p[2]*(lreff_sample - np.log10(3.))

    def errfunc(p):
        return modelfunc(p) - aimf_sample

    par, cov = leastsq(errfunc, guess)
    scat = (sum(errfunc(par)**2)/float(len(aimf_sample)))**0.5

    return par, scat


def fit_mstar_sigma_fixed_z(lmstar_sample, lsigma_sample, aimf_sample, guess=(0.3, 0.3, 0.)):

    def modelfunc(p):
        return p[0] + p[1]*(lmstar_sample - 11.5) + p[2]*(lsigma_sample - np.log10(200.))

    def errfunc(p):
        return modelfunc(p) - aimf_sample

    par, cov = leastsq(errfunc, guess)
    scat = (sum(errfunc(par)**2)/float(len(aimf_sample)))**0.5

    return par, scat


def fit_sigma_only(lsigma_sample, aimf_sample, guess=(1.3, 0.)):

    def modelfunc(p):
        return p[0] + p[1]*(lsigma_sample - 2.3)

    def errfunc(p):
        return modelfunc(p) - aimf_sample

    par, cov = leastsq(errfunc, guess)
    scat = (sum(errfunc(par)**2)/float(len(aimf_sample)))**0.5

    return par, scat


def fit_mstar_only(lmstar_sample, aimf_sample, guess=(0.3, 0.)):

    def modelfunc(p):
        return p[0] + p[1]*(lmstar_sample - 11.)

    def errfunc(p):
        return modelfunc(p) - aimf_sample

    par, cov = leastsq(errfunc, guess)
    scat = (sum(errfunc(par)**2)/float(len(aimf_sample)))**0.5

    return par, scat


if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'mstar_dep_imf_coeff0.5.dat'

f = open(filename, 'r')
galaxies = pickle.load(f)
f.close()

Ngal = len(galaxies)
lmstar2 = np.empty(Ngal)
lreff2 = np.empty(Ngal)
lsigma2 = np.empty(Ngal)
laimf2 = np.empty(Ngal)

lmstar0 = np.empty(Ngal)
lreff0 = np.empty(Ngal)
lsigma0 = np.empty(Ngal)
laimf0 = np.empty(Ngal)

for i in range(0, Ngal):
    lmstar2[i] = np.log10(galaxies[i].mstar_chab[-1])
    lreff2[i] = np.log10(galaxies[i].re[-1])
    lsigma2[i] = np.log10(galaxies[i].veldisp[-1])
    laimf2[i] = np.log10(galaxies[i].aimf[-1])

    lmstar0[i] = np.log10(galaxies[i].mstar_chab[0])
    lreff0[i] = np.log10(galaxies[i].re[0])
    lsigma0[i] = np.log10(galaxies[i].veldisp[0])
    laimf0[i] = np.log10(galaxies[i].aimf[0])

pars_re2 = fit_mstar_re_fixed_z(lmstar2, lreff2, laimf2)
pars_re0 = fit_mstar_re_fixed_z(lmstar0, lreff0, laimf0)

pars_mstar_sigma2 = fit_mstar_sigma_fixed_z(lmstar2, lsigma2, laimf2)
pars_mstar_sigma0 = fit_mstar_sigma_fixed_z(lmstar0, lsigma0, laimf0)

pars_sigma0 = fit_sigma_only(lsigma0, laimf0)
pars_mstar0 = fit_mstar_only(lmstar0, laimf0)

print 'fit of dependence on mstar-sigma, z=2:', pars_mstar_sigma2
print 'fit of dependence on mstar-sigma, z=0:', pars_mstar_sigma0
print 'fit of dependence on sigma, z=0:', pars_sigma0
print 'fit of dependence on mstar, z=0:', pars_mstar0

import pylab
pylab.scatter(lmstar2, laimf2)
pylab.scatter(lmstar0, laimf0, color='r')
pylab.xlabel('$\log{M_*}$')
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$')
pylab.show()

pylab.scatter(lsigma2, laimf2)
pylab.scatter(lsigma0, laimf0, color='r')
pylab.xlabel('$\sigma_*$')
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$')
pylab.show()

