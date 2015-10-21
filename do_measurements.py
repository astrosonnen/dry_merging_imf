import numpy as np
from scipy.optimize import leastsq
import pickle
import sys
import pylab
from matplotlib import rc
rc('text', usetex=True)


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


def do_fits(filename):

    f = open(filename, 'r')
    galaxies = pickle.load(f)
    f.close()

    Ngal = len(galaxies)
    lmstar2 = np.empty(Ngal)
    lreff2 = np.empty(Ngal)
    lsigma2 = np.empty(Ngal)
    laimf2 = np.empty(Ngal)

    lmstar1 = np.empty(Ngal)
    lreff1 = np.empty(Ngal)
    lsigma1 = np.empty(Ngal)
    laimf1 = np.empty(Ngal)

    lmstar0 = np.empty(Ngal)
    lreff0 = np.empty(Ngal)
    lsigma0 = np.empty(Ngal)
    laimf0 = np.empty(Ngal)

    for i in range(0, Ngal):
        lmstar2[i] = np.log10(galaxies[i].mstar_chab[-1])
        lreff2[i] = np.log10(galaxies[i].re[-1])
        lsigma2[i] = np.log10(galaxies[i].veldisp[-1])
        laimf2[i] = np.log10(galaxies[i].aimf[-1])

        lmstar1[i] = np.log10(galaxies[i].mstar_chab[100])
        lreff1[i] = np.log10(galaxies[i].re[100])
        lsigma1[i] = np.log10(galaxies[i].veldisp[100])
        laimf1[i] = np.log10(galaxies[i].aimf[100])

        lmstar0[i] = np.log10(galaxies[i].mstar_chab[0])
        lreff0[i] = np.log10(galaxies[i].re[0])
        lsigma0[i] = np.log10(galaxies[i].veldisp[0])
        laimf0[i] = np.log10(galaxies[i].aimf[0])

    pars_mstar_sigma2 = fit_mstar_sigma_fixed_z(lmstar2, lsigma2, laimf2)
    pars_mstar_sigma1 = fit_mstar_sigma_fixed_z(lmstar1, lsigma1, laimf1)
    pars_mstar_sigma0 = fit_mstar_sigma_fixed_z(lmstar0, lsigma0, laimf0)

    pars_sigma1 = fit_sigma_only(lsigma1, laimf1)
    pars_mstar1 = fit_mstar_only(lmstar1, laimf1)

    pars_sigma0 = fit_sigma_only(lsigma0, laimf0)
    pars_mstar0 = fit_mstar_only(lmstar0, laimf0)

    output = {'z=2': {'mstar-sigma': pars_mstar_sigma2}, \
              'z=1':{'mstar-sigma': pars_mstar_sigma1, 'sigma': pars_sigma1, 'mstar': pars_mstar1}, \
              'z=0':{'mstar-sigma': pars_mstar_sigma0, 'sigma': pars_sigma0, 'mstar': pars_mstar0}}

    return output

mstar_fits = do_fits('mstar_dep_imf_coeff0.5.dat')
vdisp_fits = do_fits('vdisp_dep_imf_coeff2.0.dat')
both_fits = do_fits('mstar-vdisp_dep_imf_coeff0.52.0.dat')

mstar_model_z2 = (mstar_fits['z=2']['mstar-sigma'][0][1], mstar_fits['z=2']['mstar-sigma'][0][2])
vdisp_model_z2 = (vdisp_fits['z=2']['mstar-sigma'][0][1], vdisp_fits['z=2']['mstar-sigma'][0][2])
both_model_z2 = (both_fits['z=2']['mstar-sigma'][0][1], both_fits['z=2']['mstar-sigma'][0][2])

mstar_model_z1 = (mstar_fits['z=1']['mstar-sigma'][0][1], mstar_fits['z=1']['mstar-sigma'][0][2])
vdisp_model_z1 = (vdisp_fits['z=1']['mstar-sigma'][0][1], vdisp_fits['z=1']['mstar-sigma'][0][2])
both_model_z1 = (both_fits['z=1']['mstar-sigma'][0][1], both_fits['z=1']['mstar-sigma'][0][2])

mstar_model_z0 = (mstar_fits['z=0']['mstar-sigma'][0][1], mstar_fits['z=0']['mstar-sigma'][0][2])
vdisp_model_z0 = (vdisp_fits['z=0']['mstar-sigma'][0][1], vdisp_fits['z=0']['mstar-sigma'][0][2])
both_model_z0 = (both_fits['z=0']['mstar-sigma'][0][1], both_fits['z=0']['mstar-sigma'][0][2])

pylab.scatter(mstar_model_z2[0], mstar_model_z2[1], color='r', s=40, label='$M_*$ model, $z=2$')
pylab.scatter(vdisp_model_z2[0], vdisp_model_z2[1], s=40, label='$\sigma$ model, $z=2$')
pylab.scatter(both_model_z2[0], both_model_z2[1], color='g', s=40, label='$M_*-\sigma$ model, $z=2$')

pylab.scatter(mstar_model_z1[0], mstar_model_z1[1], color='r', s=40, label='$M_*$ model, $z=1$', marker='s')
pylab.scatter(vdisp_model_z1[0], vdisp_model_z1[1], s=40, label='$\sigma$ model, $z=1$', marker='s')
pylab.scatter(both_model_z1[0], both_model_z1[1], color='g', s=40, label='$M_*-\sigma$ model, $z=1$', marker='s')

pylab.scatter(mstar_model_z0[0], mstar_model_z0[1], color='r', s=40, label='$M_*$ model, $z=0$', marker='^')
pylab.scatter(vdisp_model_z0[0], vdisp_model_z0[1], s=40, label='$\sigma$ model, $z=0$', marker='^')
pylab.scatter(both_model_z0[0], both_model_z0[1], color='g', s=40, label='$M_*-\sigma$ model, $z=0$', marker='^')
#pylab.legend(scatterpoints=1)
pylab.xlabel('$\partial \log{\\alpha_{\mathrm{IMF}}} / \partial \log{M_*}$', fontsize=16)
pylab.ylabel('$\partial \log{\\alpha_{\mathrm{IMF}}} / \partial \log{\sigma}$', fontsize=16)
pylab.savefig('dep_evolution_plot.png')
pylab.show()


