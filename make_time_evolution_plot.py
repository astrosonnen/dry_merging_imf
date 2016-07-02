import pickle
import pylab
import numpy as np
import do_measurements as dm
from plotters import rgb_alpha
from sonnentools import cgsconstants as cgs
from scipy.interpolate import splev
from scipy.optimize import leastsq
import pymc


mpiv = 11.5
vpiv = 2.4

f = open('/gdrive/Python/vdmodel_2013/Kvir_spline.dat', 'r')
kvir_spline = pickle.load(f)
f.close()

f = open('spiniello_table.txt', 'r')
spin_table = np.loadtxt(f)
f.close()
spin_z = 0.025
spin_aimf = spin_table[1, :]
spin_lsigma = spin_table[0, :]
spin_lsigma_err = spin_table[2, :]
spin_aimf_err = spin_table[3, :]

f = open('conroy2012_table.txt', 'r')
conroy_table = np.loadtxt(f, usecols=(5, 10, 11))
f.close()
conroy_z = 0.01
conroy_sigma = conroy_table[:, 0]
conroy_ml = conroy_table[:, 1]
conroy_mlmw = conroy_table[:, 2]
conroy_lsigma = np.log10(conroy_sigma)
conroy_lsigma_err = 0.*conroy_lsigma + 0.05

conroy_aimf = conroy_ml/conroy_mlmw
conroy_aimf = np.log10(conroy_aimf) - 0.2

conroy_aimf_err = 0.*conroy_aimf + 0.1

def fit_aimf_vs_sigma(aimf, aimf_err, sigma, sigma_err):
    apar = pymc.Uniform('a', lower=0., upper=2., value=1.)
    bpar = pymc.Uniform('b', lower=-1., upper=1., value=0.)
    spar = pymc.Uniform('scat', lower=0., upper=1., value=0.1)

    sigma_mu = pymc.Uniform('sigma_mu', lower=2., upper=3., value=sigma.mean())
    sigma_sig = pymc.Uniform('sigma_sig', lower=0., upper=1., value=sigma.std())

    pars = [apar, bpar, spar, sigma_mu, sigma_sig]

    nobj = len(aimf)

    @pymc.deterministic
    def like(p=pars):
        a, b, scat, sigma_mu, sigma_sig = p

        sumlogp = 0.

        for i in range(nobj):

            alpha = b
            E = alpha + a*sigma_sig**2/(sigma_sig**2 + sigma_err[i]**2)*(sigma[i] - sigma_mu)
            Var = a**2*sigma_sig**2 + scat**2 + aimf_err[i]**2 - (a*sigma_sig**2)**2/(sigma_sig**2 + sigma_err[i]**2)

            logpy = -0.5*(aimf[i] - E)**2/Var - 0.5*np.log(2.*np.pi*Var)

            varx = sigma_err[i]**2 + sigma_sig**2

            logpx = -0.5*(sigma[i] - sigma_mu)**2/varx - 0.5*np.log(2.*np.pi*varx)

            sumlogp += logpy + logpx

        return sumlogp

    @pymc.stochastic
    def logp(p=pars, value=0., observed=True):
        return like

    M = pymc.MCMC(pars)
    M.use_step_method(pymc.AdaptiveMetropolis, pars)
    M.sample(11000, 1000)

    chain = {}
    for par in pars:
        chain[str(par)] = M.trace(str(par))[:]

    return chain

spin_fit = fit_aimf_vs_sigma(spin_aimf, spin_aimf_err, spin_lsigma, spin_lsigma_err)
conroy_fit = fit_aimf_vs_sigma(conroy_aimf, conroy_aimf_err, conroy_lsigma, conroy_lsigma_err)

fsize = 16

bandcolor = rgb_alpha((0,255,255), 1.)
bandcolor = (0, int(round(0.7*255)), 255)

f = open('pop_mstar_model.dat', 'r')
mpop = pickle.load(f)
f.close()

f = open('pop_vdisp_model.dat', 'r')
#f = open('pop_vdisp_constsigma_model.dat', 'r')
vpop = pickle.load(f)
f.close()

f = open('pop_mstar_constsigma_model.dat', 'r')
mpop_constsigma = pickle.load(f)
f.close()

f = open('pop_vdisp_constsigma_model.dat', 'r')
vpop_constsigma = pickle.load(f)
f.close()

f = open('gargiulo_table.txt', 'r')
gar_table = np.loadtxt(f)
f.close()

gar_z = gar_table[:4, 0]
gar_re = gar_table[:4, 1]
gar_errre = gar_table[:4, 2]
gar_n = gar_table[:4, 3]
gar_sigma = gar_table[:4, 4]
gar_errsigma = gar_table[:4, 5]
gar_mstar = gar_table[:4, 6] + 0.25

ngar = len(gar_z)

f = open('Belli2014a_modified.txt', 'r')
abelli_table = np.loadtxt(f, usecols=(3, 4, 5, 6, 7, 9, 10))
f.close()

abelli_z = abelli_table[:, 0]
abelli_re = abelli_table[:, 3]
abelli_n = abelli_table[:, 4]
abelli_sigma = abelli_table[:, 1]
abelli_errsigma = abelli_table[:, 2]
abelli_mstar = abelli_table[:, 5] + 0.25
abelli_mstar_err = abelli_table[:, 6]

nabelli = len(abelli_z)

f = open('vandeSande2013.txt', 'r')
vande_table = np.loadtxt(f, usecols=(2, 3, 4, 5, 6, 9, 10, 11, 14))
f.close()

f = open('vandeSande2013.txt', 'r')
vande_samples = np.loadtxt(f, usecols=(0, ), dtype=int)
f.close()

vande_labels = ['van de Sande et al. 2013', \
		'Bezanson et al. 2013', \
		'van Dokkum et al. 2009', \
		'Onodera et al. 2012', \
		'Cappellari et al. 2009', \
		'Newman et al. 2010', \
		'van der Wel et al. 2008', \
		'Toft et al. 2012']

vande_markers = ['o', 'v', 's', '*', '<', 'D', 'p', '>']

vande_colors = [(1., 0., 0.), \
		(0., 3/5., 0.), \
		(1./10., 1./5., 0.), \
		(4/5., 0., 2/5.), \
		(0., 0., 4./5.), \
		(1/5., 3./5., 1.), \
		(3/5., 0., 3./5.), \
		(0., 1., 1.)]

nlabels = len(vande_labels)

vande_z = vande_table[:, 0]
vande_re = vande_table[:, 1][vande_z < 2.]
vande_re_err = vande_table[:, 2][vande_z < 2.]
vande_n = vande_table[:, 3][vande_z < 2.]
vande_n_err = vande_table[:, 4][vande_z < 2.]
vande_sigma = vande_table[:, 5][vande_z < 2.]
vande_mstar = vande_table[:, 8][vande_z < 2.]
vande_errsigma = 0.5*(vande_table[:, 6] + vande_table[:, 7])

vande_samples = vande_samples[vande_z < 2.]

vande_z = vande_z[vande_z < 2.]

nvande = len(vande_z)

nsamp = 10000

f = open('shetty_table.txt', 'r')
she_table = np.loadtxt(f)
f.close()

she_z = she_table[:, 1]
she_sigma = she_table[:, 4]
she_errsigma = she_table[:, 5]
she_jam = she_table[:, 6]
she_salp = she_table[:, 7]
she_aimf = np.log10(she_jam/she_salp)

nshe = len(she_z)

she_aimf_err = 0.*she_aimf + 0.1
she_lsigma = np.log10(she_sigma)
she_lsigma_err = 0.*she_lsigma + 0.05

def fit_shetty_fixed_slope(sigma, aimf, guess=(0.)):

    spin_slope = spin_fit['a'].mean()

    apar = pymc.Normal('a', mu=spin_fit['a'].mean(), tau=1./(spin_fit['a'].std())**2, value=spin_slope)

    bpar = pymc.Uniform('b', lower=-1., upper=1., value=0.)
    scat = pymc.Uniform('s', lower=0., upper=1., value=0.1)

    pars = [bpar, scat]

    @pymc.deterministic
    def like(a=apar, b=bpar, s=scat):

	model_aimf = b + a*(sigma - 2.4)
	logp = (-0.5*(aimf - model_aimf)**2/s**2 - np.log(s)).sum()

	return logp

    @pymc.stochastic
    def logp(value=0., observed=True, p=pars):
	return like

    M = pymc.MCMC(pars)
    M.use_step_method(pymc.AdaptiveMetropolis, pars)
    M.sample(11000, 1000)

    return M.trace('b')[:].mean(), M.trace('b')[:].std()

she_fit_fixedslope = fit_shetty_fixed_slope(she_lsigma, she_aimf)
print she_fit_fixedslope
she_fit = fit_aimf_vs_sigma(she_aimf, she_aimf_err, she_lsigma, she_lsigma_err)


f = open('/gdrive/projects/SL2S_hierarch/evol_nfw_alpha_msps_sigma.dat', 'r')
chain = pickle.load(f)
f.close()

samp = np.random.choice(np.arange(90000), nsamp)

burnin = 10000
mdep_samp = chain['amstar'][burnin:][samp]
sdep_samp = chain['asigma'][burnin:][samp]
zdep_samp = chain['zalpha'][burnin:][samp]
aimf_samp = chain['calpha'][burnin:][samp]

nz = 200
nson = 81

z = np.linspace(0., 2., nz)
zson = np.linspace(0., 0.8, nson)

mstar_coeffs = np.zeros((nz, 3))
mstar_pivs = np.zeros((nz, 2))
vdisp_coeffs = np.zeros((nz, 3))
vdisp_pivs = np.zeros((nz, 2))
mstar_constsigma_coeffs = np.zeros((nz, 3))
vdisp_constsigma_coeffs = np.zeros((nz, 3))

band84 = np.zeros(nson)
band16 = np.zeros(nson)
band95 = np.zeros(nson)
band05 = np.zeros(nson)
band99 = np.zeros(nson)
band01 = np.zeros(nson)

for i in range(nz):
    mpars, mpivs = dm.fit_imf_coeff_around_average(np.log10(mpop.mstar_salp[:, i]), np.log10(mpop.veldisp[:, i]), \
                                                   np.log10(mpop.aimf[:, i]))
    mstar_coeffs[i, :] = mpars
    mstar_pivs[i] = mpivs

    vpars, vpivs = dm.fit_imf_coeff_around_average(np.log10(vpop.mstar_salp[:, i]), np.log10(vpop.veldisp[:, i]), \
                                                   np.log10(vpop.aimf[:, i]))
    vdisp_coeffs[i, :] = vpars
    vdisp_pivs[i] = vpivs

    mpars, mpivs = dm.fit_imf_coeff_around_average(np.log10(mpop_constsigma.mstar_salp[:, i]), \
                                                   np.log10(mpop_constsigma.veldisp[:, i]), \
                                                   np.log10(mpop_constsigma.aimf[:, i]))
    mstar_constsigma_coeffs[i, :] = mpars

    vpars, vpivs = dm.fit_imf_coeff_around_average(np.log10(vpop_constsigma.mstar_salp[:, i]), \
                                                   np.log10(vpop_constsigma.veldisp[:, i]), \
                                                   np.log10(vpop_constsigma.aimf[:, i]))
    vdisp_constsigma_coeffs[i, :] = vpars

for i in range(nson):
    maimf_samp = aimf_samp + zdep_samp*(z[i] - 0.3)
    band84[i] = np.percentile(maimf_samp, 84)
    band16[i] = np.percentile(maimf_samp, 16)
    band95[i] = np.percentile(maimf_samp, 95)
    band05[i] = np.percentile(maimf_samp, 5)
    band99[i] = np.percentile(maimf_samp, 99.7)
    band01[i] = np.percentile(maimf_samp, 0.3)
 
fig = pylab.figure()
ax = fig.add_subplot(111)
pylab.subplots_adjust(left=0.12, right=0.95, bottom=0.10, top=0.75)

pylab.fill_between(zson, band01, band99, color=rgb_alpha(bandcolor, 0.2))
pylab.fill_between(zson, band05, band95, color=rgb_alpha(bandcolor, 0.5))
pylab.fill_between(zson, band16, band84, color=rgb_alpha(bandcolor, 1.))

pylab.plot(z, mpop.imf_coeff[:, 0], label='$M_*$ model')
pylab.plot(z, vpop.imf_coeff[:, 0], label='$\sigma$ model')

pylab.plot(zson, band16 - 700., color=rgb_alpha(bandcolor, 1.), linewidth=7, label='S15')
#pylab.plot(z, mpop_constsigma.imf_coeff[:, 0], label='$M_*$ model (const. $\sigma$)', color='b', linestyle='--')
#pylab.plot(z, vpop_constsigma.imf_coeff[:, 0], label='$\sigma$ model (const. $\sigma)', color='g', linestyle='--')

spin_aimf_sample = spin_fit['b'] + spin_fit['a']*(2.4 - spin_fit['sigma_mu'])

spin_aimf = spin_aimf_sample.mean()
spin_aimf_err = spin_aimf_sample.std()

conroy_aimf_sample = conroy_fit['b'] + conroy_fit['a']*(2.4 - conroy_fit['sigma_mu'])

conroy_aimf = conroy_aimf_sample.mean()
conroy_aimf_err = conroy_aimf_sample.std()

she_aimf_sample = she_fit['b'] + she_fit['a']*(2.4 - she_fit['sigma_mu'])

she_aimf_mean = she_aimf_sample.mean()
she_aimf_err = she_aimf_sample.std()



def get_upper_limit(z, msps, msps_err, sigma, sigma_err, re, n, re_err=None, percentile=84.):

    nobj = len(z)
    upper_limits = []
    medians = []
    nabove = 0
    nbelow = 0
    for i in range(nobj):

        if z[i] < 2.:

            zind = int(100*z[i])

            if re_err is not None:
                re_samp = np.random.normal(re[i], re_err[i], nsamp)
            else:
                re_samp = re[i]*np.ones(nsamp)

            mstar_samp = np.random.normal(msps[i], msps_err[i], nsamp)
            sigma_samp = np.random.normal(sigma[i], sigma_err[i], nsamp)

            vgshift_samp = mdep_samp*(mstar_samp - 11.5) + \
                           sdep_samp*(np.log10(sigma_samp) - 2.4)
            #vgshift_samp = spin_fit['a']*(np.log10(sigma_samp) - 2.4)

            vgamma_samp = 2.*np.log10(sigma_samp) + 10. - np.log10(cgs.G) + np.log10(re_samp) - mstar_samp + \
                          - vgshift_samp - np.log10(cgs.M_Sun) + np.log10(cgs.kpc) + np.log10(splev(n[i], kvir_spline))

            vgamma_good = vgamma_samp[vgamma_samp==vgamma_samp]

            upper_limit = np.percentile(vgamma_good, percentile)

	    if upper_limit - vgamma_good.std() > vpop.imf_coeff[zind, 0]:
		nabove += 1
	    else:
		nbelow += 1

            upper_limits.append(upper_limit)
            medians.append(np.median(vgamma_good))

    return np.array(upper_limits), np.array(medians), nabove, nbelow

gar_upper_limits, gar_med, gar_above, gar_below = get_upper_limit(gar_z, gar_mstar, np.ones(ngar)*0.15, gar_sigma, gar_errsigma, \
                                                  gar_re, gar_n, re_err=gar_errre)

abelli_upper_limits, abelli_med, abelli_above, abelli_below = get_upper_limit(abelli_z, abelli_mstar, abelli_mstar_err, abelli_sigma, abelli_errsigma, \
                                      abelli_re, abelli_n)


vande_upper_limits, vande_med, vande_above, vande_below = get_upper_limit(vande_z, vande_mstar, 0.15*np.ones(nvande), vande_sigma, \
                                                vande_errsigma, vande_re, vande_n, re_err=vande_re_err)

print 'objects above:', gar_above + abelli_above + vande_above
print 'objects below:', gar_below + abelli_below + vande_below

corr_she_aimf = she_aimf.copy()
for i in range(nshe):
    zind = 100*int(she_z[i])
    she_shift = 0.90*(np.log10(she_sigma[i]) - 2.4)
    corr_she_aimf[i] = she_shift

she_aimf -= corr_she_aimf

med_she_aimf = np.median(she_aimf)
scat_she_aimf = she_aimf.std()
med_she_z = np.median(she_z)

pylab.errorbar(spin_z, spin_aimf, yerr=spin_aimf_err, xerr=0.025, fmt='+', marker='^', color='r', label='Spiniello et al. 2014')
pylab.errorbar(conroy_z, conroy_aimf, yerr=conroy_aimf_err, fmt='o', color=(3/5., 3/5., 0.), label='Conroy \& van Dokkum 2012')
#pylab.errorbar(med_she_z, she_aimf_mean, yerr=she_aimf_err, xerr=she_z.std(), fmt='o', color=(0, 1., 0.), label='Shetty & Cappellari 2014')
pylab.errorbar(med_she_z, she_fit_fixedslope[0], yerr=she_fit_fixedslope[1], xerr=0.5*(she_z.max() - she_z.min()), fmt='o', color='k', label='Shetty & Cappellari 2014')
#pylab.errorbar(med_she_z, med_she_aimf, yerr=scat_she_aimf, xerr=she_z.std(), fmt='o', color='k', \
#               label='Shetty & Cappellari 2014')

pylab.errorbar(gar_z, gar_upper_limits, yerr=gar_upper_limits-gar_med, uplims=True, fmt='+', color=(1., 0.2, 1.), marker='h', label='Gargiulo et al. (2014)')
pylab.errorbar(abelli_z[abelli_z<2.], abelli_upper_limits, yerr=abelli_upper_limits-abelli_med, uplims=True, fmt='+', color=(0., 2./5., 1./5.), \
               marker='x', label='Belli et al. (2014)')

for i in range(nlabels):
    this_sample = vande_samples == i
    if len(vande_z[this_sample]) > 0:
	pylab.errorbar(vande_z[this_sample], vande_upper_limits[this_sample], yerr=vande_upper_limits[this_sample]-vande_med[this_sample], fmt='+', uplims=True, color=vande_colors[i], \
               markeredgecolor='none', marker=vande_markers[i], label=vande_labels[i])

pylab.xlabel('$z$', fontsize=fsize)
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}(\log{M_*}=11.5,\, \log{\sigma}=2.4)$', fontsize=fsize)
pylab.ylim(-0.7, 0.8)

box = ax.get_position()
#ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height*0.9])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.35), numpoints = 1, ncol=3, fontsize=10)
pylab.savefig('timeevol.eps')
pylab.show()

