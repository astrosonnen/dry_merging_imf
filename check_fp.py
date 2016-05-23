import pylab
import pickle
import numpy as np
import recipes
from scipy.optimize import leastsq

f = open('mstar_dep_imf_coeff0.33.dat', 'r')
#f = open('vdisp_dep_imf_coeff2.30.dat', 'r')
galaxies = pickle.load(f)
f.close()

afp = 1.185
bfp = -0.952
gfp = 0.261

snap = 0

ngal = len(galaxies)

reffs = []
mstars = []
vdisps = []
aimfs = []
xieff = np.zeros(ngal)
reff2 = []
mstar2 = []
vdisp2 = []
aimf2 = []
for i in range(ngal):
    reffs.append(galaxies[i].re[snap])
    mstars.append(galaxies[i].mstar_chab[snap])
    vdisps.append(galaxies[i].veldisp[snap])
    aimfs.append(galaxies[i].aimf[snap])

    dmstars = galaxies[i].mstar_chab[:-1] - galaxies[i].mstar_chab[1:]
    #zeff[i] = sum(galaxies[i].z[:-1]*dmstars)/sum(dmstars)
    xieff[i] = sum(galaxies[i].xieff[:-1]*dmstars)/sum(dmstars)

    reff2.append(galaxies[i].re[199])
    mstar2.append(galaxies[i].mstar_chab[199])
    vdisp2.append(galaxies[i].veldisp[199])
    aimf2.append(galaxies[i].aimf[199])


reffs = np.log10(np.array(reffs))
mstars = np.log10(np.array(mstars))
vdisps = np.log10(np.array(vdisps))
aimfs = np.log10(np.array(aimfs))
z0vdisp = np.log10(recipes.vdisp_mstar_rel_auger(mstars))
z0reff = np.log10(recipes.re_mstar_rel_auger(mstars))

reff2 = np.log10(np.array(reff2))
mstar2 = np.log10(np.array(mstar2))
vdisp2 = np.log10(np.array(vdisp2))
aimf2 = np.log10(np.array(aimf2))

reff_if_newman = np.log10(recipes.re_mstar_rel_z0(mstars))

def fit_mstar_re(m, r, guess=(0.6, 0.6)):
    def modelfunc(p):
        return p[0] + p[1]*(m - 11.)

    def errfunc(p):
        return modelfunc(p) - r

    par, cov = leastsq(errfunc, guess)

    return par


def fit_aimf_fp(m, r, v, a, guess=(0., 0.6, 0.6)):

    def v1(m, r, v):
        return afp*(v - 2.) + bfp*(m - 9.) - bfp*np.log10(2.*np.pi) - bfp*r + gfp

    def modelfunc(p):
        return p[0] + p[1]*v1(m, r, v) + p[2]*r

    def errfunc(p):
        return modelfunc(p) - a

    par, cov = leastsq(errfunc, guess)

    return par

"""
pylab.scatter(reffs, z0reff)
xlim = pylab.xlim()
xs = np.linspace(xlim[0], xlim[1])
pylab.plot(xs, xs, linestyle='--', color='k')
pylab.show()
"""


pylab.scatter(mstars, reffs)
fit = fit_mstar_re(mstars, reffs)
print fit
xlim = pylab.xlim()
xs = np.linspace(xlim[0], xlim[1])
pylab.plot(xs, fit[0] + fit[1]*(xs - 11.))
pylab.plot(xs, np.log10(recipes.re_mstar_rel_z0(xs)), linestyle='--', color='r')
pylab.show()

"""
pylab.scatter(vdisps+2., z0vdisp+2.)
xlim = pylab.xlim()
xs = np.linspace(xlim[0], xlim[1])
pylab.plot(xs, xs, linestyle='--', color='k')
pylab.show()
"""

fsize=14

def plot_fp(mstar, re, vdisp, zdim, zlabel=None):
    pylab.scatter(afp*(vdisp - 2.) + bfp*(mstar - 9.) - bfp*np.log10((2.*np.pi)) - 2.*bfp*re + gfp, re, \
                  c=zdim, cmap=pylab.cm.coolwarm)
    pylab.colorbar(label=zlabel)

    imffit = fit_aimf_fp(mstar, re, vdisp, zdim)
    print imffit

    xlim = pylab.xlim(-1, 2)
    pylab.ylim(xlim[0], xlim[1])
    xs = np.linspace(xlim[0], xlim[1])
    pylab.plot(xs, xs, linestyle='--', color='k')

pylab.subplots_adjust(wspace=0.3, hspace=0.3)

pylab.subplot(2, 2, 1)
plot_fp(mstars, reffs, vdisps, aimfs, zlabel='$\log{\\alpha_{\mathrm{IMF}}}$')
pylab.ylabel('$\log{r_e}$', fontsize=fsize)

pylab.subplot(2, 2, 2)
plot_fp(mstars, reffs, vdisps, mstars, zlabel='$\log{M_*}$')
pylab.xlabel('$1.18\log{\sigma} - 0.95\log{M_*/2\pi r_e^2} + 0.26$', fontsize=fsize)

pylab.subplot(2, 2, 3)
plot_fp(mstars, reffs, vdisps, vdisps, zlabel='$\log{\sigma}$')
pylab.xlabel('$1.18\log{\sigma} - 0.95\log{M_*/2\pi r_e^2} + 0.26$', fontsize=fsize)
pylab.ylabel('$\log{r_e}$', fontsize=fsize)

"""
pylab.subplot(2, 2, 4)
plot_fp(mstars, reff_if_newman, vdisps + 0.15, aimfs, zlabel='$\log{\\alpha_{\mathrm{IMF}}}(z=2)$')
"""

pylab.savefig('fp.png')
pylab.show()

