import numpy as np
from cgsconstants import *
from cosmolopy import density, distance
from scipy.interpolate import splrep, splev

Nz = 1001
z_grid = np.linspace(0., 5., Nz)
logt_grid = 0.*z_grid
rhoc_grid = 0.*z_grid

for i in range(0, Nz):
    logt_grid[i] = np.log10(distance.lookback_time(z_grid[i], 0., **cosmo)/yr) - 9.
    rhoc_grid[i] = density.cosmo_densities(**cosmo)[0]*distance.e_z(z_grid[i], **cosmo)*M_Sun/Mpc**3

z_spline = splrep(logt_grid[1:], z_grid[1:])
logt_spline = splrep(z_grid[1:], logt_grid[1:])
rhoc_spline = splrep(z_grid[1:], rhoc_grid[1:])


def z_form_vdisp_func(lvdisp):
    return splev(0.46 + 0.238*lvdisp, z_spline)


def z_form_mstar_func(lmstar):
    return splev(0.427 + 0.053*lmstar, z_spline)


def dt_form_vdisp_func(lvdisp):
    return 10.**(3.44 - 1.68*lvdisp + 9.)


def dt_form_mstar_func(lmstar):
    return 10.**(3.67 - 0.37*lmstar + 9.)


def lmstar_zfunc(z):
    logt = splev(z, logt_spline)
    return logt/0.053 - 0.427/0.053


def sigmasf(mstar, re, dt):
    return 0.5*mstar/(np.pi*re**2)/dt


def vdisp_mstar_rel_thomas(lmstar):  # from Thomas et al. 2005 eq. (2)
    return 10.**((lmstar - 0.63)/4.52)


def vdisp_mstar_rel_auger(lmstar):  # from Auger et al. 2010
    return 10.**(2.34 + 0.18*(lmstar - 11.))


def vdisp_mstar_rel_mason(lmstar, z):  # from Auger et al. 2010
    return 10.**(2.34 + 0.18*(lmstar - 11.) + 0.20*np.log10(1. + z))


def inverse_vdisp_mstar_rel_auger(lvdisp):
    return 10.**((lvdisp - 2.34)/0.18 + 11.)


def re_mstar_rel_z0(lmstar):  # from Newman et al. 2012, SDSS bin
    return 10.**(0.54 + 0.57*(lmstar - 11.))


def re_mstar_rel(lmstar, z):  # from Newman et al. 2012, z-dependent
    return 10.**(0.38 + 0.57*(lmstar - 11.) - 0.26*(z - 1.))


def generate_reff(lmstar_sample, z):  # draws values of Re from the mass-radius relation of Newman et al. (2012)
    return re_mstar_rel(lmstar_sample, z) + np.random.normal(0., 0.22, len(np.atleast_1d(lmstar_sample)))


def generate_veldisp_from_fp(lmstar_sample, reff_sample):
    """
    draws velocity dispersions given mstar and Re using the stellar mass fundamental plane by Hyde & Bernardi (2009b).
    does NOT add any scatter! This is because this function is meant to be used when scatter has already been added
    in generating the values of Re...
    Is this consistent with function vdisp_mstar_rel from Thomas et al. 2005? Don't know...
    :param lmstar_sample:
    :param reff_sample:
    :return:
    """

    a = 1.3989
    b = 0.3164
    c = 4.4858
    scat = 0.0894

    return 10.**(1./a*(np.log10(reff_sample) + 2.5*b*(lmstar_sample - 2.*np.log10(reff_sample) - np.log10(2.*np.pi)) - c))


def generate_veldisp_from_mstar(lmstar_sample, z):
    return 10.**(np.log10(vdisp_mstar_rel_mason(lmstar_sample, z)) + np.random.normal(0., 0.04, len(np.atleast_1d(lmstar_sample))))


def limf_func_cvd12(mstar, re, dt, coeff=(0.1, 0.3)):
    return coeff[0] + coeff[1]*(np.log10(sigmasf(mstar, re, dt)) - 1.)


def limf_func_rhoc(z_form, coeff=(0.3, 1.0)):
    rhoc = splev(z_form, rhoc_spline)
    return coeff[0] + coeff[1]*(np.log10(rhoc) + 28.)


def limf_func_mstar(lmstar, coeff=(3., 0.2)):
    #return (2./np.pi*np.arctan((lmstar - 11.2)*coeff[0]) + 0.8)*coeff[1]
    return coeff[0]*(lmstar - 11.) + coeff[1]


def limf_func_vdisp(lvdisp, coeff=(3., 0.2)):
    #return (2./np.pi*np.arctan((lvdisp - 2.3)*coeff[0]) + 0.8)*coeff[1]
    return coeff[0]*(lvdisp - 2.3) + coeff[1]


def limf_func_mhalo(lmhalo, coeff=(0.3, 0.0)):
    #return (2./np.pi*np.arctan((lmstar - 11.2)*coeff[0]) + 0.8)*coeff[1]
    return coeff[0]*(lmhalo - 12.) + coeff[1]


def satellite_imf(lmstar, z=2., recipe='SigmaSF', coeff=(0.1, 0.3), lmhalo=None):

    if recipe == 'SigmaSF':
        dt_form = dt_form_mstar_func(lmstar)
        re = re_mstar_rel_z0(lmstar)
        return 10.**limf_func_cvd12(10.**lmstar, re, dt_form, coeff)

    elif recipe == 'density':
        z_form = z_form_mstar_func(lmstar)
        return 10.**limf_func_rhoc(z_form, coeff)

    elif recipe == 'mstar' or recipe == 'mstar-wscatter':
        return 10.**(limf_func_mstar(lmstar, coeff))

    elif recipe == 'vdisp':
        return 10.**(limf_func_vdisp(np.log10(vdisp_mstar_rel_mason(lmstar, z)), coeff))

    elif recipe == 'mstar-vdisp':
        return 10.**(limf_func_mstar(lmstar, (coeff[0], coeff[1])) + \
                     limf_func_vdisp(np.log10(vdisp_mstar_rel_mason(lmstar, z)), (coeff[2], coeff[3])))

    elif recipe == 'mhalo':
        return 10.**(limf_func_mhalo(lmhalo, coeff))

    else:
        raise ValueError("recipe must be one between 'SigmaSF' and 'density'.")

