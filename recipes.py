from sonnentools.cgsconstants import *
import numpy as np
from cosmolopy import density, distance
from scipy.interpolate import splrep,splev

Nz = 1001
z_grid = np.linspace(0., 5., Nz)
logt_grid = 0.*z_grid
rhoc_grid = 0.*z_grid

for i in range(1, Nz):
    logt_grid[i] = np.log10(distance.lookback_time(z_grid[i]), **cosmo) - 9.
    rhoc_grid[i] = density.cosmo_densities(**cosmo)*distance.e_z(z_grid[i], **cosmo)

z_spline = splrep(logt_grid[1:], z_grid[1:])
logt_spline = splrep(z_grid[1:], logt_grid[1:])
rhoc_spline = splrep(z_grid[1:], rhoc_grid[1:])

def z_form_mstar_func(lmstar):
    return splev(0.427 + 0.053*lmstar, z_spline)


def z_form_vdisp_func(lvdisp):
    return splev(0.46 + 0.238*lvdisp,z_spline)


def dtform_func(lvdisp):
    return 10.**(3.44 - 1.68*lvdisp + 9.)


def lmstar_zfunc(z):
    logt = splev(z,logt_spline)
    return logt/0.053 - 0.427/0.053


def SigmaSF(mstar,re,dt):
    return 0.5*mstar/(np.pi*re**2)/dt


def vdisp_mstar_rel(lmstar): #from Thomas et al. 2005 eq. (2)
    return 10.**((lmstar - 0.63)/4.52)


def re_mstar_rel(lmstar): #from Newman et al. 2012, SDSS bin
    return 10.**(0.54 + 0.57*(lmstar - 11.))


def limf_func_CvD12(mstar, re, dt, coeff=(0.1,0.3)):
    return coeff[0] + coeff[1]*(np.log10(SigmaSF(mstar, re, dt)) - 1.)


def limf_func_rhoc(z_form, coeff=(0.3, 1.0)):
    rhoc = splev(z_form, rhoc_spline)
    return coeff[0] + coeff[1]*(np.log10(rhoc) + 28.)


def central_imf(galaxy, recipe='SigmaSF', coeff=(0.1,0.3)):
    """

    :param galaxy: ETG class object
    :param recipe: string. Allowed values are 'density' and 'SigmaSF'
    :param coeff: two-element tuple of floats
    :return: float. IMF-normalization coefficient relative to a Chabrier IMF

    This function assigns an IMF to the stellar population of a galaxy at the time of star formation.
    Two recipes for assigning the IMF are implemented.
    'SigmaSF': The IMF normalization scales with the surface density of star formation, following the idea of
    Conroy and van Dokkum (2012).
    'density': The IMF normalization scales with the critical density of the Universe at the time of star formation
    """

    if recipe == 'SigmaSF':
        return 10.**limf_func_CvD12(galaxy.mstar_chab[-1], galaxy.re[-1], galaxy.dtform, coeff[0], coeff[1])

    elif recipe == 'density':
        return 10.**limf_func_rhoc(galaxy.z_form, coeff)

    else:
        raise ValueError("recipe must be one between 'SigmaSF' and 'density'.")


def satellite_imf(lmstar, recipe='SigmaSF', coeff=(0.1,0.3)):

    if recipe == 'SigmaSF':
        dtform = dtform_func(np.log10(vdisp_mstar_rel(lmstar)))
        re = re_mstar_rel(lmstar)
        return 10.**limf_func_CvD12(10.**lmstar, re, dtform, coeff[0], coeff[1])

    elif recipe == 'mstar':
        z_form = z_form_mstar_func(lmstar)
        return 10.**limf_func_rhoc(z_form, coeff[0], coeff[1])

    else:
         raise ValueError("recipe must be one between 'SigmaSF' and 'mstar'.")

