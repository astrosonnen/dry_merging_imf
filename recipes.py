from sonnentools.cgsconstants import *
import numpy as np
import cosmology
from scipy.interpolate import splrep,splev

Nz = 1001
z_grid = np.linspace(0.,5.,Nz)
logt_grid = 0.*z_grid
rhoc_grid = 0.*z_grid

for i in range(1,Nz):
    logt_grid[i] = np.log10(cosmology.lookback(z_grid[i])) - 9.
    rhoc_grid[i] = cosmology.rhoc(z_grid[i])

z_spline = splrep(logt_grid[1:],z_grid[1:])
logt_spline = splrep(z_grid[1:],logt_grid[1:])
rhoc_spline = splrep(z_grid[1:],rhoc_grid[1:])

def z_form_mstar_func(lmstar):
    return splev(0.427 + 0.053*lmstar,z_spline)
    #return 1.9856 + 0.8713*(lmstar -11.) + 0.3451*(lmstar - 11.)**2

def z_form_sigma_func(lsigma):
    return splev(0.46 + 0.238*lsigma,z_spline)

def dtform_func(lsigma):
    return 10.**(3.44 - 1.68*lsigma + 9.)

def lmstar_zfunc(z):
    logt = splev(z,logt_spline)
    return logt/0.053 - 0.427/0.053

def SigmaSF(mstar,re,dt):
    return 0.5*mstar/(np.pi*re**2)/dt


def sigma_from_mstar(lmstar): #from Thomas et al. 2005 eq. (2)
    return 10.**((lmstar - 0.63)/4.52)


def re_from_mstar(lmstar): #from Newman et al. 2012, SDSS bin
    return 10.**(0.54 + 0.57*(lmstar - 11.))


def limf_func_CvD12(mstar,re,dt,a=0.1,b=0.3):
    return a + b*(np.log10(SigmaSF(mstar,re,dt)) - 1.)
    

def limf_func_rhoc(z_form,a=0.3,b=1.):
    #rhoc = cosmology.rhoc(z_form_func(lmstar))
    rhoc = splev(z_form,rhoc_spline)
    return a + b*(np.log10(rhoc) + 28.)



def central_imf(galaxy, recipe='sigma', coeff=(0.1,0.3)):
    """

    :param galaxy: ETG class object
    :param recipe: string. Allowed values are 'density' and 'SigmaSF'
    :param coeff: two-element tuple of floats
    :return: float. IMF-normalization coefficient relative to a Chabrier IMF

    This function assigns an IMF to the stellar population of a galaxy at the time of star formation.
    Two recipes for assigning the IMF are implemented.
    'density': We use
    """

    if recipe == 'SigmaSF':
        return 10.**limf_func_CvD12(galaxy.mstar_chab[-1], galaxy.re[-1], galaxy.dtform, coeff[0], coeff[1])

    elif recipe == 'density':
        return 10.**limf_func_rhoc(galaxy.z_form, coeff[0], coeff[1])

    else:
        raise ValueError("recipe must be one between 'SigmaSF' and 'density'.")


def satellite_imf(lmstar, recipe='SigmaSF', coeff=(0.1,0.3)):

    if recipe == 'SigmaSF':
        dtform = dtform_func(np.log10(sigma_from_mstar(lmstar)))
        re = re_from_mstar(lmstar)
        return 10.**limf_func_CvD12(10.**lmstar,re,dtform,coeff[0],coeff[1])

    elif recipe == 'mstar':
        z_form = z_form_mstar_func(lmstar)
        return 10.**limf_func_rhoc(z_form,coeff[0],coeff[1])

    else:
         raise ValueError("recipe must be one between 'sigma' and 'mstar'.")

