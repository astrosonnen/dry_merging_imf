import numpy as np
from cgsconstants import *
from scipy.interpolate import splrep,splev
import shmrs
import recipes
from scipy.integrate import quad

A = 0.0104
xitilde = 9.72e-3
alpha = 0.133
beta = -1.995
gamma = 0.263
etap = 0.0993

betaR = 0.6

apar = 1.11
bpar = 1.10
Mdot = 46.1/yr

Ngrid = 1001
lmstar_grid = np.linspace(9., 13., Ngrid)


def izfunc(z, z_0):
    return apar*(z - z_0) - (apar-1.)*np.log((1.+z)/(1.+z_0))


class ETG:

    def __init__(
                self, z_0=0., mstar_chab_0=1e11, mhalo_0=1e13, re_0=5., sigma_0=200., z_form=None):

        self.z_0 = z_0
        self.mstar_chab_0 = mstar_chab_0
        self.mhalo_0 = mhalo_0
        self.re_0 = re_0
        self.sigma_0 = sigma_0
        self.z_form = z_form
        self.dt_form = None
        self.imf_form = None
        self.z = None
        self.mhalo = None
        self.mstar_chab = None
        self.mstar_true = None
        self.re = None
        self.dmstar_chab_dz = None
        self.dmstar_true_dz = None
        self.mstardlnre_dz = None

    def get_sf_history(self, sf_recipe='veldisp'):

        if sf_recipe == 'veldisp':
            self.z_form = recipes.z_form_vdisp_func(np.log10(self.sigma_0))
            self.dt_form = recipes.dt_form_vdisp_func(np.log10(self.sigma_0))

        elif sf_recipe == 'mstar':
            self.z_form = recipes.z_form_mstar_func(np.log10(self.mstar_chab_0))
            self.z_form = recipes.z_form_mstar_func(np.log10(self.mstar_chab_0))

        else:
            raise ValueError("sf_recipe must be one between 'sigma' and 'mstar'.")


    def evolve(self, ximin=0.03, dz=0.01, z_up=2., imf_recipe='SigmaSF', imf_coeff=(0.1, 0.3)):

        self.z = np.arange(0., z_up, dz)

        nz = len(self.z)

        i_0 = int(self.z_0/(self.z[-1]-self.z[0])*(nz-1))
        i_form = int(self.z_form/(self.z[-1]-self.z[0])*(nz-1))

        self.mhalo = np.empty(nz, dtype='float')
        self.mhalo = 1e12*(
                        (self.mhalo_0/1e12)**(1.-bpar) - (1-bpar)/H0*Mdot/1e12*izfunc(self.z, self.z_0))**(1./(1.-bpar))

        self.dmstar_chab_dz = 0.*self.z
        self.dmstar_true_dz = 0.*self.z
        self.mstardlnre_dz = 0.*self.z

        self.mstar_chab = 0.*self.z + self.mstar_chab_0
        self.mstar_true = 0.*self.z
        self.re = 0.*self.z + self.re_0

        for i in range(0, nz):

            lmhalo_grid = shmrs.mhfunc(lmstar_grid, self.z[i])
            lmhalo_spline = splrep(lmhalo_grid, lmstar_grid)

            def rfunc(mhalo):
                return 10.**(splev(np.log10(mhalo), lmhalo_spline) - np.log10(mhalo))

            # calculates what is the minimum mass of galaxies that have already formed stars at this redshift.
            lmstar_form = recipes.lmstar_zfunc(self.z[i])
            ximh_min = shmrs.mhfunc(lmstar_form, self.z[i])
            ximinmin = 10.**ximh_min/self.mhalo[i]

            if ximinmin > ximin:
                ximin_eff = ximinmin
            else:
                ximin_eff = ximin

            if ximin_eff < 1:
                imz = quad(
                    lambda xi: rfunc(xi*self.mhalo[i])*xi**(beta+1.)*np.exp((xi/xitilde)**gamma), ximin_eff, 1.)[0]

                self.dmstar_chab_dz[i] = -A*imz*self.mhalo[i]*(self.mhalo[i]/1e12)**alpha*(1.+self.z[i])**etap

                imtz = quad(
                    lambda xi: recipes.satellite_imf(
                        np.log10(xi*self.mhalo[i]*rfunc(xi*self.mhalo[i])), imf_recipe, imf_coeff)* \
                               rfunc(xi*self.mhalo[i])*xi**(beta+1.)*np.exp((xi/xitilde)**gamma), ximin_eff, 1.)[0]

                self.dmstar_true_dz[i] = -A*imtz*self.mhalo[i]*(self.mhalo[i]/1e12)**alpha*(1.+self.z[i])**etap

                ire = quad(lambda xi: (2.-np.log(1.+xi**(2.-betaR))/np.log(1.+xi))*rfunc(xi*self.mhalo[i])* \
                               xi**(beta+1.)*np.exp((xi/xitilde)**gamma), ximin_eff, 1.)[0]

                self.mstardlnre_dz[i] = -A*ire*self.mhalo[i]*(self.mhalo[i]/1e12)**alpha*(1.+self.z[i])**etap

            else:
                pass
                # no mergers at this z. Galaxies below this mass have not formed stars.
                # but then our central galaxy shouldn't even exist!!

        for i in range(i_0-1, -1, -1):
            self.mstar_chab[i] = self.mstar_chab[i+1] - 0.5*(self.dmstar_chab_dz[i] + self.dmstar_chab_dz[i+1])*dz
            self.re[i] = self.re[i+1]*(1. - 0.5*(self.mstardlnre_dz[i]/self.mstar_chab[i] + \
                                                 self.mstardlnre_dz[i+1]/self.mstar_chab[i+1])*dz)

        for i in range(i_0+1, nz):
            self.mstar_chab[i] = self.mstar_chab[i-1] + 0.5*(self.dmstar_chab_dz[i] + self.dmstar_chab_dz[i-1])*dz
            self.re[i] = self.re[i-1]*(1. + 0.5*(self.mstardlnre_dz[i]/self.mstar_chab[i] + \
                                                 self.mstardlnre_dz[i-1]/self.mstar_chab[i-1])*dz)

        if imf_recipe == 'SigmaSF':
            self.imf_form = 10.**recipes.limf_func_cvd12(self.mstar_chab[i_form], self.re[i_form], self.dt_form, imf_coeff)
        elif imf_recipe == 'density':
            self.imf_form = 10.**recipes.limf_func_rhoc(self.z_form, imf_coeff)
        else:
            raise ValueError("recipe must be one between 'SigmaSF' and 'density'.")

        self.mstar_true = self.mstar_chab[i_form]*self.imf_form + 0.*self.z

        for i in range(i_form-1, -1, -1):
            self.mstar_true[i] = self.mstar_true[i+1] + 0.5*(self.dmstar_true_dz[i] + self.dmstar_true_dz[i+1])*dz
