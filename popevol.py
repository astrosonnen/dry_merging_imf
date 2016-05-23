import numpy as np
from cgsconstants import *
from scipy.interpolate import splrep,splev
import shmrs
import recipes
from scipy.integrate import quad
import do_measurements

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


class population:

    def __init__(
                self, z_0=2., nobj=100, mstar_chab_0=None, mhalo_0=None, veldisp_0=None):

        self.z_0 = z_0
        self.nobj = nobj
        self.mstar_chab_0 = mstar_chab_0
        self.mstar_salp_0 = mstar_chab_0*10.**0.25
        self.mhalo_0 = mhalo_0
        self.veldisp_0 = veldisp_0
        self.z = None
        self.mhalo = None
        self.mstar_chab = None
        self.mstar_true = None
        self.veldisp = None
        self.aimf = None
        self.xieff = None
        self.rfuncatxieff = None
        self.dlnsigma2_dz = None

    def evolve(self, ximin=0.03, dz=0.01, z_low=0., imf_recipe='mstar', imf_coeff=(0.1, 0.3), vdisp_coeff=(2.40, 0.18)):

        self.z = np.arange(z_low, self.z_0, dz)
        nz = len(self.z)

        self.mhalo = np.zeros((self.nobj, nz), dtype='float')
        self.mstar_chab = 0.*self.mhalo
        self.mstar_salp = 0.*self.mhalo
        self.mstar_true = 0.*self.mhalo
        self.veldisp = 0.*self.mhalo
        self.aimf = 0.*self.mhalo
        self.dlnsigma2_dz = 0.*self.mhalo

        self.mstar_chab[:, -1] = self.mstar_chab_0
        self.mstar_salp[:, -1] = self.mstar_salp_0
        self.veldisp[:, -1] = self.veldisp_0

        for i in range(self.nobj):
            self.mhalo[i, :] = 1e12*((self.mhalo_0[i]/1e12)**(1.-bpar) - (1-bpar)/H0*Mdot/1e12*izfunc(self.z, self.z_0))**(1./(1.-bpar))
            if imf_recipe == 'mstar':
                self.aimf[i, -1] = 10.**recipes.limf_func_mstar(np.log10(self.mstar_salp_0[i]), imf_coeff)
            elif imf_recipe == 'vdisp':
                self.aimf[i, -1] = 10.**recipes.limf_func_vdisp(np.log10(self.veldisp_0[i]), imf_coeff)
            elif imf_recipe == 'mhalo':
                self.aimf[i, -1] = 10.**recipes.limf_func_mhalo(np.log10(self.mhalo_0[i]), imf_coeff)

        self.mstar_true[:, -1] += self.mstar_salp_0 * self.aimf[:, -1]

        self.xieff = 0.*self.mhalo
        self.rfuncatxieff = 0.*self.mhalo

        for i in range(nz-2, -1, -1):

            lmhalo_grid = shmrs.mhfunc(lmstar_grid, self.z[i])
            lmhalo_spline = splrep(lmhalo_grid, lmstar_grid)

            def rfunc(mhalo):
                return 10.**(splev(np.log10(mhalo), lmhalo_spline) - np.log10(mhalo))

            for j in range(self.nobj):
                imz = quad(
                    lambda xi: rfunc(xi*self.mhalo[j, i])*xi**(beta+1.)*np.exp((xi/xitilde)**gamma), ximin, 1.)[0]
                imxiz = quad(
                    lambda xi: rfunc(xi*self.mhalo[j, i])*xi**(beta+2.)*np.exp((xi/xitilde)**gamma), ximin, 1.)[0]

                self.xieff[j, i] = imxiz/imz
                self.rfuncatxieff[j, i] = rfunc(self.xieff[j, i]*self.mhalo[j, i])

                dmstar_chab_dz = -A*imz*self.mhalo[j, i]*(self.mhalo[j, i]/1e12)**alpha*(1.+self.z[i])**etap

                def integrand(xi):
                    thing = recipes.satellite_imf(np.log10(xi*self.mhalo[j, i]*rfunc(xi*self.mhalo[j, i])) + 0.25, \
                                                     z=self.z[i], recipe=imf_recipe, coeff=imf_coeff, \
                                                     lmhalo=np.log10(xi*self.mhalo[j, i]))* \
                           rfunc(xi*self.mhalo[j, i])*10.**0.25*xi**(beta+1.)*np.exp((xi/xitilde)**gamma)
                    return thing

                imtz = quad(integrand, ximin, 1.)[0]

                dmstar_true_dz = -A*imtz*self.mhalo[j,i]*(self.mhalo[j, i]/1e12)**alpha*(1.+self.z[i])**etap

                def epsilon(xi):
                    mstar_sat = rfunc(self.mhalo[j, i]*xi)*xi*self.mhalo[j, i]
                    vdisp_sat = 10.**(vdisp_coeff[0] + vdisp_coeff[1]*(np.log10(mstar_sat) - 11.))
                    return (vdisp_sat/self.veldisp[j, i+1])**2

                isigma2 = quad(lambda xi: ((1.+xi*epsilon(xi))/(1.+xi) - 1)*xi**beta*np.exp((xi/xitilde)**gamma), \
                               ximin, 1.)[0]

                self.dlnsigma2_dz[j, i] = -A*isigma2*(self.mhalo[j, i]/1e12)**alpha*(1.+self.z[i])**etap

                self.mstar_chab[j, i] = self.mstar_chab[j, i+1] - dmstar_chab_dz*dz
                self.mstar_salp[j, i] = self.mstar_chab[j, i]*10.**0.25
                self.mstar_true[j, i] = self.mstar_true[j, i+1] - dmstar_true_dz*dz

                self.veldisp[j, i] = self.veldisp[j, i+1]*(1. - 0.5*self.dlnsigma2_dz[j, i]*dz)

            # now fits for the vdisp - mstar relation and imf - mstar relation
            fit_vdisp_coeff = do_measurements.fit_mstar_only(np.log10(self.mstar_chab[:, i]), np.log10(self.veldisp[:, i]), \
                                                 guess=vdisp_coeff)[0]

            vdisp_coeff = fit_vdisp_coeff
            self.aimf[:, i] = self.mstar_true[:, i] / self.mstar_salp[:, i]

            if imf_recipe == 'mstar':
                imf_coeff = do_measurements.fit_mstar_only(np.log10(self.mstar_salp[:, i]), np.log10(self.aimf[:, i]), \
                                                   guess=imf_coeff)[0]
            elif imf_recipe == 'vdisp':
                 imf_coeff = do_measurements.fit_sigma_only(np.log10(self.veldisp[:, i]), np.log10(self.aimf[:, i]), \
                                                   guess=imf_coeff)[0]
            elif imf_recipe == 'mhalo':
                  imf_coeff = do_measurements.fit_mhalo_only(np.log10(self.mhalo[:, i]), np.log10(self.aimf[:, i]), \
                                                   guess=imf_coeff)[0]

            print i, imf_coeff, fit_vdisp_coeff



