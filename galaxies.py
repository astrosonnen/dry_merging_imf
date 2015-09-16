import numpy as np
from sonnentools.cgsconstants import *
from scipy.interpolate import splrep,splev
import shmrs
import imf_func
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
lmstar_grid = np.linspace(9.,13.,Ngrid)

def Izfunc(z,z_0):
    return apar*(z - z_0) - (apar-1.)*np.log((1.+z)/(1.+z_0))



class etg:

    def __init__(self,mstar_chab_0=1e11,mhalo_0=1e13,re_0=5.,sigma_0=200.,zform=None,dtform=None):
        self.mstar_chab_0 = mstar_chab_0
        self.mhalo_0 = mhalo_0
        self.re_0 = re_0
        self.sigma_0 = sigma_0
        self.zform = zform
	self.dtform = None
        self.imf_form = None
        self.z = None
        self.mhalo = None
        self.mstar_chab = None
        self.mstar_true = None
	self.re = None

    def get_zform(self):
        #self.zform = imf_func.zform_mstar_func(np.log10(self.mstar_chab_0))
        self.zform = imf_func.zform_sigma_func(np.log10(self.sigma_0))

    def get_dtform(self):
	self.dtform = imf_func.dtform_func(np.log10(self.sigma_0))


    def make_z_grid(self,dz=0.001,zup=None):
        if zup is not None:
            self.z = np.arange(0.,zup,dz)
        else:
            if self.zform is None:
                self.get_zform()

            self.z = np.arange(0.,self.zform,dz)


    def evolve_back(self,ximin=0.03,zup=None,dz=0.001,usesigma=True,coeff=(0.1,0.3)):
        if self.z == None:
            self.make_z_grid(dz,zup)
        else:
            dz = self.z[1] - self.z[0]

        Nz = len(self.z)


        self.mhalo = 1e12*((self.mhalo_0/1e12)**(1.-bpar) - (1-bpar)/H0*Mdot/1e12*Izfunc(self.z,0.))**(1./(1.-bpar))

        self.mstar_chab = 0.*self.mhalo + self.mstar_chab_0
        self.mstar_true = 0.*self.mhalo + self.mstar_chab_0

	self.re = 0.*self.mhalo + self.re_0

        for i in range(1,Nz):
            lmhalo_grid = shmrs.mhfunc(lmstar_grid,self.z[i])
            lmhalo_spline = splrep(lmhalo_grid,lmstar_grid)

            rfunc = lambda mhalo: 10.**(splev(np.log10(mhalo),lmhalo_spline) - np.log10(mhalo))

            lmstar_form = imf_func.lmstar_zfunc(self.z[i])
            ximh_min = shmrs.mhfunc(lmstar_form,self.z[i])
            ximinmin = 10.**ximh_min/self.mhalo[i]
            if ximinmin > ximin:
                ximin_eff = ximinmin
            else:
                ximin_eff = ximin

            if ximin_eff < 1:
                IMz = quad(lambda xi: rfunc(xi*self.mhalo[i])*xi**(beta+1.)*np.exp((xi/xitilde)**gamma),ximin_eff,1.)[0]

                dmdz = -A*IMz*self.mhalo[i]*(self.mhalo[i]/1e12)**alpha*(1.+self.z[i])**etap

                IMtz = quad(lambda xi: imf_func.satellite_imf(np.log10(xi*self.mhalo[i]*rfunc(xi*self.mhalo[i])),usesigma,coeff)*rfunc(xi*self.mhalo[i])*xi**(beta+1.)*np.exp((xi/xitilde)**gamma),ximin_eff,1.)[0]

                dmtdz = -A*IMtz*self.mhalo[i]*(self.mhalo[i]/1e12)**alpha*(1.+self.z[i])**etap

		Ire = quad(lambda xi: (2.-np.log(1.+xi**(2.-betaR))/np.log(1.+xi))*rfunc(xi*self.mhalo[i])*xi**(beta+1.)*np.exp((xi/xitilde)**gamma),ximin_eff,1.)[0]

		dlnredz = -A*Ire*self.mhalo[i]/(self.mstar_chab[i-1] + dz*dmdz)*(self.mhalo[i]/1e12)**alpha*(1.+self.z[i])**etap

            else:
                dmdz = 0.
                dmtdz = 0.
		dlnredz = 0.

            self.mstar_chab[i] = self.mstar_chab[i-1] + dz*dmdz
            self.mstar_true[i] = self.mstar_true[i-1] + dz*dmtdz
	    self.re[i] = self.re[i-1]*(1. + dlnredz*dz)


        if zup is None:
            self.imf_form = imf_func.central_imf(self,usesigma,coeff)
            self.mstar_true += self.mstar_chab[-1]*self.imf_form - self.mstar_true[-1]



   
