import numpy as np
from scipy.interpolate import splrep,splev

N = 1001

lmstars = np.linspace(9.,13.,N)

def mhalo_dist(lmhalo,z):    #halo mass distribution from Tinker et al. 2008. Not normalized.
    lmref = 10. - np.log10(0.7)
    logsigminus = -0.64 + 0.11*(lmhalo - lmref) + 0.012*(lmhalo - lmref)**2 + 0.0006*(lmhalo - lmref)**3
    dlogsigminusdlogm = 0.11 + 2.*0.012*(lmhalo - lmref) + 3.*0.0006*(lmhalo - lmref)**2
    sigma = 10.**-logsigminus

    A = 0.186*(1.+z)**(-0.14)
    a = 1.47*(1.+z)**(-0.06)
    Delta = 200.
    alpha = 10.**(-(0.75/(np.log10(Delta/75.)))**1.2)
    b = 2.57*(1.+z)**alpha
    c = 1.19

    return ((sigma/b)**-a + 1)*np.exp(-c/sigma**2)/10.**lmhalo*dlogsigminusdlogm



def mhfunc(lmstar,z):#,mstar00,mstar0a,mstar0a2,m10,m1a,beta0,betaa,delta0,deltaa,gamma0):

    mstar00 = 10.78
    mstar0a = 0.36
    mstar0a2 = 0.
    m10 = 12.40
    m1a = 0.38
    beta0 = 0.45
    betaa = 0.026
    delta0 = 0.56
    deltaa = 0.
    gamma0 = 0.82
    gammaa = 1.86

    a = 1./(1. + z)
    logm1 = m10 + m1a*(a - 1.)
    logms0 = mstar00 + mstar0a*(a - 1.)
    beta = beta0 + betaa*(a-1.)
    delta = delta0 + deltaa*(a-1.)
    gamma = gamma0 + gammaa*(a-1.)

    return logm1 + beta*(lmstar-logms0) + (10.**(lmstar-logms0))**delta/(1.+(10.**(lmstar-logms0))**(-gamma)) - 0.5

def mstarfunc_z0(lmhalo):
    lmstars = np.linspace(9.,12.,1001)
    lmhalos = mhfunc(lmstars,0.)
    lmstar_spline = splrep(lmhalos,lmstars)
    return splev(lmhalo,lmstar_spline)
    

def rstarh(lmhalo,z):

    lmhalos = mhfunc(lmstars,z)#,mstar00,mstar0a,mstar0a2,m10,m1a,beta0,betaa,delta0,deltaa,gamma0)

    lmh_spline = splrep(lmhalos,lmstars)

    return 10.**(splev(lmhalo,lmh_spline) - lmhalo)


