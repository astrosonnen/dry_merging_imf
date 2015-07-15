import numpy as np
import pylab
import shmrs
import galaxies
import read_outputs
from matplotlib import rc
rc('text',usetex=True)

Ngal = 10
zref = 1.8

lmhalos = np.random.normal(12.7,0.3,Ngal)
lmstars = shmrs.mstarfunc_z0(lmhalos)

gals = read_outputs.read_evolution()

mygals = []

for i in range(0,Ngal):
    print i,lmhalos[i],lmstars[i]
    #gal = galaxies.etg(mstar_chab_0=10.**lmstars[i],mhalo_0 = 10.**lmhalos[i])
    gal = galaxies.etg(mstar_chab_0=gals[i+1]['mstar_chab'][0],mhalo_0 =gals[i+1]['Mh'][0])
    gal.evolve_back()

    if i==0:
        pylab.plot(gal.z,np.log10(gal.mstar_chab),color='r',label='Chab')
        pylab.plot(gal.z,np.log10(gal.mstar_true),color='g',label='True')
    else:
        pylab.plot(gal.z,np.log10(gal.mstar_chab),color='r')
        pylab.plot(gal.z,np.log10(gal.mstar_true),color='g')

    #gals.append(gal)
    mygals.append(gal)

pylab.xlabel('$z$',fontsize=20)
pylab.ylabel('$\log{M_*}$',fontsize=20)
pylab.ylim(10.,12.)
pylab.xticks(fontsize=20)
pylab.yticks(fontsize=20)
pylab.legend()
pylab.savefig('mstar_vs_z.png')
pylab.show()

done = False
for gal in mygals:
    zind = abs(gal.z-zref).argmin()
    if not done:
        pylab.scatter(np.log10(gal.mstar_chab[0]),np.log10(gal.mstar_true[0])-np.log10(gal.mstar_chab[0]),color='r',label='$z=0$')
        pylab.scatter(np.log10(gal.mstar_chab[zind]),np.log10(gal.mstar_true[zind])-np.log10(gal.mstar_chab[zind]),color='g',label='$z=1.8$')
        pylab.scatter(np.log10(gal.mstar_chab[-1]),np.log10(gal.mstar_true[-1]) - np.log10(gal.mstar_chab[-1]),color='b',label='$z=z_{\mathrm{form}}$')

        done = True
    else:
        pylab.scatter(np.log10(gal.mstar_chab[0]),np.log10(gal.mstar_true[0])-np.log10(gal.mstar_chab[0]),color='r')
        pylab.scatter(np.log10(gal.mstar_chab[zind]),np.log10(gal.mstar_true[zind])-np.log10(gal.mstar_chab[zind]),color='g')
        pylab.scatter(np.log10(gal.mstar_chab[-1]),np.log10(gal.mstar_true[-1]) - np.log10(gal.mstar_chab[-1]),color='b')


pylab.xlabel('$\log{M_*^{\mathrm{(Chab)}}}$',fontsize=20)
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$',fontsize=20)
pylab.xlim(10.,12.)
pylab.xticks(fontsize=20)
pylab.yticks(fontsize=20)
pylab.legend()
pylab.savefig('aimf_vs_mstar.png')
pylab.show()


