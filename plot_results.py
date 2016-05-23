import numpy as np
import pylab
import pickle
import sys
from matplotlib import rc
rc('text',usetex=True)


if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'output.dat'

f = open(filename,'r')
galaxies = pickle.load(f)
f.close()

Ngal = len(galaxies)

for i in range(0,Ngal):

    gal = galaxies[i]

    if i==0:
        pylab.plot(gal.z,np.log10(gal.mstar_chab),color='r',label='Chab')
        pylab.plot(gal.z,np.log10(gal.mstar_true),color='g',label='True')
    else:
        pylab.plot(gal.z,np.log10(gal.mstar_chab),color='r')
        pylab.plot(gal.z,np.log10(gal.mstar_true),color='g')

pylab.xlabel('$z$',fontsize=20)
pylab.ylabel('$\log{M_*}$',fontsize=20)
pylab.ylim(10.,12.)
pylab.xticks(fontsize=20)
pylab.yticks(fontsize=20)
pylab.legend()
pylab.savefig('mstar_vs_z.png')
pylab.show()

done = False
for gal in galaxies:
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


