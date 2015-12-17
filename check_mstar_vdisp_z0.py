import pickle
import pylab
import numpy as np
import recipes

f = open('vdisp_dep_imf_coeff2.30.dat', 'r')
gals = pickle.load(f)
f.close()

ms0 = []
vd0 = []
ms2 = []
vd2 = []

for gal in gals:
    ms0.append(np.log10(gal.mstar_chab[0]))
    vd0.append(np.log10(gal.veldisp[0]))
    ms2.append(np.log10(gal.mstar_chab[-1]))
    vd2.append(np.log10(gal.veldisp[-1]))
    pylab.plot(np.log10(gal.mstar_chab), np.log10(gal.veldisp), color='gray', alpha=0.5)

pylab.scatter(ms2, vd2, color='b', label='$z=2$')
pylab.scatter(ms0, vd0, color='r', label='$z=0$')

pylab.xlim(10.,12.5)

xs = np.linspace(10., 12.5)
mas2 = np.log10(recipes.vdisp_mstar_rel_mason(xs, 2.))
mas0 = np.log10(recipes.vdisp_mstar_rel_mason(xs, 0.))

pylab.plot(xs, mas2, color='b', linestyle='--', label='Mason+15 ($z=2$)')
pylab.plot(xs, mas0, color='r', linestyle='--', label='Mason+15 ($z=0$)')

pylab.legend(scatterpoints=1, loc='upper left')

pylab.xlabel('$\log{M_*}$', fontsize=16)
pylab.ylabel('$\log{\sigma}$', fontsize=16)
pylab.savefig('mstar-vdisp_evolution.png')
pylab.show()



