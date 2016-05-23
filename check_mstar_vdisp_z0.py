import pickle
import pylab
import numpy as np
import recipes

#f = open('pop_vdisp_model_-0.22_1.60.dat', 'r')
f = open('pop_mstar_model_0.11_0.28.dat', 'r')
#f = open('pop_mhalo_model_0.00_0.30.dat', 'r')
pop = pickle.load(f)
f.close()

nobj = pop.mhalo.shape[0]

ms0 = []
vd0 = []
ms2 = []
vd2 = []

for i in range(nobj):
    ms0.append(np.log10(pop.mstar_chab[i, 0]))
    vd0.append(np.log10(pop.veldisp[i, 0]))
    ms2.append(np.log10(pop.mstar_chab[i, -1]))
    vd2.append(np.log10(pop.veldisp[i, -1]))
    pylab.plot(np.log10(pop.mstar_chab[i, :]), np.log10(pop.veldisp[i, :]), color='gray', alpha=0.5)

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



