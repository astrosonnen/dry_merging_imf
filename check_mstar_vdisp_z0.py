import pickle
import pylab
import numpy as np
import recipes
from plotters import rgb_alpha
from matplotlib import rc
rc('text', usetex=True)

#f = open('pop_vdisp1.10_model.dat', 'r')
#f = open('pop_vdisp_constsigma_model.dat', 'r')
f = open('pop_vdisp_model.dat', 'r')
pop = pickle.load(f)
f.close()

nobj = pop.mhalo.shape[0]

ms0 = []
vd0 = []
ms2 = []
vd2 = []

for i in range(nobj):
    ms0.append(np.log10(pop.mstar_salp[i, 20]))
    vd0.append(np.log10(pop.veldisp[i, 20]))
    ms2.append(np.log10(pop.mstar_salp[i, -1]))
    vd2.append(np.log10(pop.veldisp[i, -1]))
    #pylab.plot(np.log10(pop.mstar_salp[i, :]), np.log10(pop.veldisp[i, :]), color='gray', alpha=0.5)

pylab.xlim(10.5,12.5)

xs = np.linspace(10.5, 12.5)

auger_line = 2.34 + 0.18*(xs - 11.25)

mas_best = auger_line + 0.20*(np.log10(3.) - np.log10(1.2))
mas_up = auger_line + 0.27*(np.log10(3.) - np.log10(1.2))
mas_dw = auger_line + 0.13*(np.log10(3.) - np.log10(1.2))

pylab.fill_between(xs, mas_dw, mas_up, color=(0.7, 0.7, 0.7))
pylab.plot(xs, mas_best, color='b', label='Mason+15 ($z=2$)')
pylab.plot(xs, auger_line, color='r', linestyle='--', label='Auger+10 ($z=0.2$)')

pylab.scatter(ms2, vd2, color='b', label='$z=2$', marker='s', s=30)
pylab.scatter(ms0, vd0, color='r', label='$z=0.2$', s=30)


pylab.legend(scatterpoints=1, loc='upper left')

pylab.xlabel('$\log{M_*^{\mathrm{Salp}}}$', fontsize=16)
pylab.ylabel('$\log{\sigma}$', fontsize=16)
pylab.savefig('mstar-vdisp_evolution.eps')
pylab.show()

