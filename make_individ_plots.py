import pickle
import numpy as np
import pylab
#import do_measurements as dm
#from matplotlib import rc
#rc('text', usetex=True)

fsize=20
lsize=24
msize=100

model = 'vdisp'
# does the fit at each timestep, to follow the time evolution on the scaling relations in better detail
f = open('pop_%s_model_-0.22_1.60.dat'%model, 'r')
pop = pickle.load(f)
f.close()

colors = ['b', 'g', 'r', 'y', 'c', 'm', 'k', 'purple', 'orange', 'gray', 'brown']
nobj = len(colors)

#pylab.subplot(4, 1, 1)
for i in range(0, nobj):
    pylab.plot(pop.z, np.log10(pop.mstar_chab[i, :]), color=colors[i])
pylab.ylabel('$\log{M_*}$')
pylab.savefig(model + '_model_mstar_evol.png')
pylab.close()

#pylab.subplot(4, 1, 2)
for i in range(0, nobj):
    pylab.plot(pop.z, np.log10(pop.mhalo[i, :]), color=colors[i])
pylab.ylabel('$\log{M_h}$')
pylab.savefig(model + '_model_mhalo_evol.png')
pylab.close()

#pylab.subplot(4, 1, 3)
for i in range(0, nobj):
    pylab.plot(pop.z, pop.veldisp[i, :], color=colors[i])
pylab.ylabel('$\sigma$')
pylab.savefig(model + '_model_vdisp_evol.png')
pylab.close()

#pylab.subplot(4, 1, 4)
for i in range(0, nobj):
    pylab.plot(pop.z, np.log10(pop.aimf[i, :]), color=colors[i])
pylab.xlabel('$z$')
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$')
pylab.savefig(model + '_model_aimf_evol.png')
pylab.close()


#pylab.show()


