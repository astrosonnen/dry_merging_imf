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
coeff = '2.30'
# does the fit at each timestep, to follow the time evolution on the scaling relations in better detail
f = open('%s_dep_imf_coeff%s.dat'%(model, coeff), 'r')
galaxies = pickle.load(f)
f.close()

f = open('%s_dep_imf_coeff%s_it1.dat'%(model, coeff), 'r')
it1 = pickle.load(f)
f.close()

f = open('%s_dep_imf_coeff%s_it2.dat'%(model, coeff), 'r')
it2 = pickle.load(f)
f.close()

f = open('%s_dep_imf_coeff%s_it3.dat'%(model, coeff), 'r')
it3 = pickle.load(f)
f.close()

colors = ['b', 'g', 'r', 'y', 'c', 'm', 'k', 'purple', 'orange', 'gray', 'brown']
nobj = len(colors)

#pylab.subplot(4, 1, 1)
for i in range(0, nobj):
    pylab.plot(galaxies[i].z, np.log10(galaxies[i].mstar_chab), color=colors[i])
pylab.ylabel('$\log{M_*}$')
pylab.savefig(model + '_model_mstar_evol.png')
pylab.close()

#pylab.subplot(4, 1, 2)
for i in range(0, nobj):
    pylab.plot(galaxies[i].z, np.log10(galaxies[i].mhalo), color=colors[i])
pylab.ylabel('$\log{M_h}$')
pylab.savefig(model + '_model_mhalo_evol.png')
pylab.close()

#pylab.subplot(4, 1, 3)
for i in range(0, nobj):
    pylab.plot(galaxies[i].z, galaxies[i].veldisp, color=colors[i])
for i in range(nobj):
    pylab.plot(galaxies[i+nobj].z, galaxies[i+nobj].veldisp, color=colors[i])
pylab.ylabel('$\sigma$')
pylab.savefig(model + '_model_vdisp_evol.png')
pylab.close()

#pylab.subplot(4, 1, 4)
for i in range(0, nobj):
    pylab.plot(galaxies[i].z, np.log10(galaxies[i].aimf), color=colors[i])
    pylab.plot(it1[i].z, np.log10(it1[i].aimf), color=colors[i])
    pylab.plot(it2[i].z, np.log10(it2[i].aimf), color=colors[i])
    pylab.plot(it3[i].z, np.log10(it3[i].aimf), color=colors[i])
pylab.xlabel('$z$')
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$')
pylab.savefig(model + '_model_aimf_evol.png')
pylab.close()


#pylab.show()


