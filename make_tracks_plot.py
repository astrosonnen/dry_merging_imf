import pickle
import numpy as np
import pylab
import do_measurements as dm
from plotters import probcontour
from matplotlib import rc
rc('text', usetex=True)

snaps = [199, 100, 30]
markers = ['s', '^', 'o']
labels = ['$z=2$', '$z=1$', '$z=0.3$']
colors = ['b', 'g', 'r']
nsnap = len(snaps)

f = open('/gdrive/projects/SL2S_hierarch/evol_nfw_alpha_only_msps_sigma.dat', 'r')
chain = pickle.load(f)
f.close()
burnin = 10000
mdep_chain = chain['amstar'][burnin:]
sdep_chain = chain['asigma'][burnin:]

fsize=20
lsize=24
msize=100

# does the fit at each timestep, to follow the time evolution on the scaling relations in better detail
f = open('pop_mstar_model.dat', 'r')
pop = pickle.load(f)
f.close()

Ngal = 100
lmstar = np.empty(Ngal)
lreff = np.empty(Ngal)
lsigma = np.empty(Ngal)
laimf = np.empty(Ngal)

fig = pylab.figure()
ax = fig.add_subplot(111)
pylab.subplots_adjust(left=0.15, right=0.99, bottom=0.15, top=0.99)


probcontour(mdep_chain, sdep_chain, style=(0., 0.7, 1.))
probcontour(mdep_chain, sdep_chain, style='black')

mspar = []
vdpar = []
for j in range(0, 200):

    lmstar = np.log10(pop.mstar_chab[:, j])
    lsigma = np.log10(pop.veldisp[:, j])
    laimf = np.log10(pop.aimf[:, j])

    """
    pars, scat = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf)

    mspar.append(pars[1])
    vdpar.append(pars[2])
    """
    mspar.append(pop.imf_coeff[j, 1])
    vdpar.append(pop.imf_coeff[j, 2])

pylab.plot(mspar[30:], vdpar[30:], color='k')
for i in range(0, nsnap):
    pylab.scatter(mspar[snaps[i]], vdpar[snaps[i]], color=colors[i], marker=markers[i], s=msize)



f = open('pop_vdisp_model.dat', 'r')
pop = pickle.load(f)
f.close()

Ngal = 100
lmstar = np.empty(Ngal)
lsigma = np.empty(Ngal)
laimf = np.empty(Ngal)

mspar = []
vdpar = []
for j in range(0, 200):

    lmstar = np.log10(pop.mstar_chab[:, j])
    lsigma = np.log10(pop.veldisp[:, j])
    laimf = np.log10(pop.aimf[:, j])

    """
    pars, scat = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf, guess=(0.,0.,2.3))
    mspar.append(pars[1])
    vdpar.append(pars[2])
    """

    mspar.append(pop.imf_coeff[j, 1])
    vdpar.append(pop.imf_coeff[j, 2])

pylab.plot(mspar[30:], vdpar[30:], color='k')

for i in range(0, nsnap):
    pylab.scatter(mspar[snaps[i]], vdpar[snaps[i]], color=colors[i], marker=markers[i], s=msize)


"""
f = open('pop_mhalo_model_0.00_0.30.dat', 'r')
pop = pickle.load(f)
f.close()

Ngal = 100
lmstar = np.empty(Ngal)
lreff = np.empty(Ngal)
lsigma = np.empty(Ngal)
laimf = np.empty(Ngal)

mspar = []
vdpar = []
for j in range(0,200):

    lmstar = np.log10(pop.mstar_chab[:, j])
    lsigma = np.log10(pop.veldisp[:, j])
    laimf = np.log10(pop.aimf[:, j])


    pars, scat = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf)
    mspar.append(pars[1])
    vdpar.append(pars[2])

pylab.plot(mspar, vdpar, color='k')

for i in range(0, nsnap):
    pylab.scatter(mspar[snaps[i]], vdpar[snaps[i]], color=colors[i], marker=markers[i], s=msize, label=labels[i])
"""

pylab.xlabel('$a_*$', fontsize=lsize)# ($\partial \log{\\alpha_{\mathrm{IMF}}} / \partial \log{M_*}$)', fontsize=fsize)
pylab.ylabel('$a_\sigma$', fontsize=lsize) #(\partial \log{\\alpha_{\mathrm{IMF}}} / \partial \log{\sigma}$)', fontsize=fsize)
pylab.xticks(fontsize=fsize)
pylab.yticks(fontsize=fsize)
xticks = ax.xaxis.get_major_ticks()
xticks[0].label1.set_visible(False)
xticks[-2].label1.set_visible(False)

yticks = ax.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
yticks[-1].label1.set_visible(False)
pylab.text(-0.02, 1.6, '$\sigma$ model', fontsize=fsize)
pylab.text(0.27, -0.2, '$M_*$ model', fontsize=fsize)
#pylab.text(0.17, 1.6, 'Halo model', fontsize=fsize)


pylab.legend(scatterpoints=1, loc='upper right', fontsize=fsize)
pylab.savefig('tracks.eps')

pylab.show()


