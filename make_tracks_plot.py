import pickle
import numpy as np
import pylab
import do_measurements as dm
from matplotlib import rc
rc('text', usetex=True)


snaps = [199, 100, 0]
markers = ['s', '^', 'o']
labels = ['$z=2$', '$z=1$', '$z=0$']
colors = ['b', 'g', 'r']
nsnap = len(snaps)

msize=50

# does the fit at each timestep, to follow the time evolution on the scaling relations in better detail
f = open('mstar_dep_imf_coeff0.5.dat', 'r')
galaxies = pickle.load(f)
f.close()

Ngal = len(galaxies)
lmstar = np.empty(Ngal)
lreff = np.empty(Ngal)
lsigma = np.empty(Ngal)
laimf = np.empty(Ngal)

fig = pylab.figure()
ax = fig.add_subplot(111)
pylab.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99)

mspar = []
vdpar = []
for j in range(0,200):

    for i in range(0, Ngal):
        lmstar[i] = np.log10(galaxies[i].mstar_chab[j])
        lreff[i] = np.log10(galaxies[i].re[j])
        lsigma[i] = np.log10(galaxies[i].veldisp[j])
        laimf[i] = np.log10(galaxies[i].aimf[j])

    pars, scat = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf)
    mspar.append(pars[1])
    vdpar.append(pars[2])

pylab.plot(mspar, vdpar, color='k')
for i in range(0, nsnap):
    pylab.scatter(mspar[snaps[i]], vdpar[snaps[i]], color=colors[i], marker=markers[i], s=msize)



f = open('vdisp_dep_imf_coeff2.0.dat', 'r')
galaxies = pickle.load(f)
f.close()

Ngal = len(galaxies)
lmstar = np.empty(Ngal)
lreff = np.empty(Ngal)
lsigma = np.empty(Ngal)
laimf = np.empty(Ngal)

mspar = []
vdpar = []
for j in range(0,200):

    for i in range(0, Ngal):
        lmstar[i] = np.log10(galaxies[i].mstar_chab[j])
        lreff[i] = np.log10(galaxies[i].re[j])
        lsigma[i] = np.log10(galaxies[i].veldisp[j])
        laimf[i] = np.log10(galaxies[i].aimf[j])

    pars, scat = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf)
    mspar.append(pars[1])
    vdpar.append(pars[2])

pylab.plot(mspar, vdpar, color='k')

for i in range(0, nsnap):
    pylab.scatter(mspar[snaps[i]], vdpar[snaps[i]], color=colors[i], marker=markers[i], s=msize)


#f = open('mstar-vdisp_dep_imf_coeff0.31.5.dat', 'r')
f = open('mstar-vdisp_dep_imf_coeff-0.31.5.dat', 'r')
#f = open('mhalo_dep_imf_coeff0.3.dat', 'r')
#f = open('mstar-wscatter_dep_imf_coeff0.5.dat', 'r')
galaxies = pickle.load(f)
f.close()

Ngal = len(galaxies)
lmstar = np.empty(Ngal)
lreff = np.empty(Ngal)
lsigma = np.empty(Ngal)
laimf = np.empty(Ngal)

mspar = []
vdpar = []
for j in range(0,200):

    for i in range(0, Ngal):
        lmstar[i] = np.log10(galaxies[i].mstar_chab[j])
        lreff[i] = np.log10(galaxies[i].re[j])
        lsigma[i] = np.log10(galaxies[i].veldisp[j])
        laimf[i] = np.log10(galaxies[i].aimf[j])

    pars, scat = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf)
    mspar.append(pars[1])
    vdpar.append(pars[2])

pylab.plot(mspar, vdpar, color='k')

for i in range(0, nsnap):
    pylab.scatter(mspar[snaps[i]], vdpar[snaps[i]], color=colors[i], marker=markers[i], s=msize, label=labels[i])

pylab.xlabel('$\partial \log{\\alpha_{\mathrm{IMF}}} / \partial \log{M_*}$', fontsize=16)
pylab.ylabel('$\partial \log{\\alpha_{\mathrm{IMF}}} / \partial \log{\sigma}$', fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
xticks = ax.xaxis.get_major_ticks()
xticks[0].label1.set_visible(False)
xticks[-2].label1.set_visible(False)

yticks = ax.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
yticks[-1].label1.set_visible(False)
pylab.text(-0.05, 2.1, '$\sigma$ model', fontsize=16)
pylab.text(0.45, 0.3, '$M_*$ model', fontsize=16)
pylab.text(0.25, 1.7, '$M_H$ model', fontsize=16)

pylab.legend(scatterpoints=1, loc='upper right')
pylab.savefig('tracks.eps')

pylab.show()


