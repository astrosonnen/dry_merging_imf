import pickle
import numpy as np
import pylab
import do_measurements as dm

snaps = [199, 100, 0]
markers = ['s', '^', 'o']
labels = ['$z=2$', '$z=1$', '$z=0$']
colors = ['b', 'g', 'r']
nsnap = len(snaps)

fsize=20
lsize=24
msize=100

# does the fit at each timestep, to follow the time evolution on the scaling relations in better detail
f = open('mstar_dep_imf_coeff0.33.dat', 'r')
#f = open('noscatter_vdisp_dep_imf_coeff2.3.dat', 'r')
#f = open('noscatter_mstar_dep_imf_coeff0.5.dat', 'r')
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



#f = open('vdisp_dep_imf_coeff2.3.dat', 'r')
f = open('true_noscatter_vdisp_dep_imf_coeff2.30.dat', 'r')
#f = open('vdisp_dep_imf_coeff2.0_mhmax13.dat', 'r')
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

    pars, scat = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf, guess=(0.,0.,2.3))
    mspar.append(pars[1])
    vdpar.append(pars[2])

pylab.plot(mspar, vdpar, color='k')

for i in range(0, nsnap):
    pylab.scatter(mspar[snaps[i]], vdpar[snaps[i]], color=colors[i], marker=markers[i], s=msize)


f = open('mstar-vdisp_dep_imf_coeff0.21.5.dat', 'r')
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
pylab.text(-0.02, 2.1, '$\sigma$ model', fontsize=fsize)
pylab.text(0.27, -0.2, '$M_*$ model', fontsize=fsize)
pylab.text(0.17, 1.6, 'Hybrid model', fontsize=fsize)

pylab.legend(scatterpoints=1, loc='upper right', fontsize=fsize)
pylab.savefig('tracks.eps')

pylab.show()


