import pickle
import numpy as np
import pylab
import do_measurements as dm
from matplotlib import rc
rc('text', usetex=True)


f = open('mstar_dep_imf_coeff0.5.dat', 'r')
ms_galaxies = pickle.load(f)
f.close()

ngal = 100

f = open('vdisp_dep_imf_coeff2.0.dat', 'r')
vd_galaxies = pickle.load(f)
f.close()

snaps = [0, 100, 199]
markers = ['s', '^', 'o']
labels = ['$z=2$', '$z=1$', '$z=0$']
colors = ['b', 'g', 'r']
nsnap = len(snaps)

fig = pylab.figure()
ax1 = fig.add_subplot(2, 2, 1)
pylab.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99, hspace=0., wspace=0.)
for i in range(0, nsnap):
    mchab_ms = np.zeros(ngal)
    vdisp_ms = np.zeros(ngal)
    aimf_ms = np.zeros(ngal)
    for j in range(0, ngal):
	mchab_ms[j] = np.log10(ms_galaxies[j].mstar_chab[snaps[i]])
	aimf_ms[j] = np.log10(ms_galaxies[j].aimf[snaps[i]])
    pylab.scatter(mchab_ms, aimf_ms, color=colors[i], s=30, marker=markers[i], label=labels[i])
pylab.text(12., 0.5, '$M_*$ model')
pylab.xlabel('$\log{M_*}$', fontsize=14)
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.legend(scatterpoints=1, fontsize=14, loc='upper left')
#pylab.xticks(fontsize=14)
pylab.yticks(fontsize=14)
yticks = ax1.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
yticks[-2].label1.set_visible(False)
pylab.tick_params(axis='x', labelbottom='off')

pylab.subplot(2,2,2)
for i in range(0, nsnap):
    vdisp_ms = np.zeros(ngal)
    aimf_ms = np.zeros(ngal)
    for j in range(0, ngal):
	vdisp_ms[j] = np.log10(ms_galaxies[j].veldisp[snaps[i]])
	aimf_ms[j] = np.log10(ms_galaxies[j].aimf[snaps[i]])
    pylab.scatter(vdisp_ms, aimf_ms, color=colors[i], s=30, marker=markers[i], label=labels[i])

pylab.xlabel('$\log{\sigma}$', fontsize=14)
#pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.xticks(fontsize=14)
#pylab.yticks(fontsize=16)
pylab.tick_params(axis='y', labelleft='off')
pylab.tick_params(axis='x', labelbottom='off')
pylab.legend(scatterpoints=1, fontsize=14, loc='upper left')


ax3 = pylab.subplot(2, 2, 3)

for i in range(0, nsnap):
    mchab_vd = np.zeros(ngal)
    vdisp_vd = np.zeros(ngal)
    aimf_vd = np.zeros(ngal)
    for j in range(0, ngal):
	mchab_vd[j] = np.log10(vd_galaxies[j].mstar_chab[snaps[i]])
	aimf_vd[j] = np.log10(vd_galaxies[j].aimf[snaps[i]])
    pylab.scatter(mchab_vd, aimf_vd, color=colors[i], s=30, marker=markers[i], label=labels[i])
pylab.xlabel('$\log{M_*}$', fontsize=14)
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.legend(scatterpoints=1, fontsize=14, loc='upper left')
pylab.xticks(fontsize=14)
pylab.yticks(fontsize=14)
yticks = ax3.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
yticks[-2].label1.set_visible(False)

xticks = ax3.yaxis.get_major_ticks()
xticks[0].label1.set_visible(False)
xticks[-2].label1.set_visible(False)

pylab.show()





# does the fit at each timestep, to follow the time evolution on the scaling relations in better detail
f = open('mstar_dep_imf_coeff0.5.dat', 'r')
galaxies = pickle.load(f)
f.close()

Ngal = len(galaxies)
lmstar = np.empty(Ngal)
lreff = np.empty(Ngal)
lsigma = np.empty(Ngal)
laimf = np.empty(Ngal)

for j in range(0,200):

    for i in range(0, Ngal):
        lmstar[i] = np.log10(galaxies[i].mstar_chab[j])
        lreff[i] = np.log10(galaxies[i].re[j])
        lsigma[i] = np.log10(galaxies[i].veldisp[j])
        laimf[i] = np.log10(galaxies[i].aimf[j])

    pars_mstar_sigma = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf)

    pylab.scatter(pars_mstar_sigma[0][1], pars_mstar_sigma[0][2], color='r')

f = open('vdisp_dep_imf_coeff2.0.dat', 'r')
galaxies = pickle.load(f)
f.close()

Ngal = len(galaxies)
lmstar = np.empty(Ngal)
lreff = np.empty(Ngal)
lsigma = np.empty(Ngal)
laimf = np.empty(Ngal)

for j in range(0,200):

    for i in range(0, Ngal):
        lmstar[i] = np.log10(galaxies[i].mstar_chab[j])
        lreff[i] = np.log10(galaxies[i].re[j])
        lsigma[i] = np.log10(galaxies[i].veldisp[j])
        laimf[i] = np.log10(galaxies[i].aimf[j])

    pars_mstar_sigma = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf)

    pylab.scatter(pars_mstar_sigma[0][1], pars_mstar_sigma[0][2], color='b')

f = open('mstar-vdisp_dep_imf_coeff0.31.5.dat', 'r')
galaxies = pickle.load(f)
f.close()

Ngal = len(galaxies)
lmstar = np.empty(Ngal)
lreff = np.empty(Ngal)
lsigma = np.empty(Ngal)
laimf = np.empty(Ngal)

for j in range(0,200):

    for i in range(0, Ngal):
        lmstar[i] = np.log10(galaxies[i].mstar_chab[j])
        lreff[i] = np.log10(galaxies[i].re[j])
        lsigma[i] = np.log10(galaxies[i].veldisp[j])
        laimf[i] = np.log10(galaxies[i].aimf[j])

    pars_mstar_sigma = dm.fit_mstar_sigma_fixed_z(lmstar, lsigma, laimf)

    pylab.scatter(pars_mstar_sigma[0][1], pars_mstar_sigma[0][2], color='m')

pylab.show()




mstar_fits = dm.do_fits('mstar_dep_imf_coeff0.5.dat')
vdisp_fits = dm.do_fits('vdisp_dep_imf_coeff2.0.dat')
both_fits = dm.do_fits('mstar-vdisp_dep_imf_coeff0.52.0.dat')

mstar_model_z2 = (mstar_fits['z=2']['mstar-sigma'][0][1], mstar_fits['z=2']['mstar-sigma'][0][2])
vdisp_model_z2 = (vdisp_fits['z=2']['mstar-sigma'][0][1], vdisp_fits['z=2']['mstar-sigma'][0][2])
both_model_z2 = (both_fits['z=2']['mstar-sigma'][0][1], both_fits['z=2']['mstar-sigma'][0][2])

mstar_model_z1 = (mstar_fits['z=1']['mstar-sigma'][0][1], mstar_fits['z=1']['mstar-sigma'][0][2])
vdisp_model_z1 = (vdisp_fits['z=1']['mstar-sigma'][0][1], vdisp_fits['z=1']['mstar-sigma'][0][2])
both_model_z1 = (both_fits['z=1']['mstar-sigma'][0][1], both_fits['z=1']['mstar-sigma'][0][2])

mstar_model_z0 = (mstar_fits['z=0']['mstar-sigma'][0][1], mstar_fits['z=0']['mstar-sigma'][0][2])
vdisp_model_z0 = (vdisp_fits['z=0']['mstar-sigma'][0][1], vdisp_fits['z=0']['mstar-sigma'][0][2])
both_model_z0 = (both_fits['z=0']['mstar-sigma'][0][1], both_fits['z=0']['mstar-sigma'][0][2])

pylab.scatter(mstar_model_z2[0], mstar_model_z2[1], color='r', s=40, label='$M_*$ model, $z=2$')
pylab.scatter(vdisp_model_z2[0], vdisp_model_z2[1], s=40, label='$\sigma$ model, $z=2$')
pylab.scatter(both_model_z2[0], both_model_z2[1], color='g', s=40, label='$M_*-\sigma$ model, $z=2$')

pylab.scatter(mstar_model_z1[0], mstar_model_z1[1], color='r', s=40, label='$M_*$ model, $z=1$', marker='s')
pylab.scatter(vdisp_model_z1[0], vdisp_model_z1[1], s=40, label='$\sigma$ model, $z=1$', marker='s')
pylab.scatter(both_model_z1[0], both_model_z1[1], color='g', s=40, label='$M_*-\sigma$ model, $z=1$', marker='s')

pylab.scatter(mstar_model_z0[0], mstar_model_z0[1], color='r', s=40, label='$M_*$ model, $z=0$', marker='^')
pylab.scatter(vdisp_model_z0[0], vdisp_model_z0[1], s=40, label='$\sigma$ model, $z=0$', marker='^')
pylab.scatter(both_model_z0[0], both_model_z0[1], color='g', s=40, label='$M_*-\sigma$ model, $z=0$', marker='^')
#pylab.legend(scatterpoints=1)
pylab.xlabel('$\partial \log{\\alpha_{\mathrm{IMF}}} / \partial \log{M_*}$', fontsize=16)
pylab.ylabel('$\partial \log{\\alpha_{\mathrm{IMF}}} / \partial \log{\sigma}$', fontsize=16)
pylab.savefig('dep_evolution_plot.png')
pylab.show()


