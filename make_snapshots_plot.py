import pickle
import numpy as np
import pylab
import do_measurements as dm
from matplotlib import rc
rc('text', usetex=True)


#f = open('noscatter_mstar_dep_imf_coeff0.5.dat', 'r')
#f = open('mstar_dep_imf_coeff0.4.dat', 'r')
f = open('mstar_dep_imf_coeff0.3.dat', 'r')
#f = open('mstar_dep_imf_coeff0.5_mhmax13.dat', 'r')
#f = open('mstar-wscatter_dep_imf_coeff0.5.dat', 'r')
ms_galaxies = pickle.load(f)
f.close()

ngal = 100

f = open('vdisp_dep_imf_coeff2.3.dat', 'r')
#f = open('vdisp_dep_imf_coeff2.0_mhmax13.dat', 'r')
#f = open('mhalo_dep_imf_coeff0.3.dat', 'r')
#f = open('mstar-vdisp_dep_imf_coeff-0.31.5.dat', 'r')
vd_galaxies = pickle.load(f)
f.close()

snaps = [199, 100, 0]
markers = ['s', '^', 'o']
labels = ['$z=2$', '$z=1$', '$z=0$']
colors = ['b', 'g', 'r']
nsnap = len(snaps)

son15_a = 0.20
son15_a_err = 0.04
son15_b = 0.05
son15_b_err = 0.02
Nsamp = 1000
son15_a_samp = np.random.normal(son15_a,son15_a_err,Nsamp)
son15_b_samp = np.random.normal(son15_b,son15_b_err,Nsamp)

pos15_a = 1.30
pos15_a_err = 0.23
pos15_b = -0.14
pos15_b_err = 0.03
pos15_piv = np.log10(200.)
pos15_a_samp = np.random.normal(pos15_a,pos15_a_err,Nsamp)
pos15_b_samp = np.random.normal(pos15_b,pos15_b_err,Nsamp)

xlimm = (10.5,12.5)
xlimv = (2.2,2.6)

xsm = np.linspace(xlimm[0], xlimm[1],51)

#makes observational band
m84 = 0.*xsm
m16 = 0.*xsm
for i in range(0,51):
    aimf_samp = son15_b_samp + son15_a_samp*(xsm[i] - 11.5)
    m84[i] = np.percentile(aimf_samp, 84.)
    m16[i] = np.percentile(aimf_samp, 16.)

#makes observational band

xsv = np.linspace(xlimv[0], xlimv[1], 51)
v84 = 0.*xsv
v16 = 0.*xsv
for i in range(0,51):
    aimf_samp = pos15_b_samp + pos15_a_samp*(xsv[i] - pos15_piv)
    v84[i] = np.percentile(aimf_samp, 84.)
    v16[i] = np.percentile(aimf_samp, 16.)



fig = pylab.figure()
ax1 = fig.add_subplot(2, 2, 1)
pylab.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99, hspace=0., wspace=0.)
pylab.fill_between(xsm, m16, m84, color='c')

fitpars = []
for i in range(0, nsnap):
    mchab_ms = np.zeros(ngal)
    vdisp_ms = np.zeros(ngal)
    aimf_ms = np.zeros(ngal)
    for j in range(0, ngal):
	mchab_ms[j] = np.log10(ms_galaxies[j].mstar_chab[snaps[i]])
	aimf_ms[j] = np.log10(ms_galaxies[j].aimf[snaps[i]])
    pylab.scatter(mchab_ms, aimf_ms, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_mstar_only(mchab_ms, aimf_ms)
    print scat
    fitpars.append(par)

xlim = pylab.xlim()
ylim = pylab.ylim()

xs = np.linspace(xlim[0], xlim[1])

for i in range(0, nsnap):
    par = fitpars[i]
    pylab.plot(xs, par[0] + (xs - 11.)*par[1], linestyle='--', color=colors[i])

pylab.ylim(ylim)
pylab.text(11.8, -0.1, '$M_*$ model', fontsize=14)
#pylab.xlabel('$\log{M_*}$', fontsize=14)
pylab.xlim(10.5,12.5)
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.legend(scatterpoints=1, fontsize=14, loc='upper left')
#pylab.xticks(fontsize=14)
pylab.yticks(fontsize=14)
yticks = ax1.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
yticks[-2].label1.set_visible(False)
pylab.tick_params(axis='x', labelbottom='off')

fitpars = []
pylab.subplot(2,2,2)
pylab.fill_between(xsv, v16, v84, color='c')
for i in range(0, nsnap):
    vdisp_ms = np.zeros(ngal)
    aimf_ms = np.zeros(ngal)
    for j in range(0, ngal):
	vdisp_ms[j] = np.log10(ms_galaxies[j].veldisp[snaps[i]])
	aimf_ms[j] = np.log10(ms_galaxies[j].aimf[snaps[i]])
    pylab.scatter(vdisp_ms, aimf_ms, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_sigma_only(vdisp_ms, aimf_ms)
    print scat
    fitpars.append(par)

xlim = pylab.xlim()
ylim = pylab.ylim()
xs = np.linspace(xlim[0], xlim[1], 51)

for i in range(0, nsnap):
    par = fitpars[i]
    pylab.plot(xs, par[0] + (xs - 2.3)*par[1], linestyle='--', color=colors[i])
#pylab.xlim(xlim)
pylab.xlim(xlimv[0], xlimv[1])
pylab.ylim(ylim)



pylab.xlabel('$\log{\sigma}$', fontsize=14)
#pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.xticks(fontsize=14)
#pylab.yticks(fontsize=16)
pylab.text(2.5, -0.1, '$M_*$ model', fontsize=14)
pylab.tick_params(axis='y', labelleft='off')
pylab.tick_params(axis='x', labelbottom='off')


alphas = [0.2, 0.5, 1.0]

ax3 = pylab.subplot(2, 2, 3)

pylab.fill_between(xsm, m16, m84, color='c')

fitpars = []
for i in range(0, nsnap):
    mchab_vd = np.zeros(ngal)
    vdisp_vd = np.zeros(ngal)
    aimf_vd = np.zeros(ngal)
    for j in range(0, ngal):
	mchab_vd[j] = np.log10(vd_galaxies[j].mstar_chab[snaps[i]])
	aimf_vd[j] = np.log10(vd_galaxies[j].aimf[snaps[i]])
    pylab.scatter(mchab_vd, aimf_vd, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_mstar_only(mchab_vd, aimf_vd)
    print scat
    fitpars.append(par)

xlim = pylab.xlim()
ylim = pylab.ylim()
xs = np.linspace(xlim[0], xlim[1], 51)

for i in range(0, nsnap):
    par = fitpars[i]
    pylab.plot(xs, par[0] + (xs - 11.)*par[1], linestyle='--', color=colors[i])



pylab.ylim(ylim)

pylab.xlabel('$\log{M_*}$', fontsize=14)
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.xlim(10.5,12.5)
pylab.xticks(fontsize=14)
pylab.yticks(fontsize=14)
yticks = ax3.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
yticks[-2].label1.set_visible(False)
pylab.text(11.8, -0.3, '$\sigma$ model', fontsize=14)

xticks = ax3.xaxis.get_major_ticks()
#xticks[0].label1.set_visible(False)
xticks[-1].label1.set_visible(False)


ax4 = pylab.subplot(2, 2, 4)
pylab.fill_between(xsv, v16, v84, color='c')

fitpars = []
for i in range(0, nsnap):
    mchab_vd = np.zeros(ngal)
    vdisp_vd = np.zeros(ngal)
    aimf_vd = np.zeros(ngal)
    for j in range(0, ngal):
	vdisp_vd[j] = np.log10(vd_galaxies[j].veldisp[snaps[i]])
	aimf_vd[j] = np.log10(vd_galaxies[j].aimf[snaps[i]])
    pylab.scatter(vdisp_vd, aimf_vd, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_sigma_only(vdisp_vd, aimf_vd)
    print scat
    fitpars.append(par)

xlim = pylab.xlim()
ylim = pylab.ylim()
xs = np.linspace(xlim[0], xlim[1], 51)

for i in range(0, nsnap):
    par = fitpars[i]
    pylab.plot(xs, par[0] + (xs - 2.3)*par[1], linestyle='--', color=colors[i])
#pylab.xlim(xlim)
pylab.xlim(xlimv[0], xlimv[1])
pylab.ylim(ylim)


pylab.xlabel('$\log{\sigma}$', fontsize=14)
#pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.xticks(fontsize=14)
#pylab.yticks(fontsize=14)
pylab.text(2.5, -0.3, '$\sigma$ model', fontsize=14)

pylab.tick_params(axis='y', labelleft='off')
xticks = ax4.xaxis.get_major_ticks()
xticks[0].label1.set_visible(False)
xticks[-2].label1.set_visible(False)


pylab.savefig('snapshots.png')


pylab.show()





