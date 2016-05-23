import pickle
import numpy as np
import pylab
import do_measurements as dm
from plotters import rgb_alpha

bandcolor = rgb_alpha((0,255,255), 1.)

f = open('pop_mstar_model_0.11_0.26.dat', 'r')
ms_pop = pickle.load(f)
f.close()

ngal = 100

f = open('pop_vdisp_model_-0.06_1.00.dat', 'r')
vd_pop = pickle.load(f)
f.close()

snaps = [199, 100, 30]
markers = ['s', '^', 'o']
labels = ['$z=2$', '$z=1$', '$z=0.3$']
colors = ['b', 'g', 'r']
nsnap = len(snaps)

fitdir = '/gdrive/projects/SL2S_hierarch/'

f = open(fitdir+'evol_nfw_alpha_only_msps_only.dat', 'r')
msps_chain = pickle.load(f)
f.close()

Nsamp = 1000
indices = np.arange(20000, 100000)
samp = np.random.choice(indices, Nsamp)
mspsfit_a_samp = msps_chain['amstar'][samp]
mspsfit_b_samp = msps_chain['calpha'][samp]

f = open(fitdir+'evol_nfw_alpha_only_sigma_only.dat', 'r')
sigma_chain = pickle.load(f)
f.close()

sigmafit_piv = np.log10(250.)
sigmafit_a_samp = sigma_chain['asigma'][samp]
sigmafit_b_samp = sigma_chain['calpha'][samp]

xlimm = (10.5,12.5)
xlimv = (2.2,2.6)

xsm = np.linspace(xlimm[0], xlimm[1],51)

#makes observational band
m84 = 0.*xsm
m16 = 0.*xsm
for i in range(0,51):
    aimf_samp = mspsfit_b_samp + mspsfit_a_samp*(xsm[i] - 11.5)
    m84[i] = np.percentile(aimf_samp, 84.)
    m16[i] = np.percentile(aimf_samp, 16.)

#makes observational band

xsv = np.linspace(xlimv[0], xlimv[1], 51)
v84 = 0.*xsv
v16 = 0.*xsv
for i in range(0,51):
    aimf_samp = sigmafit_b_samp + sigmafit_a_samp*(xsv[i] - sigmafit_piv)
    v84[i] = np.percentile(aimf_samp, 84.)
    v16[i] = np.percentile(aimf_samp, 16.)



fig = pylab.figure()
ax1 = fig.add_subplot(2, 2, 1)
pylab.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99, hspace=0., wspace=0.)
pylab.fill_between(xsm, m16, m84, color=bandcolor)

fitpars = []
for i in range(0, nsnap):
    mchab_ms = np.log10(ms_pop.mstar_salp[:, snaps[i]])
    aimf_ms = np.log10(ms_pop.aimf[:, snaps[i]])
    pylab.scatter(mchab_ms, aimf_ms, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_mstar_only(mchab_ms, aimf_ms)
    print i, par, scat
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
#yticks[0].label1.set_visible(False)
yticks[1].label1.set_visible(False)
yticks[-2].label1.set_visible(False)
pylab.tick_params(axis='x', labelbottom='off')

fitpars = []
pylab.subplot(2,2,2)
pylab.fill_between(xsv, v16, v84, color=bandcolor)
for i in range(0, nsnap):
    vdisp_ms = np.log10(ms_pop.veldisp[:, snaps[i]])
    aimf_ms = np.log10(ms_pop.aimf[:, snaps[i]])
    pylab.scatter(vdisp_ms, aimf_ms, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_sigma_only(vdisp_ms, aimf_ms)
    print scat
    fitpars.append(par)

xlim = pylab.xlim()
#ylim = pylab.ylim()
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

pylab.fill_between(xsm, m16, m84, color=bandcolor)

fitpars = []
for i in range(0, nsnap):
    mchab_vd = np.log10(vd_pop.mstar_salp[:, snaps[i]])
    aimf_vd = np.log10(vd_pop.aimf[:, snaps[i]])
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
yticks[-1].label1.set_visible(False)
pylab.text(11.8, -0.3, '$\sigma$ model', fontsize=14)

xticks = ax3.xaxis.get_major_ticks()
#xticks[0].label1.set_visible(False)
xticks[-1].label1.set_visible(False)


ax4 = pylab.subplot(2, 2, 4)
pylab.fill_between(xsv, v16, v84, color=bandcolor)

fitpars = []
for i in range(0, nsnap):
    vdisp_vd = np.log10(vd_pop.veldisp[:, snaps[i]])
    aimf_vd = np.log10(vd_pop.aimf[:, snaps[i]])
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
xticks[-1].label1.set_visible(False)


pylab.savefig('snapshots.eps')


pylab.show()





