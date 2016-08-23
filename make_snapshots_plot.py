import pickle
import numpy as np
import pylab
import do_measurements as dm
from plotters import rgb_alpha
from matplotlib import rc
rc('text', usetex=True)

bandcolor = rgb_alpha((0,255,255), 1.)

f = open('pop_mstar_model.dat', 'r')
ms_pop = pickle.load(f)
f.close()

ngal = 100

msgals = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
vdgals = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

f = open('pop_vdisp_model.dat', 'r')
vd_pop = pickle.load(f)
f.close()

snaps = [199, 100, 30]
markers = ['s', '^', 'o']
labels = ['$z=2$', '$z=1$', '$z=0.3$']
colors = ['b', 'g', 'r']
nsnap = len(snaps)

fitdir = '/gdrive/projects/SL2S_hierarch/'

f = open(fitdir+'evol_nfw_alpha_msps_sigma.dat', 'r')
all_chain = pickle.load(f)
f.close()

Nsamp = 1000
indices = np.arange(20000, 100000)
samp = np.random.choice(indices, Nsamp)
astar_samp = all_chain['amstar'][samp]
asigma_samp = all_chain['asigma'][samp]
b_samp = all_chain['calpha'][samp]

xlimm = (10.3,12.5)
xlimv = (2.15,2.75)
ylim = (-0.15, 0.45)

xsm = np.linspace(xlimm[0], xlimm[1], 51)
#ysm = ms_pop.vdisp_coeff[30, 0] + ms_pop.vdisp_coeff[30, 1]*(xsm - 11.)
ysm = all_chain['vdisp_mu'].mean() + all_chain['vdisp_mdep'].mean()*(xsm - 11.5)

# makes observational band
m84 = 0.*xsm
m16 = 0.*xsm
for i in range(0,51):
    aimf_samp = b_samp + astar_samp*(xsm[i] - 11.5) + asigma_samp*(ysm[i] - 2.4)
    m84[i] = np.percentile(aimf_samp, 84.)
    m16[i] = np.percentile(aimf_samp, 16.)

# makes observational band

xsv = np.linspace(xlimv[0], xlimv[1], 51)
#ysv = 11. + (xsv - ms_pop.vdisp_coeff[30, 0])/ms_pop.vdisp_coeff[30, 1]
ysv = 11.5 + (xsv - all_chain['vdisp_mu'].mean())/all_chain['vdisp_mdep'].mean()

v84 = 0.*xsv
v16 = 0.*xsv
for i in range(0,51):
    aimf_samp = b_samp + astar_samp*(ysv[i] - 11.5) + asigma_samp*(xsv[i] - 2.4)
    v84[i] = np.percentile(aimf_samp, 84.)
    v16[i] = np.percentile(aimf_samp, 16.)

fig = pylab.figure()
ax1 = fig.add_subplot(2, 2, 1)
pylab.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99, hspace=0., wspace=0.)
pylab.fill_between(xsm, m16, m84, color=bandcolor)

xs = np.linspace(xlimm[0], xlimm[1])

mfitpars = []
tablines = []
for i in range(nsnap):
    mchab_ms = np.log10(ms_pop.mstar_salp[:, snaps[i]])
    aimf_ms = np.log10(ms_pop.aimf[:, snaps[i]])
    pylab.scatter(mchab_ms, aimf_ms, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_mstar_only(mchab_ms, aimf_ms)
    print i, par, scat
    mfitpars.append((par[0], par[1], scat))
    pylab.plot(xs, par[0] + par[1]*(xs - 11.5), color=colors[i], linestyle='--')

for gal in msgals:
    pylab.plot(np.log10(ms_pop.mstar_salp[gal, 30:]), np.log10(ms_pop.aimf[gal, 30:]), color='k')

pylab.text(11.8, 0., '$M_*$ model', fontsize=14)
#pylab.xlabel('$\log{M_*}$', fontsize=14)
pylab.xlim(xlimm)
pylab.ylim(ylim)

pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.legend(scatterpoints=1, fontsize=14, loc='upper left')
#pylab.xticks(fontsize=14)
pylab.yticks(fontsize=14)
yticks = ax1.yaxis.get_major_ticks()
#yticks[0].label1.set_visible(False)
#yticks[1].label1.set_visible(False)
#yticks[-2].label1.set_visible(False)
pylab.tick_params(axis='x', labelbottom='off')

xs = np.linspace(xlimv[0], xlimv[1])

vfitpars = []
pylab.subplot(2,2,2)
pylab.fill_between(xsv, v16, v84, color=bandcolor)
for i in range(nsnap):
    vdisp_ms = np.log10(ms_pop.veldisp[:, snaps[i]])
    aimf_ms = np.log10(ms_pop.aimf[:, snaps[i]])
    pylab.scatter(vdisp_ms, aimf_ms, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_sigma_only(vdisp_ms, aimf_ms)
    print scat
    vfitpars.append((par[0], par[1], scat))
    pylab.plot(xs, par[0] + par[1]*(xs - 2.4), color=colors[i], linestyle='--')

for i in range(nsnap):
    mpar = mfitpars[i]
    vpar = vfitpars[i]
    tablines.append("''$M_*$ model (%s)'' & (%3.2f, %3.2f, %3.2f) & (%3.2f, %3.2f, %3.2f)\\\\\n"%(labels[i], mpar[1], mpar[0], mpar[2], vpar[1], vpar[0], vpar[2]))

for gal in msgals:
    pylab.plot(np.log10(ms_pop.veldisp[gal, 30:]), np.log10(ms_pop.aimf[gal, 30:]), color='k')

pylab.xlim(xlimv[0], xlimv[1])
pylab.ylim(ylim)

#pylab.xlabel('$\log{\sigma}$', fontsize=14)
#pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.xticks(fontsize=14)
#pylab.yticks(fontsize=16)
pylab.text(2.5, 0., '$M_*$ model', fontsize=14)
pylab.tick_params(axis='y', labelleft='off')
pylab.tick_params(axis='x', labelbottom='off')


alphas = [0.2, 0.5, 1.0]

ax3 = pylab.subplot(2, 2, 3)

pylab.fill_between(xsm, m16, m84, color=bandcolor)

xs = np.linspace(xlimm[0], xlimm[1])

mfitpars = []
for i in range(nsnap):
    mchab_vd = np.log10(vd_pop.mstar_salp[:, snaps[i]])
    aimf_vd = np.log10(vd_pop.aimf[:, snaps[i]])
    pylab.scatter(mchab_vd, aimf_vd, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_mstar_only(mchab_vd, aimf_vd)
    mfitpars.append((par[0], par[1], scat))
    pylab.plot(xs, par[0] + par[1]*(xs - 11.5), color=colors[i], linestyle='--')

for gal in vdgals:
    pylab.plot(np.log10(vd_pop.mstar_salp[gal, 30:]), np.log10(vd_pop.aimf[gal, 30:]), color='k')

pylab.xlabel('$\log{M_*^{\mathrm{Salp}}}$', fontsize=14)
pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.xlim(xlimm[0], xlimm[1])
pylab.ylim(ylim)
pylab.xticks(fontsize=14)
pylab.yticks(fontsize=14)
yticks = ax3.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)
yticks[-1].label1.set_visible(False)
pylab.text(11.8, 0., '$\sigma$ model', fontsize=14)

xticks = ax3.xaxis.get_major_ticks()
#xticks[0].label1.set_visible(False)
xticks[-1].label1.set_visible(False)


ax4 = pylab.subplot(2, 2, 4)
pylab.fill_between(xsv, v16, v84, color=bandcolor)

xs = np.linspace(xlimv[0], xlimv[1])

vfitpars = []
for i in range(0, nsnap):
    vdisp_vd = np.log10(vd_pop.veldisp[:, snaps[i]])
    aimf_vd = np.log10(vd_pop.aimf[:, snaps[i]])
    pylab.scatter(vdisp_vd, aimf_vd, color=colors[i], s=30, marker=markers[i], label=labels[i])
    par, scat = dm.fit_sigma_only(vdisp_vd, aimf_vd)
    print scat
    vfitpars.append((par[0], par[1], scat))
    pylab.plot(xs, par[0] + par[1]*(xs - 2.4), color=colors[i], linestyle='--')

for i in range(nsnap):
    mpar = mfitpars[i]
    vpar = vfitpars[i]
    tablines.append("''$\sigma$ model (%s)'' & (%3.2f, %3.2f, %3.2f) & (%3.2f, %3.2f, %3.2f)\\\\\n"%(labels[i], mpar[1], mpar[0], mpar[2], vpar[1], vpar[0], vpar[2]))

f = open('table1.tex', 'w')
f.writelines(tablines)
f.close()

for gal in vdgals:
    pylab.plot(np.log10(vd_pop.veldisp[gal, 30:]), np.log10(vd_pop.aimf[gal, 30:]), color='k')

pylab.xlabel('$\log{\sigma}$', fontsize=14)
#pylab.ylabel('$\log{\\alpha_{\mathrm{IMF}}}$', fontsize=16)
pylab.xlim(xlimv[0], xlimv[1])
pylab.ylim(ylim)
pylab.xticks(fontsize=14)
#pylab.yticks(fontsize=14)
pylab.text(2.5, 0., '$\sigma$ model', fontsize=14)

pylab.tick_params(axis='y', labelleft='off')
xticks = ax4.xaxis.get_major_ticks()
xticks[0].label1.set_visible(False)
xticks[-1].label1.set_visible(False)


pylab.savefig('snapshots.eps')
#pylab.savefig('puremstar_constsigma_snapshots.png')


pylab.show()





