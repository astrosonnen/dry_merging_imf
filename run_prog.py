import numpy as np
import pylab
import shmrs
import galaxies
from matplotlib import rc
rc('text', usetex=True)

Ngal = 10   # sample size
zref = 1.8

mstar_scat = 0.1    # scatter in logmstar around the stellar-to-halo mass relation

# generates halo masses, then assigns stellar masses, then makes a cut in stellar mass and draws a
# sample of Ngal galaxies

Nhalo = 10000
lmhmin = 12.
lmhmax = 14.

z_0 = 0.

lmhalos = shmrs.generate_halos(Nhalo, lmhmin, lmhmax, z=z_0)

lmstars = shmrs.mstarfunc(lmhalos, z_0) + np.random.normal(0., mstar_scat, Nhalo)

selection = lmstars > 11.

lmhalos_cut = lmhalos[selection]
lmstars_cut = lmstars[selection]

indices = np.arange(len(lmhalos_cut))

sample = np.random.choice(indices, Ngal)

lmstars_sample = lmstars_cut[sample]
lmhalos_sample = lmhalos_cut[sample]




