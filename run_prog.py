import numpy as np
import pylab
import shmrs
import galaxies
from matplotlib import rc
rc('text',usetex=True)

Ngal = 10
zref = 1.8

# generates halo masses, then assigns stellar masses, then makes a cut in stellar mass and draws a
# sample of Ngal galaxies

Nhalo = 10000
lmhmin = 12.
lmhmax = 14.

z_0 = 0.

lmhalos = shmrs.generate_halos(Nhalo, lmhmin, lmhmax, z=z_0)

lmstars = shmrs.mstarfunc(lmhalos, z_0)

selection = lmstars > 11.

print len(lmhalos[selection])
