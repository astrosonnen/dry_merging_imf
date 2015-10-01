import numpy as np
import pylab
import shmrs
import galaxies
import recipes
import sys
import pickle


Ngal = 10   # sample size

z_0 = 2.

if len(sys.argv) > 1:
    outname = sys.argv[1]
else:
    outname = 'output.dat'

lmhalos = shmrs.generate_halos(10000, z=2.)
lmstars = shmrs.generate_mstar(lmhalos, z=2., scat=0.18)
selection = lmstars > 10.8
indices = np.arange(10000)
ind_sample = np.random.choice(indices[selection], Ngal)

lmhalo_sample = lmhalos[ind_sample]
lmstar_sample = lmstars[ind_sample]

centrals = []
for i in range(0, Ngal):
    central = galaxies.ETG(z_0=z_0, mstar_chab_0=10.**lmstar_sample[i], mhalo_0=10.**lmhalo_sample[i], \
                           re_0=5., sigma_0=200.)

    central.z_form = 2.
    central.evolve(z_low=0., z_up = z_0, imf_recipe='mstar', imf_coeff=(3., 0.2))

    centrals.append(central)

f = open(outname, 'w')
pickle.dump(centrals, f)
f.close()
