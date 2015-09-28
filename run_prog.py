import numpy as np
import pylab
import shmrs
import galaxies
import recipes
from matplotlib import rc
rc('text', usetex=True)

Ngal = 10   # sample size

# generates halo masses, then assigns stellar masses, then makes a cut in stellar mass and draws a
# sample of Ngal galaxies

Nhalo = 10000
lmhmin = 12.
lmhmax = 14.

z_0 = 0.

lmhalos = shmrs.generate_halos(Nhalo, lmhmin, lmhmax, z=z_0)

lmstars = shmrs.generate_mstar(lmhalos, z_0)

selection = lmstars > 11.

lmhalo_cut = lmhalos[selection]
lmstar_cut = lmstars[selection]

indices = np.arange(len(lmhalo_cut))

sample = np.random.choice(indices, Ngal)

lmstar_sample = lmstar_cut[sample]
lmhalo_sample = lmhalo_cut[sample]

reff_sample = recipes.generate_reff(lmstar_sample)
veldisp_sample = recipes.generate_veldisp_from_fp(lmstar_sample, reff_sample)

gal0 = galaxies.ETG(z_0=0., mstar_chab_0=1e11, mhalo_0=4.4846e12, re_0=5., sigma_0=250.)
gal0.get_sf_history()
gal0.evolve(dz = 0.001, z_up = gal0.z_form)

centrals = []
for i in range(0, Ngal):
    central = galaxies.ETG(z_0=z_0, mstar_chab_0=10.**lmstar_sample[i], mhalo_0=10.**lmhalo_sample[i], \
                           re_0=reff_sample[i], sigma_0=veldisp_sample[i])

    central.get_sf_history()

    central.evolve(z_up = central.z_form)

    centrals.append(central)
