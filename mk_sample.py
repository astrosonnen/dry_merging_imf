import numpy as np
import shmrs
import popevol
import recipes
import pickle

Ngal = 100   # sample size

z_0 = 2.

imf_recipe = 'mstar'
imf_coeff = (0.11, 0.26)
#imf_recipe = 'vdisp'
#imf_coeff = (-0.06, 1.00)
#imf_recipe = 'mhalo'
#imf_coeff = (0.0, 0.3)

vdisp_coeff = (2.48, 0.20)

outname = 'pop_%s_model_%3.2f_%3.2f.dat'%(imf_recipe, imf_coeff[0], imf_coeff[1])

lmhalos = shmrs.generate_halos(100000, z=2., lmhmax=13.5)
lmstars = shmrs.generate_mstar(lmhalos, z=2., scat=0.18)
selection = lmstars > 10.5
indices = np.arange(100000)
ind_sample = np.random.choice(indices[selection], Ngal)

lmhalo_sample = lmhalos[ind_sample]
lmstar_sample = lmstars[ind_sample]
print len(lmstar_sample[lmstar_sample>12])
reff_sample = recipes.generate_reff(lmstar_sample, z_0)
#vdisp_sample = recipes.generate_veldisp_from_mstar(lmstar_sample, z_0)
vdisp_sample = 10.**(vdisp_coeff[0] + vdisp_coeff[1]*(lmstar_sample - 11.) + np.random.normal(0., 0.04, Ngal))

aimf_z2_sample = 0.*lmstar_sample

pop = popevol.population(z_0=z_0, nobj=Ngal, mstar_chab_0=10.**lmstar_sample, mhalo_0=10.**lmhalo_sample, \
                         veldisp_0=vdisp_sample)

pop.evolve(z_low=0., imf_recipe=imf_recipe, imf_coeff=imf_coeff, vdisp_coeff=vdisp_coeff)

f = open(outname, 'w')
pickle.dump(pop, f)
f.close()
