import numpy as np
import shmrs
import popevol_noimfevol
import recipes
import pickle

Ngal = 100 # sample size

z_0 = 2.
ximin = 0.03

#modelname = 'mstar'
#imf_coeff = (0.21, 0.26, 0.)

#modelname = 'vdisp_n1000'
modelname = 'noimfevol'
imf_coeff = (0.0, 0., 1.20)

vdisp_coeff = (2.48, 0.20)

outname = 'pop_%s_model.dat'%(modelname)

lmhalos = shmrs.generate_halos(100000, z=2., lmhmax=13.5)
lmstars = shmrs.generate_mstar(lmhalos, z=2., scat=0.18)
selection = lmstars > 10.3
indices = np.arange(100000)
ind_sample = np.random.choice(indices[selection], Ngal)

lmhalo_sample = lmhalos[ind_sample]
lmstar_sample = lmstars[ind_sample]
print len(lmstar_sample[lmstar_sample>12])
reff_sample = recipes.generate_reff(lmstar_sample, z_0)
#vdisp_sample = recipes.generate_veldisp_from_mstar(lmstar_sample, z_0)
vdisp_sample = 10.**(vdisp_coeff[0] + vdisp_coeff[1]*(lmstar_sample - 11.) + np.random.normal(0., 0.04, Ngal))

aimf_z2_sample = 0.*lmstar_sample

pop = popevol_noimfevol.population(z_0=z_0, nobj=Ngal, mstar_chab_0=10.**lmstar_sample, mhalo_0=10.**lmhalo_sample, \
                         veldisp_0=vdisp_sample)

pop.evolve(z_low=0., imf_recipe='mstar-vdisp', imf_coeff=imf_coeff, vdisp_coeff=vdisp_coeff, ximin=ximin)

f = open(outname, 'w')
pickle.dump(pop, f)
f.close()
