import numpy as np
import pylab
import shmrs
import galaxies
import recipes
import pickle
import do_measurements

z_0 = 2.

nit = 3

imf_recipe = 'vdisp'
filename = '%s_dep_imf_coeff2.30.dat'%imf_recipe

f = open(filename, 'r')
objects = pickle.load(f)
f.close()

ngal = len(objects)
nstep = len(objects[0].z)

# fits for the imf coefficients

def get_imf_coeffs(objlist):
    imf_coeffs = []
    for j in range(nstep):
        mstars = np.zeros(ngal)
	sigmas = np.zeros(ngal)
        aimfs = np.zeros(ngal)
        for i in range(ngal):
            mstars[i] = np.log10(objlist[i].mstar_chab[j])
            sigmas[i] = np.log10(objlist[i].veldisp[j])
            aimfs[i] = np.log10(objlist[i].aimf[j])

        coeff_imf = do_measurements.fit_sigma_only(sigmas, aimfs)[0]
        imf_coeff = (coeff_imf[1], coeff_imf[0])
        imf_coeffs.append(imf_coeff)
    return imf_coeffs


imf_coeffs = get_imf_coeffs(objects)

for j in range(nit):

    centrals = []
    for i in range(ngal):
        print j, i
        central = galaxies.ETG(z_0=z_0, mstar_chab_0=objects[i].mstar_chab_0, mhalo_0=objects[i].mhalo_0, \
                               re_0=objects[i].re_0, sigma_0=objects[i].sigma_0)

        central.z_form = z_0
        central.evolve(z_low=0., z_up = z_0, imf_recipe=imf_recipe, imf_coeffs=imf_coeffs)

        centrals.append(central)

    imf_coeffs = get_imf_coeffs(centrals)

    f = open(filename.replace('.dat', '_it%d.dat'%(j+1)), 'w')
    pickle.dump(centrals, f)
    f.close()

