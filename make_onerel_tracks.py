import numpy as np
import pylab
import shmrs
import galaxies
import recipes
import pickle
import do_measurements

fsize = 20

z_0 = 2.

nit = 3

mstar_filename = 'pop_mstar_model_0.11_0.28.dat'
vdisp_filename = 'pop_vdisp_model_-0.18_1.50.dat'

f = open(mstar_filename, 'r')
mstar_pop = pickle.load(f)
f.close()

f = open(vdisp_filename, 'r')
vdisp_pop = pickle.load(f)
f.close()

ngal = 100
nstep = len(mstar_pop.z)

# fits for the imf coefficients

def get_mstar_imf_coeffs(pop):
    imf_coeffs = []
    for j in range(nstep):
	mstars = np.log10(pop.mstar_chab[:, j])
	aimfs = np.log10(pop.aimf[:, j])

        coeff_imf = do_measurements.fit_mstar_only(mstars, aimfs)[0]
        imf_coeff = (coeff_imf[1], coeff_imf[0])
        imf_coeffs.append(imf_coeff)
    return imf_coeffs


def get_vdisp_imf_coeffs(pop):
    imf_coeffs = []
    for j in range(nstep):
	sigmas = np.log10(pop.veldisp[:, j])
	aimfs = np.log10(pop.aimf[:, j])

        coeff_imf = do_measurements.fit_sigma_only(sigmas, aimfs)[0]
        imf_coeff = (coeff_imf[1], coeff_imf[0])
        imf_coeffs.append(imf_coeff)
    return imf_coeffs


mstar_model_mstar_coeffs = np.array(get_mstar_imf_coeffs(mstar_pop))
vdisp_model_mstar_coeffs = np.array(get_mstar_imf_coeffs(vdisp_pop))
mstar_model_vdisp_coeffs = np.array(get_vdisp_imf_coeffs(mstar_pop))
vdisp_model_vdisp_coeffs = np.array(get_vdisp_imf_coeffs(vdisp_pop))

pylab.subplot(2, 1, 1)
pylab.plot(mstar_pop.z, mstar_model_mstar_coeffs[:, 0], label='$M_*$ model')
pylab.plot(vdisp_pop.z, vdisp_model_mstar_coeffs[:, 0], label='$\sigma$ model')
pylab.ylabel('$a_*$', fontsize=fsize)
pylab.legend(loc = 'upper left', fontsize=fsize)
pylab.title('Eq. 3 fit', fontsize=fsize)

pylab.subplot(2, 1, 2)
pylab.plot(mstar_pop.z, mstar_model_mstar_coeffs[:, 1])
pylab.plot(vdisp_pop.z, vdisp_model_mstar_coeffs[:, 1])
pylab.ylabel('$b_*$', fontsize=fsize)
pylab.xlabel('$z$', fontsize=fsize)
pylab.savefig('eq3_fit_plot.png')
pylab.show()


pylab.subplot(2, 1, 1)
pylab.plot(mstar_pop.z, mstar_model_vdisp_coeffs[:, 0], label='$M_*$ model')
pylab.plot(vdisp_pop.z, vdisp_model_vdisp_coeffs[:, 0], label='$\sigma$ model')
pylab.ylabel('$a_\sigma$', fontsize=fsize)
pylab.legend(loc = 'upper left', fontsize=fsize)
pylab.title('Eq. 4 fit', fontsize=fsize)

pylab.subplot(2, 1, 2)
pylab.plot(mstar_pop.z, mstar_model_vdisp_coeffs[:, 1])
pylab.plot(vdisp_pop.z, vdisp_model_vdisp_coeffs[:, 1])
pylab.ylabel('$b_\sigma$', fontsize=fsize)
pylab.xlabel('$z$', fontsize=fsize)
pylab.savefig('eq4_fit_plot.png')
pylab.show()

