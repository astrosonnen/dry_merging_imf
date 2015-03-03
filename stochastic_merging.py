#This code generates ensembles of halos drawn from a halo mass function, then merges them stochastically.

import numpy as np

def generate_ensemble(mh_min=1e10,mh_max=10.**13.5,alpha=2.): #draws halos from a power-law halo mass function with exponent alpha and minimum mass mh_min. The halos are drawn one-by-one until the sum of all halo masses exceeds mh_max.

    ensemble = 0.

    return ensemble

def evolve_ensemble(ensemble,zf=0.,zi=2.,dz=0.001): #evolves the ensemble by merging its halos at each timestep.

    zs = np.arange(zf,zi,dz)
    Ns = len(zs)

    for i in range(0,Ns):
        pass

    new_ensemble = ensemble

    return new_ensemble



