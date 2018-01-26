#!/usr/bin/python

# Import the necessary libraries
import astropy.io.fits as pyfits
import gc
import numpy as np
import sys
import os
import time
from glob import glob

from astropy.cosmology import LambdaCDM

h, O_matter, O_lambda = [0.7, 0.29, 0.71]
cosmo = LambdaCDM(H0=h*100, Om0=O_matter, Ode0=O_lambda)
troughZ = np.array([0.1539674, 0.24719192, 0.33112174, 0.42836386])
am_to_rad = np.pi/(60.*180.)

data = np.load('/data2/brouwer/shearprofile/trough_results_final/slics_mockZ_nomask/Redshift_bins_covariance.npy')
data_x = np.array([np.array(data[0][x+1][0]) for x in range(len(data[0]))])

#data_R = np.array([ data_x[x]*am_to_rad * (cosmo.angular_diameter_distance(troughZ[x]).to('Mpc')).value \
#    for x in range(len(data_x))])

data_R = np.array([ data_x[x]*am_to_rad * (cosmo.comoving_distance(troughZ[x]).to('Mpc')).value \
    for x in range(len(data_x))])

covariance_tot = np.array((data[2][1]).values())

print(covariance_tot[0])
print()
print(data[2][1][0])

#print(data_cov)
#print(data_y)

print(np.shape(data))

# 0: R-bins
# 1: ESD profiles, 1-4: Redshifts, 0-9: Percentiles
# 2: Cov-matrices, 1-4: Redshifts, 0-9: Percentiles
# 3: Cor-matrices, 1-4: Redshifts, 0-9: Percentiles
