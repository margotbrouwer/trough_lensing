#!/usr/bin/python

"Module to determine the local overdensity delta_r within a sphere of radius r."

# Import the necessary libraries
import astropy.io.fits as pyfits
import gc
import numpy as np
import sys
import os
import time
from glob import glob

from astropy import constants as const, units as u
from astropy.coordinates import SkyCoord
from collections import Counter
from astropy.cosmology import LambdaCDM

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

import trough_modules_all as utils

theta = 5. # in arcmin

# Names of the GAMA fields
fieldnames = ['G9', 'G12', 'G15']

# Boundaries of the GAMA fields
coordsG9 = [[129., 141.], [-2.,3.]]
coordsG12 = [[174., 186.], [-3.,2.]]
coordsG15 = [[211.5, 223.5], [-2.,3.]]
fieldboundaries = np.array([coordsG9,coordsG12,coordsG15])

# Path to the KiDS fields
path_gamacat = '/data2/brouwer/MergedCatalogues/'
gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'

# Importing the GAMA coordinates
galRA, galDEC, galZ, rmag, rmag_abs = utils.import_gamacat(path_gamacat, gamacatname)

galmask = (rmag_abs < -19.7)&(rmag <= 19.8)
galRA, galDEC, galZ, rmag, rmag_abs = galRA[galmask], galDEC[galmask], galZ[galmask], rmag[galmask], rmag_abs[galmask]

dz = 0.01
zmin = 0.1
zmax = 0.5

zbins = np.arange(zmin, zmax, dz)
zcenters = (zbins + dz/2.)
zbins = np.append(zbins, zmax)
print('Z-bins:', zbins)



cosmo = LambdaCDM(H0=70., Om0=0.315, Ode0=0.685)
Dcbins = (cosmo.comoving_distance(zbins).to('kpc')).value
covollist = cosmo.comoving_volume(zbins).to('kpc3').value

Ngals = np.zeros(len(zcenters))

for z in range(len(zcenters)):
    
    zmask = (zbins[z] < galZ) & (galZ < zbins[z+1])
    print('Zbin: %g, galaxies: %g'%(zcenters[z],  np.sum(zmask)))
    Ngals[z] = np.sum(zmask)

kpc_am = (cosmo.kpc_comoving_per_arcmin(zcenters)).value
area = np.array([np.pi * (theta*kpc_am[z])**2. for z in range(len(zcenters))])
density = Ngals/area

    
plt.plot(zcenters, density)
plt.ylabel('Density')
plt.show()

plt.plot(zcenters, Ngals)
plt.ylabel('Ngals')
plt.show()

plt.plot(zcenters, area)
plt.ylabel('Area')
plt.show()
