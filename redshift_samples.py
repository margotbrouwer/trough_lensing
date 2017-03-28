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

# Defining the circle size and redshift bins
theta = 5. # in arcmin

dz = 1e-5
zmin = 0.
zmax = 0.3

zbins = np.arange(zmin, zmax, dz)
zcenters = (zbins + dz/2.)
zbins = np.append(zbins, zmax)
print('Z-bins:', zbins)


# Path to the GAMA fields
path_gamacat = '/data2/brouwer/MergedCatalogues/'
gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'

# Importing the GAMA coordinates
galRA, galDEC, galZ, rmag, rmag_abs = utils.import_gamacat(path_gamacat, gamacatname)

galmask = (rmag_abs < -21.)&(galZ < zmax)
galRA, galDEC, galZ, rmag, rmag_abs = galRA[galmask], galDEC[galmask], galZ[galmask], rmag[galmask], rmag_abs[galmask]


# Calculating the number of galaxies...
Ngals = np.array([np.sum(galZ<zbins[z]) for z in range(len(zbins))]) # below each redshift bin
Ngals_high = len(galZ)-Ngals # above each redshift bin

# Calculating the volume of the cone at each redshift bin
cosmo = LambdaCDM(H0=70., Om0=0.315, Ode0=0.685)

Dcbins = (cosmo.comoving_distance(zbins).to('Mpc')).value # Comoving distance to each bin limit
Mpc_am = (cosmo.kpc_comoving_per_arcmin(zbins).to('Mpc/arcmin')).value # Comoving distance per arcmin at each bin limit
areabins = np.pi * (theta * Mpc_am)**2. # Comoving area of the circle at each bin limit
#covolbins = cosmo.comoving_volume(zbins).to('kpc3').value

covolbins = 1./3. * areabins * Dcbins # Comoving cone volume below each bin limit
covolbins_high = covolbins[-1] - covolbins # Comoving cone volume above each bin limit

density = Ngals/covolbins # Density below the redshift limit
density_high = Ngals_high/covolbins_high # Density above the redshift limit


densmask = (0.2<zbins)&(zbins<0.25)
densplus = np.sum(zbins<0.2)

plt.plot(zbins, density, label='Low-Z sample')
plt.plot(zbins, density_high, label='High-Z sample')
plt.ylabel('Density')
#plt.show()
plt.close()

plt.plot(zbins, Ngals, label='Low-Z sample')
plt.plot(zbins, Ngals_high, label='High-Z sample')
plt.ylabel('Number of galaxies')
#plt.show()
plt.close()

plt.plot(zbins, covolbins, label='Low-Z sample')
plt.plot(zbins, covolbins_high, label='High-Z sample')
plt.ylabel('Volume')
#plt.show()
plt.close()

plt.hist2d(galZ, rmag_abs)
plt.colorbar()
plt.invert_yaxis()
plt.show()

idxdensity = densplus + (np.abs(density[densmask] - density_high[densmask])).argmin()
idxNgals = (np.abs(Ngals - Ngals_high)).argmin()
idxcovol = (np.abs(covolbins - covolbins_high)).argmin()

print()
print('Density (#/Mpc^3): number=%g, Zlim=%g, lowsamp=%g, highsamp=%g'%(idxdensity, zbins[idxdensity], density[idxdensity], density_high[idxdensity]))
print('Ngals (#): number=%g, Zlim=%g, lowsamp=%g, highsamp=%g'%(idxNgals, zbins[idxNgals], Ngals[idxNgals], Ngals_high[idxNgals]))
print('Volume (Mpc^3): number=%g, Zlim=%g, lowsamp=%g, highsamp=%g'%(idxcovol, zbins[idxcovol], covolbins[idxcovol], covolbins_high[idxcovol]))
