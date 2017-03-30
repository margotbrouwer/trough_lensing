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
from astropy.cosmology import LambdaCDM, z_at_value

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

import trough_modules_all as utils

# Defining the circle size and redshift bins
thetalow = 10. # in arcmin
am_to_rad = np.pi/(60.*180.)

zmin = 0.05
zmax = 0.3

zlims = np.array([zmin, zmax])

# Path to the GAMA fields
path_gamacat = '/data2/brouwer/MergedCatalogues/'
gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'

# Importing the GAMA coordinates
galRA, galDEC, galZ, rmag, rmag_abs = utils.import_gamacat(path_gamacat, gamacatname)

galmask = (rmag_abs < -21.)&(galZ < zmax)
galRA, galDEC, galZ, rmag, rmag_abs = galRA[galmask], galDEC[galmask], galZ[galmask], rmag[galmask], rmag_abs[galmask]

# Calculating the volume of the cone at each redshift bin
cosmo = LambdaCDM(H0=70., Om0=0.315, Ode0=0.685)

Dcmin, Dcmax = [ (cosmo.comoving_distance(z).to('Mpc')).value for z in [zmin, zmax] ] # Comoving distance to each bin limit
Dclim = Dcmin + (Dcmax - Dcmin)/2
zlim = z_at_value(cosmo.comoving_distance, Dclim*u.Mpc)

zlims = [zmin, zlim, zmax]
Dclims = [Dcmin, Dclim, Dcmax]

print()
print('Z(min,lim,max):', zlims)
print('Dc(min,lim,max):', Dclims, 'Mpc')

tanthetahigh = np.tan(thetalow*am_to_rad) * np.sqrt((Dclims[1]**3. - Dclims[0]**3.)/(Dclims[2]**3. - Dclims[1]**3.))
thetahigh = np.arctan(tanthetahigh)/am_to_rad
print('theta(low,high):', thetalow, thetahigh, 'arcmin')

V1low = 1./3.*np.pi * Dclims[0]**3. * np.tan(thetalow*am_to_rad)**2.
V2low = 1./3.*np.pi * Dclims[1]**3. * np.tan(thetalow*am_to_rad)**2.
V2high = 1./3.*np.pi * Dclims[1]**3. * np.tan(thetahigh*am_to_rad)**2.
V3high = 1./3.*np.pi * Dclims[2]**3. * np.tan(thetahigh*am_to_rad)**2.

Vlow = V2low - V1low
Vhigh = V3high - V2high

print('Volume(low,high):', Vlow, Vhigh, 'Mpc^3')

quit()

Mpc_am = (cosmo.kpc_comoving_per_arcmin(zlims).to('Mpc/arcmin')).value # Comoving distance per arcmin at each bin limit
areabins = np.pi * (theta * Mpc_am)**2. # Comoving area of the circle at each bin limit
#covolbins = cosmo.comoving_volume(zbins).to('kpc3').value

covolbins = 1./3. * areabins * Dclims # Comoving cone volume below each bin limit
covolbins_high = covolbins[-1] - covolbins # Comoving cone volume above each bin limit

density = Ngals/covolbins # Density below the redshift limit
density_high = Ngals_high/covolbins_high # Density above the redshift limit


