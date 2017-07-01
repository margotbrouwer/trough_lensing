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
thetalow = np.array([5., 10., 15., 20.]) # in arcmin
thetamin, thetamax = np.array([2., 100.])

am_to_rad = np.pi/(60.*180.)

zmin = 0.1
zmax = 0.3
zgama = 0.5

zlims = np.array([zmin, zmax])

# Path to the GAMA fields
path_gamacat = '/data2/brouwer/MergedCatalogues/'
gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'

# Importing the GAMA coordinates
galRA, galDEC, galZ, rmag, rmag_abs = utils.import_gamacat(path_gamacat, gamacatname)

gamamask = (rmag_abs < -19.67)
zmean = np.mean(galZ[gamamask])
print('Mean GAMA redshift:', zmean)

galmask = (rmag_abs < -21.)& (zmin < galZ)&(galZ < zmax) & (rmag < 19.8)
galRA, galDEC, galZ, rmag, rmag_abs = galRA[galmask], galDEC[galmask], galZ[galmask], rmag[galmask], rmag_abs[galmask]

# Calculating the volume of the cone at each redshift bin
cosmo = LambdaCDM(H0=70., Om0=0.25, Ode0=0.75)

Dcmin, Dcmax, Dcgama, Dcmean = [ (cosmo.comoving_distance(z).to('Mpc')).value for z in [zmin, zmax, zgama, zmean] ] # Comoving distance to each bin limit
Dclim = Dcmin + (Dcmax - Dcmin)/2
zlim = z_at_value(cosmo.comoving_distance, Dclim*u.Mpc)

zlims = [zmin, zlim, zmax]
Dclims = [Dcmin, Dclim, Dcmax]

print()
print('Z(min,lim,max):', zlims)
print('Dc(min,lim,max):', Dclims, 'Mpc')
print('Dc(gama_max):', Dcgama, 'Mpc')
print('L(low,high):', np.diff(Dclims), 'Mpc')

# Equal volume
#tanthetahigh = np.tan(thetalow*am_to_rad) * np.sqrt((Dclims[1]**3. - Dclims[0]**3.)/(Dclims[2]**3. - Dclims[1]**3.))

# Equal radius at Dlim/Dmax
#tanthetahigh = np.tan(thetalow*am_to_rad) * (Dclims[1] / Dclims[2])

# Equal radius at the center of the sample
lowmask = (galZ < zlims[1])
highmask = (zlims[1] < galZ)

Nlow, Nhigh = [sum(lowmask), sum(highmask)]

Zlow, Zhigh = np.array([ np.mean(galZ[mask]) for mask in [lowmask, highmask] ])
Dclow, Dchigh = np.array([ np.mean((cosmo.comoving_distance(galZ[mask]).to('Mpc')).value) for mask in [lowmask, highmask] ])
Dalow, Dahigh, Dagama, Damean = np.array([Dclow/(1+Zlow), Dchigh/(1+Zhigh), Dcgama/(1+zgama), Dcmean/(1+zmean)])
Vclow, Vchigh = np.array([ np.mean((cosmo.comoving_volume(galZ[mask]).to('Mpc^3')).value) for mask in [lowmask, highmask] ])

thetahigh = thetalow * (Dalow / Dahigh)

print('mean Z(low,high):', [Zlow, Zhigh])
print('mean Dc(low,high):', [Dclow, Dchigh], 'Mpc')
print('mean Da(low,high):', [Dalow, Dahigh], 'Mpc')
print()
print('theta(low):', thetalow, 'arcmin')
print('theta(high):', thetahigh, 'arcmin')
print('R(low):', (thetalow*am_to_rad) * Dalow, 'Mpc')
print('R(high):', (thetahigh*am_to_rad) * Dahigh, 'Mpc')
print('R(min,max):', (thetamin*am_to_rad) * Damean, (thetamax*am_to_rad) * Damean, 'Mpc')

V1low = 1./3.*np.pi * Dclims[0]**3. * (thetalow*am_to_rad)**2.
V2low = 1./3.*np.pi * Dclims[1]**3. * (thetalow*am_to_rad)**2.
V2high = 1./3.*np.pi * Dclims[1]**3. * (thetahigh*am_to_rad)**2.
V3high = 1./3.*np.pi * Dclims[2]**3. * (thetahigh*am_to_rad)**2.

Vlow = V2low - V1low
Vhigh = V3high - V2high

print('Cone volume(low):', Vlow, 'Mpc^3')
print('Cone volume(high):', Vhigh, 'Mpc^3')
print()
print('Galaxy number(low,high):', Nlow, Nhigh)
print('Volume(low,high):', Vclow, Vchigh, 'Mpc^3')
print('Density(low):', Nlow/Vclow, Nhigh/Vchigh, 'Mpc^-3')

quit()

Mpc_am = (cosmo.kpc_comoving_per_arcmin(zlims).to('Mpc/arcmin')).value # Comoving distance per arcmin at each bin limit
areabins = np.pi * (theta * Mpc_am)**2. # Comoving area of the circle at each bin limit
#covolbins = cosmo.comoving_volume(zbins).to('kpc3').value

covolbins = 1./3. * areabins * Dclims # Comoving cone volume below each bin limit
covolbins_high = covolbins[-1] - covolbins # Comoving cone volume above each bin limit

density = Ngals/covolbins # Density below the redshift limit
density_high = Ngals_high/covolbins_high # Density above the redshift limit


