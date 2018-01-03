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


h, O_matter, O_lambda = [0.7, 0.25, 0.75]
cosmo = LambdaCDM(H0=h*100, Om0=O_matter, Ode0=O_lambda)
micecor = 5*np.log10(h) # Correction on MICE absmag_r (= -0.7745)


# Select catalog and redshift binning
cat = 'mice'
zmax = 0.6
Nbins = 5


# Defining the circle size and redshift bins
thetalist = np.array([5., 10., 15., 20.]) # in arcmin
if 'mice' in cat:
    thetalow = np.array([thetalist[3]]) # in arcmin
else:
    thetalow = np.array([thetalist[2]]) # in arcmin

thetamin, thetamax = np.array([2., 100.])

am_to_rad = np.pi/(60.*180.)
zmin = 0.1
zgama = 0.5


# Path to the KiDS fields
if cat == 'kids':

    # Path to the KiDS fields
    path_kidscat = '/data2/brouwer/KidsCatalogues'
    kidscatname = '/KiDS_DR3_GAMA-like_Maciek_revised_1905.fits'

    # Importing the KiDS galaxies
    galRA, galDEC, galZB, galZ, galTB, mag_auto, ODDS, umag, gmag, rmag, imag = \
    utils.import_kidscat(path_kidscat, kidscatname)
    rmag_abs = utils.calc_absmag(rmag, galZ, gmag, imag, h, O_matter, O_lambda)
    
    galDc = (cosmo.comoving_distance(galZ).to('Mpc')).value
    gama_rlim = 20.2

# Path to the GAMA fields
if cat == 'gama':

    path_gamacat = '/data2/brouwer/MergedCatalogues/'
    gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'

    # Importing the GAMA coordinates
    galRA, galDEC, galZ, rmag, rmag_abs = \
    utils.import_gamacat(path_gamacat, gamacatname)
    
    galDc = (cosmo.comoving_distance(galZ).to('Mpc')).value
    gama_rlim = 19.8

if cat == 'mice':

    # Path to the Mice field
    path_mockcat = '/data2/brouwer/MergedCatalogues'
    mockcatname = 'mice_gama_highZ_catalog.fits'
    
    # Importing the Mice galaxies
    galRA, galDEC, galZ, galDc, rmag, rmag_abs, e1, e2, galmass = \
    utils.import_mockcat(path_mockcat, mockcatname)
    rmag_abs = rmag_abs + micecor # Correct absolute magnitudes
    
    #gama_rlim = 20.2
    gama_rlim = np.inf

gamamask = (rmag_abs < -19.7) & (rmag < gama_rlim)

zmean = np.mean(galZ[gamamask])
Dcmean = np.mean(galDc[gamamask])
Damean = np.mean((galDc/(1+galZ))[gamamask])
#Damean = (cosmo.angular_diameter_distance(zmean).to('Mpc')).value

print('Mean GAMA redshift:', zmean)

Dcmin, Dcmax, Dcgama = [ (cosmo.comoving_distance(z).to('Mpc')).value for z in [zmin, zmax, zgama] ] # Comoving distance to each bin limit
Dclen = (Dcmax - Dcmin)/Nbins

# Calculating the comoving length for each Z-bin
zlims = np.array([zmin])
Dclims = np.array([Dcmin])

for N in range(Nbins):

    Dclim = Dcmin + (N+1.)*Dclen
    zlim = z_at_value(cosmo.comoving_distance, Dclim*u.Mpc)

    Dclims = np.append(Dclims, Dclim)
    zlims = np.append(zlims, zlim)

print()
print('Z-limits:', zlims)
print('Dc-limits:', Dclims, 'Mpc')
print('Dc(gama_max):', Dcgama, 'Mpc')
print('L:', np.diff(Dclims), 'Mpc')

Nhigh, Zhigh, Dchigh, Dahigh, Vchigh = [np.array([])]*5
for N in range(Nbins):

    print(N+1)
    
    # Masking the GAMA galaxies
    galmask = (rmag_abs < -21.) & (rmag < gama_rlim) & (zlims[N] < galZ) & (galZ < zlims[N+1])
    Z = galZ[galmask]
    
    Nhigh = np.append( Nhigh, sum(galmask) )
    Zhigh = np.append( Zhigh, np.mean(Z) )

    if 'mice' in cat:
        Dc = galDc[galmask]
    else:
        Dc = (cosmo.comoving_distance(Z).to('Mpc')).value
    
    Dchigh = np.append( Dchigh, np.mean(Dc) )
    Dahigh = np.append( Dahigh, np.mean(Dc/(1+Z)) )


# Equal comoving projected radius at all redshifts
Dclow = Dchigh[0]
Alow = 20.*40.

thetahigh = [thetalow * (Dclow / Dchigh[N]) for N in range(Nbins)]
Ahigh = np.array([Alow * (Dclow / Dchigh[N])**2 for N in range(Nbins)])

Rlow = (thetalow*am_to_rad) * Dclow
Rhigh = [(thetahigh[N]*am_to_rad) * Dchigh[N] for N in range(Nbins)]

print('mean Z:', Zhigh)
print('mean Dc:', Dchigh, 'Mpc')
print('mean Da:', Dahigh, 'Mpc')
print()
print('theta(low):', thetalow, 'arcmin')
print('theta(high):', thetahigh, 'arcmin')
print('R(low):', Rlow, 'Mpc')
print('R(high):', Rhigh, 'Mpc')
print('A(high):', Ahigh, 'degree^2')
print('W(high):', Ahigh/20., 'degree')
print()
print('R(min,max):', (thetamin*am_to_rad) * Dcmean, (thetamax*am_to_rad) * Dcmean, 'Mpc')


quit()
for N in range(Nbins):
    Vlow = 1./3.*np.pi * Dclims[0]**3. * (thetalow*am_to_rad)**2.
    Vhigh = 1./3.*np.pi * Dclims[1]**3. * (thetahigh*am_to_rad)**2.

    Vlow = V2low - V1low
    Vhigh = V3high - V2high

    print('Cone volume(low):', Vlow, 'Mpc^3')
    print('Cone volume(high):', Vhigh, 'Mpc^3')
    print()
    print('Galaxy number(low,high):', Nlow, Nhigh)
    print('Volume(low,high):', Vclow, Vchigh, 'Mpc^3')
    print('Density(low):', Nlow/Vclow, Nhigh/Vchigh, 'Mpc^-3')
