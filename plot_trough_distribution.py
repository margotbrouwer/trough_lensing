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

# Make use of TeX
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})

colors = ['red', 'orange', 'cyan', 'blue']

# Radii theta of circular regions (in deg)
thetalist = np.array([5., 10., 15., 20.])/60.

path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogues'
troughcatname = 'trough_catalog_gama_absmag_0.04deg_complex.fits'

# Full directory & name of the trough catalogue
troughcatfile = '%s/%s'%(path_troughcat, troughcatname)
troughcat = pyfits.open(troughcatfile, memmap=True)[1].data

# List of the observables of all sources in the KiDS catalogue
RA = troughcat['RA']
DEC = troughcat['DEC']

Ngaltheta = [troughcat['Ngaltheta%g'%(theta*60)] for theta in thetalist]
rhotheta = [troughcat['rhotheta%g'%(theta*60)] for theta in thetalist]
delta = [troughcat['delta%g'%(theta*60)] for theta in thetalist]
Pmasktheta = [troughcat['Pmasktheta%g'%(theta*60)] for theta in thetalist]

fig = plt.figure(figsize=(5,4))

for t in range(len(thetalist)):
    
    Ngal = (rhotheta[t])[Pmasktheta[t]>0.8]
    n, bins, patches = plt.hist(Ngal, bins=30, histtype='step', range=[0., 0.8], \
                label=r"$\theta_{\rm A} = %g'$"%(thetalist[t]*60.), alpha=1., color=colors[t])
    plt.axvline( x=np.mean((rhotheta[t])[Pmasktheta[t]>0.8]), ls = '--', color=colors[t])
    
    plt.xlabel(r'Galaxy number density $\rho(\theta_{\rm A})$ (${\rm arcmin}^{-2}$)', fontsize=14)
    plt.ylabel(r'Number of apertures', fontsize=14)

    #plt.yscale('log')
    #plt.yscale('log')
    plt.legend()

plotfilename = '/data2/brouwer/shearprofile/trough_results_Apr/Plots/trough_density_distribution'

plt.tight_layout()

# Save plot
for ext in ['pdf']:
    plotname = '%s.%s'%(plotfilename, ext)
    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written: ESD profile plot:', plotname)

plt.show()
plt.clf
