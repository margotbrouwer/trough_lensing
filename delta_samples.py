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

# Radii theta of circular regions (in deg)
thetalist = np.array([5., 10., 15., 20.])
#thetalist = np.array([5.])/60.

Nbins = 20

# Path to the trough catalogues
path_troughcat = '/data2/brouwer/MergedCatalogues/'
troughfilename = 'trough_catalogues/trough_catalog_gama_absmag.fits'

# Full directory & name of the corresponding KiDS catalogue
troughcatfile = '%s/%s'%(path_troughcat, troughfilename)
troughcat = pyfits.open(troughcatfile, memmap=True)[1].data

deltalimits = np.zeros([len(thetalist), Nbins+1])
deltacenters = np.zeros([len(thetalist), Nbins])

for theta in range(len(thetalist)):
    
    deltalist = troughcat['delta%i'%thetalist[theta]]
    deltalist = deltalist[np.isfinite(deltalist)]
    
    sorted_delta = np.sort(deltalist)
    limits = np.linspace(0, len(deltalist)-1, Nbins+1)
    limits = np.array(limits, dtype=int)
    
    deltalimits[theta] = sorted_delta[limits]
    deltacenters[theta] = np.array([np.mean(sorted_delta[limits[i]:limits[i+1]]) for i in range(Nbins)])

print(deltalimits)

filename = '/data2/brouwer/shearprofile/KiDS-GAMA/brouwer/configs_margot/deltalimits.txt'
field_header = 'Limits for stacking troughs as a function of delta (with equal number of troughs per bin)'
np.savetxt(filename, deltalimits.T, delimiter='    ', header = field_header)
print('Written:', filename)

filename = '/data2/brouwer/shearprofile/KiDS-GAMA/brouwer/configs_margot/deltacenters.txt'
field_header = 'Centers of the delta-bins (with equal number of troughs per bin)'
np.savetxt(filename, deltacenters.T, delimiter='    ', header = field_header)
print('Written:', filename)
