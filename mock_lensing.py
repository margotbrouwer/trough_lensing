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

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

import trough_modules_all as utils
import treecorr


# Select the galaxy catalogue for trough selection (kids, gama or mice)
cat = 'mice'

# Name of the pre-defined galaxy selection
selection = 'all'

# Select mask type (nomask or complex)
masktype = 'nomask'

outputname = 'treecorr_test.txt'

# Import trough catalog
path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'
troughcatname = 'trough_catalog_%s_%s_%s.fits'%(cat, selection, masktype)

theta = 5
Rmin = 2 # in arcmin
Rmax = 100 # in arcmin

paramnames = np.array(['Ptheta%g'%theta, 'Pmasktheta%g'%theta ])
maskvals = np.array([[0., 0.2], [0.8, np.inf]])

troughRA, troughDEC, troughZ, paramlists = utils.import_troughcat(path_troughcat, troughcatname, paramnames)


# Selecting the troughs
masklist = np.ones(len(troughRA))
for m in range(len(maskvals)):
    troughmask = (maskvals[m,0] <= paramlists[m]) & (paramlists[m] < maskvals[m,1])
    masklist[np.logical_not(troughmask)] = 0

troughmask = (masklist == 1)
RA, DEC, Z = troughRA[troughmask], troughDEC[troughmask], troughZ[troughmask]

print('Selected:', len(RA), 'of', len(troughRA), 'galaxies', '(', float(len(RA))/float(len(troughRA)), 'percent )')

if 'mice' in cat:
    # Import source catalog

    # Path to the Mice field
    path_mockcat = '/data2/brouwer/MergedCatalogues'
    mockcatname = 'mice_catalog.fits'

    # Importing the Mice sources
    galRA, galDEC, galZ, rmag, rmag_abs, e1, e2 = \
    utils.import_mockcat(path_mockcat, mockcatname)

else:
    # Import source catalog

    # Path to the KiDS field
    path_srccat = '/data2/brouwer/KidsCatalogues/G9'
    srccatname = 'KiDS_G9_reweight_5x5x5_BLIND_PF.cat'

    # Importing the Mice sources
    galRA, galDEC, galZ, e1, e2, weight = \
    utils.import_srccat(path_srccat, srccatname)


# Masking the sources
srcmask = (0.1 < galZ) & (galZ < 0.9) & (rmag > 20.)
galRA, galDEC, galZ, rmag, e1, e2  = \
galRA[srcmask], galDEC[srcmask], galZ[srcmask], rmag[srcmask], e1[srcmask], e2[srcmask]


# Calculating the shear profile
troughcat = treecorr.Catalog(ra=RA, dec=DEC, ra_units='deg', dec_units='deg')
galcat = treecorr.Catalog(ra=galRA, dec=galDEC, ra_units='deg', dec_units='deg', g1=e1, g2=e2)

config = {'min_sep': Rmin, 'max_sep': Rmax, 'nbins': 20, 'sep_units': 'arcmin', 'verbose': 2}

ng = treecorr.NGCorrelation(config)
ng.process(troughcat,galcat)   # Compute the cross-correlation.
ng.write(outputname)     # Write out to a file.

shearfile = np.loadtxt(outputname).T

Rbins, gamma_t, gamma_x = [shearfile[0], shearfile[3], shearfile[4]]

plt.plot(Rbins, gamma_t)
plt.axhline(y=0., ls=':', color='black')

plt.xscale('log')
plt.axis([2,100,-1.4e-3,1.5e-3])
plt.ylim(-1.4e-3,1.5e-3)

plt.show()
