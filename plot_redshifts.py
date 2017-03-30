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


# Path to the GAMA fields
path_gamacat = '/data2/brouwer/MergedCatalogues/'
gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'

# Importing the GAMA coordinates
galRA, galDEC, galZ, rmag, rmag_abs = utils.import_gamacat(path_gamacat, gamacatname)


galmask = (-24.9 < rmag_abs) & (rmag_abs < -12.)
galRA, galDEC, galZ, rmag, rmag_abs = galRA[galmask], galDEC[galmask], galZ[galmask], rmag[galmask], rmag_abs[galmask]

plt.axhline(y=-19.7, ls='--', color='red', label= r'Fiducial sample $(M_r<-19.7)$')

plt.axhline(y=-21, ls='--', color='black', label = r'Volume limited sample $(M_r<-21)$')
plt.axvline(x=0.2324, ls='--', color='black', label = r'$Z_{\rm lim} = 0.2324$')
plt.axvline(x=0.3, ls='--', color='black', label = r'$Z_{\rm max} = 0.3$')

plt.hist2d(galZ, rmag_abs, bins=100, norm=LogNorm())
plt.colorbar()

plt.xlim(0,0.4)
plt.ylim(-25,-13)

plt.legend(loc='lower right')

plt.gca().invert_yaxis()

ext = 'pdf'
plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/redshift_limit'
plotname = '%s.%s'%(plotfilename, ext)

plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()
