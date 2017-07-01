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

#colors = ['red', 'orange', 'cyan', 'blue']
colors = ['#d7191c', '#fdae61', '#92c5de', '#0571b0']


# Path to the GAMA fields
path_gamacat = '/data2/brouwer/MergedCatalogues/'
gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'

# Importing the GAMA coordinates
galRA, galDEC, galZ, rmag, rmag_abs = utils.import_gamacat(path_gamacat, gamacatname)
galmask = (-24.9 < rmag_abs) & (rmag_abs < -14.)


# Path to the KiDS fields
path_kidscat = '/data2/brouwer/KidsCatalogues'
kidscatname = 'KiDS_DR3_GAMA-like_Maciek_revised_1905.fits'

# Importing the KiDS coordinates
kidsRA, kidsDEC, kidsZB, kidsZ, kidsTB, mag_auto, ODDS, umag_kids, gmag_kids, rmag_kids, imag_kids = \
utils.import_kidscat(path_kidscat, kidscatname)


# Plot absmag-Z diagram of GAMA
fig = plt.figure(figsize=(6,5))

plt.axhline(y=-19.67, ls='-', linewidth=2., color='black', label= r'Fiducial $(M_{\rm r}<-19.67)$', zorder=3)
plt.axhline(y=-21, ls='-', linewidth=2., color='blue', label = r'Volume limited $(M_{\rm r}<-21)$', zorder=2)

plt.axvline(x=0.1, ls='-', linewidth=2., color='green', label = r'$ \{ z_{\rm min}, z_{\rm lim}, z_{\rm max} \} $', zorder=1)
plt.axvline(x=0.198, ls='-', linewidth=2., color='green', zorder=1)
plt.axvline(x=0.3, ls='-', linewidth=2., color='green', zorder=1)


plt.hist2d(galZ[galmask], rmag_abs[galmask], bins=100, norm=LogNorm())
plt.colorbar()

plt.xlim(0,0.5)
plt.ylim(-25,-14)

plt.xlabel(r'Redshift $z$', fontsize=14)
plt.ylabel(r'Absolute magnitude $M_{\rm r}$', fontsize=14)

plt.legend(loc='lower right')

plt.gca().invert_yaxis()

ext = 'pdf'
plotfilename = '/data2/brouwer/shearprofile/trough_results_July/Plots/redshift_limit'
plotname = '%s.%s'%(plotfilename, ext)

plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()
plt.close()


# Plot redshift histograms of KiDS and GAMA
fig = plt.figure(figsize=(5,4))

galsamps = [rmag[rmag<=19.8], rmag_kids]
galnames = ['GAMA galaxies', 'GAMA-like KiDS']

for t in range(2):
    n, bins, patches = plt.hist(galsamps[t], bins=30, histtype='step', range=[15., 20.5], \
                normed=1., label=galnames[t], alpha=1., color=colors[t])
    #plt.axvline( x=np.mean(galsamps[t]), ls = '--', color=colors[t])


plt.xlabel(r'Magnitude $m_{\rm r}$', fontsize=14)
plt.ylabel(r'Number of galaxies (normalized)', fontsize=14)

#plt.axis([0.,0.6,-3,9])
plt.ylim(0.,1.)

#plt.yscale('log')
#plt.yscale('log')
plt.legend(loc='upper left')

plotfilename = '/data2/brouwer/shearprofile/trough_results_Apr/Plots/magnitude_distribution'

plt.tight_layout()

# Save plot
for ext in ['pdf']:
    plotname = '%s.%s'%(plotfilename, ext)
    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written: ESD profile plot:', plotname)

plt.show()
plt.clf

print(len(kidsZ))
