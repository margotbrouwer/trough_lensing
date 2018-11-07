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
import matplotlib.style
matplotlib.style.use('classic')

import trough_modules_all as utils

# Make use of TeX
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})

#colors = ['red', 'orange', 'cyan', 'blue']
colors = ['#d7191c', '#fdae61', '#92c5de', '#0571b0']

h, O_matter, O_lambda = [0.7, 0.25, 0.75]
micecor = 5*np.log10(h) # Correction on MICE absmag_r (= -0.7745)

# Path to the GAMA fields
path_gamacat = '/data2/brouwer/MergedCatalogues/'
gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'

# Importing the GAMA coordinates
gamaRA, gamaDEC, gamaZ, rmag_gama, rmag_abs_gama = utils.import_gamacat(path_gamacat, gamacatname)
gamamask = (gamaZ < 0.5) & (rmag_abs_gama < -18.9 + micecor) & (rmag_gama < 19.8)


# Path to the KiDS fields
path_kidscat = '/data2/brouwer/KidsCatalogues'
kidscatname = 'KiDS_DR3_GAMA-like_Maciek_revised_1905.fits'

# Importing the KiDS coordinates
kidsRA, kidsDEC, kidsZB, kidsZ, kidsTB, mag_auto, ODDS, umag_kids, gmag_kids, rmag_kids, imag_kids = \
utils.import_kidscat(path_kidscat, kidscatname)
rmag_abs_kids = utils.calc_absmag(rmag_kids, kidsZ, gmag_kids, imag_kids, h, O_matter, O_lambda)

kidsmask = (kidsZ < 0.5) & (rmag_abs_kids < -18.9 + micecor) & (rmag_kids < 20.2)

# Path to the MICE mocks
path_micecat = '/data2/brouwer/MergedCatalogues'
micecatname = 'mice_gama_catalog.fits'

# Importing the MICE galaxies
miceRA, miceDEC, miceZ, miceDc, rmag_mice, rmag_abs_mice, e1_mice, e2_mice, galmass_mice = \
utils.import_mockcat(path_micecat, micecatname)
rmag_abs_slics = rmag_abs_mice + micecor

micemask = (miceZ < 0.5) & (rmag_abs_mice < -18.9 + micecor) & (rmag_mice < 20.2)

# Importing SLICS redshifts
path_slicsdata = 'Vasiliy_SN_files/nofz_GAMA_LOS74_199.txt'
slicsdata = np.loadtxt(path_slicsdata).T
zbins_slics = slicsdata[0]
zhist_slics = slicsdata[1]

# Importing KiDS source redshifts
path_sourcedata = 'Vasiliy_SN_files/kids_data_nz.txt'
sourcedata = np.loadtxt(path_sourcedata).T

zbins_sources = np.linspace(0.,3.5,70)
zhist_sources = np.array(sourcedata[0])

dz = 3.5/70.
zhist_sources = zhist_sources/np.sum(zhist_sources*dz)

np.savetxt('kids450_nz.txt', np.array([zbins_sources, zhist_sources]).T, fmt = ['%g', '%.10e'], \
delimiter = '      ', header='KiDS-450: z        n(z)')

# Plot redshift histograms of KiDS, GAMA, MICE and SLICS
print('MICE average:', np.mean(miceZ[micemask]), (np.mean(miceZ[micemask])-np.mean(gamaZ[gamamask]))/np.mean(gamaZ[gamamask]))
print('GL-KiDS average:', np.mean(kidsZ[kidsmask]), (np.mean(kidsZ[kidsmask])-np.mean(gamaZ[gamamask]))/np.mean(gamaZ[gamamask]))
print('GAMA average:', np.mean(gamaZ[gamamask]))

lw = 1.
fig = plt.figure(figsize=(10,4))

galsamps = [gamaZ[gamamask], kidsZ[kidsmask], miceZ[micemask]]
galnames = [r'GAMA galaxies', r'GL-KiDS galaxies', r'GL-MICE mocks', r'SLICS mocks']

for t in range(len(galsamps)):
    n, bins, patches = plt.hist(galsamps[t], bins=60, histtype='step', range=[0., 0.5], \
                normed=1., alpha=1., color=colors[t], linewidth=lw)
    #plt.axvline( x=np.mean(galsamps[t]), ls = '--', color=colors[t])

plt.plot(zbins_slics, zhist_slics, drawstyle='steps', color=colors[t+1], linewidth=lw)
plt.plot(zbins_sources[zbins_sources<1.4], zhist_sources[zbins_sources<1.4], ls='--', color='black', linewidth=lw)

plt.xlabel(r'Redshift $z$', fontsize=14)
plt.ylabel(r'Number of galaxies (normalized)', fontsize=14)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = [Patch(facecolor='none', edgecolor=colors[t], linewidth=lw, label=galnames[t]) for t in range(len(galnames))]
legend_elements.append(Line2D([0], [0], ls='--', color='black', linewidth=lw, label=r'KiDS source $n(z_{\rm s})$'))
plt.legend(handles=legend_elements, loc='best')

plotfilename = '/data2/brouwer/shearprofile/Lensing_results/trough_results_final/Plots/redshift_distribution'
plt.tight_layout()

# Save plot
for ext in ['pdf']:
    plotname = '%s.%s'%(plotfilename, ext)
    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written: ESD profile plot:', plotname)

plt.show()
plt.clf
