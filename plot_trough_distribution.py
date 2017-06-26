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

#colors = ['blue', 'cyan', 'red', 'orange']
colors = ['#d7191c', '#fdae61', '#92c5de', '#0571b0']


# Radii theta of circular regions (in deg)
thetalist = np.array([5., 10., 15., 20.])/60.

selections = ['kids_19p8_complex', 'mice_all_nomask']

fig = plt.figure(figsize=(5,4))
for s in range(len(selections)):
    
    path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'
    
    # Full directory & name of the trough catalogue
    troughcatfile = '%s/trough_catalog_%s.fits'%(path_troughcat, selections[s])
    troughcat = pyfits.open(troughcatfile, memmap=True)[1].data

    # List of the observables of all sources in the KiDS catalogue
    RA = troughcat['RA']
    DEC = troughcat['DEC']

    Ngaltheta = [troughcat['Ngaltheta%g'%(theta*60)] for theta in thetalist]
    rhotheta = [troughcat['rhotheta%g'%(theta*60)] for theta in thetalist]
    delta = [troughcat['delta%g'%(theta*60)] for theta in thetalist]
    Pmasktheta = [troughcat['Pmasktheta%g'%(theta*60)] for theta in thetalist]

    for t in range(len(thetalist)):
        
        Ngal = (rhotheta[t])[Pmasktheta[t]>0.8]
        
        if ('kids' in selections[s]) or ('gama' in selections[s]):
        # Plot the observed histograms
            hist, bin_edges, patches = plt.hist(Ngal, bins=25, range=[0., 0.7], normed=1, \
            label=r"$\theta_{\rm A} = %g'$"%(thetalist[t]*60.), histtype='step', color=colors[t])
            
            plt.axvline( x=np.mean((rhotheta[t])[Pmasktheta[t]>0.8]), ls = '--', color=colors[t])
        
        else:
        # Plot the mock histograms
            hist, bin_edges = np.histogram(Ngal, bins=25, range=[0., 0.7], normed=1)
            bin_centers = bin_edges[0:-1] + np.diff(bin_edges)/2.
            
            plt.plot(bin_centers, hist, ls = '-', color=colors[t], alpha=1.)
        
        

plt.xlabel(r'Galaxy number density $n_{\rm g}(\theta_{\rm A})$ (${\rm arcmin}^{-2}$)', fontsize=14)
plt.ylabel(r'Normalized number of apertures', fontsize=14)

#plt.yscale('log')
#plt.yscale('log')
plt.legend()

plotfilename = '/data2/brouwer/shearprofile/trough_results_June/Plots/trough_density_distribution'

plt.tight_layout()

# Save plot
for ext in ['pdf']:
    plotname = '%s.%s'%(plotfilename, ext)
    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written: ESD profile plot:', plotname)

plt.show()
plt.clf
