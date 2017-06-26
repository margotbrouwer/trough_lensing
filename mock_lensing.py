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
from astropy.cosmology import LambdaCDM
from collections import Counter

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

import trough_modules_all as utils
import treecorr

# Import constants
G = const.G.to('pc3/Msun s2')
c = const.c.to('pc/s')

# Configuration

cat = 'mice' # Select the galaxy catalogue for trough selection (kids, gama or mice)
selection = 'all' # Name of the pre-defined galaxy selection
masktype = 'nomask' # Select mask type (nomask or complex)

h, O_matter, O_lambda = [0.7, 0.25, 0.75]
cosmo = LambdaCDM(H0=h*100, Om0=O_matter, Ode0=O_lambda)
Runit = 'arcmin'

# Defining the trough radii
thetalist = np.array([5., 10., 15., 20.]) # in arcmin
if 'lowZ' in selection:
    thetalist = np.array([10.])
    Runit = 'Mpc'
if 'highZ' in selection:
    thetalist = np.array([6.326])
    Runit = 'Mpc'

# Trough selection

# Percentiles
thetanum = 3

if 'arcmin' in Runit:
    Rmin = 2
    Rmax = 100
    Nbins = 20
    
    ymin, ymax = [-1.4e-3,1.5e-3]
    
    dperc = 0.05
    percnames = ['0','0p05','0p1','0p15','0p2','0p25','0p3','0p35','0p4',\
    '0p45','0p5','0p55','0p6','0p65','0p7','0p75','0p8','0p85','0p9','0p95','1']

if 'pc' in Runit:
    Rmin = 0.5
    Rmax = 20
    Nbins = 10
    
    ymin, ymax = [-1.4,1.5]
    
    dperc = 0.1
    percnames = ['0','0p1','0p2','0p3','0p4','0p5','0p6','0p7','0p8','0p9','1']
    Rlist = [1.8825293] # Physical size of the troughs (in Mpc)

# Defining the radial bins
Rbins = np.logspace(np.log10(Rmin), np.log10(Rmax), Nbins+1)
Rcenters = Rbins[0:-1] + np.diff(Rbins)/2.

# Defining the percentile bins
percmin = 0.
percmax = 1.

perclist = np.arange(percmin, percmax, dperc)
perccenters = perclist+dperc/2.
perclist = np.append(perclist, percmax)

Npercs = len(percnames)-1

"""
# Percentiles
theta = thetalist[thetanum]
paramnames_tot = np.array([ ['Pmasktheta%g'%theta, 'Ptheta%g'%theta] for p in range(Npercs) ])
maskvals_tot = np.array([ [[0.8, np.inf], [perclist[p], perclist[p+1]]] for p in range(Npercs) ])

"""

# Sizes
paramnames_tot = np.array([ ['Pmasktheta%g'%theta, 'Ptheta%g'%theta] for theta in thetalist])
maskvals_tot = np.array([ [[0.8, np.inf], [0., 0.2]] for theta in thetalist ])

#"""


# Import trough catalog
path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'
troughcatname = 'trough_catalog_%s_%s_%s.fits'%(cat, selection, masktype)
troughRA, troughDEC, troughZ, paramlists = utils.import_troughcat(path_troughcat, troughcatname, [])
troughZ = troughZ[0]
troughDc = (cosmo.comoving_distance(troughZ).to('pc')).value
troughDa = troughDc / (1+troughZ)

if 'pc' in Runit:
    # Translate radial distances from Mpc to arcmin
    Rarcmin = np.degrees(Rmin/(troughDa/10**6.))*60.
    Rarcmax = np.degrees(Rmax/(troughDa/10**6.))*60.
else:
    Rarcmin = Rmin
    Rarcmax = Rmax

# Import source catalog
if 'mice' in cat:
    # Path to the Mice field
    path_mockcat = '/data2/brouwer/MergedCatalogues'
    mockcatname = 'mice_catalog.fits'

    # Importing the Mice sources
    galRA, galDEC, galZ, rmag, rmag_abs, e1, e2 = \
    utils.import_mockcat(path_mockcat, mockcatname)

else:
    # Path to the KiDS field
    path_srccat = '/data2/brouwer/KidsCatalogues/G9'
    srccatname = 'KiDS_G9_reweight_5x5x5_BLIND_PF.cat'

    # Importing the Mice sources
    galRA, galDEC, galZ, e1, e2, weight = \
    utils.import_srccat(path_srccat, srccatname)

# Masking the sources
srcmask = (0.1 < galZ) & (galZ < 0.9) & (rmag > 20.) & (rmag_abs > -19.3)
if 'pc' in Runit:
    # Only use sources behind the lens
    zmask = (galZ > troughZ)
    srcmask = srcmask * zmask

galRA, galDEC, galZ, rmag, e1, e2  = \
galRA[srcmask], galDEC[srcmask], galZ[srcmask], rmag[srcmask], e1[srcmask], e2[srcmask]

if 'pc' in Runit:
    # Calculate source distances
    print('Calculating source distances...')

    galbins = np.arange(0, len(galRA), 1e5)
    galbins = np.append(galbins, len(galRA))

    galDc = np.zeros(len(galRA))
    for i in range(len(galbins)-1):
        print([galbins[i], galbins[i+1]])
        galDc[galbins[i]:galbins[i+1]] = (cosmo.comoving_distance(galZ[galbins[i]:galbins[i+1]]).to('pc')).value
        
    print('Calculated source distances')

    # Calculate Sigma_crit for every source (there is only one lens distance)
    DlsoDs = (galDc - troughDc)/galDc
    Sigma_crit = (c.value**2)/(4*np.pi*G.value) * 1/(troughDa*DlsoDs)


for p in range(len(paramnames_tot)):

    paramnames = paramnames_tot[p]
    maskvals = maskvals_tot[p]
    
    # Import trough parameters
    path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'
    troughcatname = 'trough_catalog_%s_%s_%s.fits'%(cat, selection, masktype)
    foo, foo, foo, paramlists = utils.import_troughcat(path_troughcat, troughcatname, paramnames)
    troughweights = np.ones(len(troughRA))

    
    # Selecting the troughs
    masklist = np.ones(len(troughRA))
    filename_var = ''
    
    # Applying the mask for each parameter
    for m in range(len(maskvals)):
        troughmask = (maskvals[m,0] <= paramlists[m]) & (paramlists[m] < maskvals[m,1])
        masklist[np.logical_not(troughmask)] = 0
        
        filename_var = '%s~%s_%g_%g'%(filename_var, paramnames[m], maskvals[m,0], maskvals[m,1])
        
    filename_var = filename_var.replace('.', 'p')
    filename_var = filename_var.replace('-', 'm')
    filename_var = filename_var.replace('~', '-')
    filename_var = filename_var.split('-', 1)[1]
    
    troughmask = (masklist == 1)
    RA, DEC = troughRA[troughmask], troughDEC[troughmask]
    
    print()
    print(filename_var)
    print('Selected:', len(RA), 'of', len(troughRA), 'lenses', '(', float(len(RA))/float(len(troughRA))*100., 'percent )')
    print()
    
    if 'pc' in Runit:
        
        # Calculating the ESD profile
        galcat = treecorr.Catalog(ra=galRA, dec=galDEC, ra_units='deg', dec_units='deg', g1=e1*Sigma_crit, g2=e2*Sigma_crit, w=1/Sigma_crit**2.)
        
    else:
        # Calculating the shear profile
        galcat = treecorr.Catalog(ra=galRA, dec=galDEC, ra_units='deg', dec_units='deg', g1=e1, g2=e2)
    
    troughcat = treecorr.Catalog(ra=RA, dec=DEC, ra_units='deg', dec_units='deg', w=troughweights)
    
    config = {'min_sep': Rarcmin, 'max_sep': Rarcmax, 'nbins': Nbins, 'sep_units': 'arcmin', 'verbose': 2}
    ng = treecorr.NGCorrelation(config)
    ng.process(troughcat,galcat)   # Compute the cross-correlation.

    output_temp = 'treecorr_temp.txt'
    ng.write(output_temp)     # Write out to a file.
    shearfile = np.loadtxt(output_temp).T

    Rbins, gamma_t, gamma_x, gamma_error, Nsrc = \
    [shearfile[0], shearfile[3], shearfile[4], np.sqrt(shearfile[5]), shearfile[7]]

    path_output = '/data2/brouwer/shearprofile/trough_results_May/No_bins_%s_%s_%s'%(cat, selection, masktype)
    filename_output = '%s/%s.txt'%(path_output, filename_var)

    if not os.path.isdir(path_output):
        os.makedirs(path_output)
        print 'Creating new folder:', path_output
    
    bias = np.ones(len(gamma_t))
    gamma_error = np.zeros(len(gamma_t))
    utils.write_stack(filename_output, Rcenters, Runit, gamma_t, gamma_x, \
        gamma_error, bias, h, Nsrc)

    """
    
    # Plot the resulting shear profile
    plt.plot(Rcenters, gamma_t)
    plt.axhline(y=0., ls=':', color='black')
    #plt.axvline(x=thetalist[p], ls=':', color='black')
    
    plt.xscale('log')
    
    plt.axis([Rmin,Rmax,ymin,ymax])
    plt.ylim(ymin, ymax)

    plt.show()
    
    """
