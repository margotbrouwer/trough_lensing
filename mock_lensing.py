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
inf = np.inf
   
h, O_matter, O_lambda = [0.7, 0.25, 0.75]
cosmo = LambdaCDM(H0=h*100, Om0=O_matter, Ode0=O_lambda)


# Defining the trough radii
thetalist = np.array([5., 10., 15., 20.]) # in arcmin
#thetalist = np.array([5.]) # in arcmin

# Defining mock patches
ijlist = np.array([ [ [i, j] for i in range(4) ] for j in range(4) ])
ijlist = np.reshape(ijlist, [16,2])

#Nruns = len(ijlist) # Mock patches
#Nruns = len(thetalist) # Theta
#Nruns = 5
Nruns = 1

# Configuration

for ij in np.arange(0, Nruns):
    
    # Number of the current run
    ijnum = ij+1

    # Select the galaxy catalogue for trough selection
    cat = 'mice'
    
    # Name of the pre-defined galaxy selection
    selection = 'all'
    #selection = 'lowZ'
    #selection = 'miceZ-%g'%ijnum
    
    # Select mask type
    if Nruns > 14:
        masktype = 'nomask-%g'%ijnum
        i, j = ijlist[ij]
        thetanum = 0
    else:
        masktype = 'nomask-1'
        thetanum = ij
    if 'miceZ' in selection:
        masktype = 'nomask-Z'
    
    
    # Import trough catalog
    path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'
    troughcatname = 'trough_catalog_%s_%s_%s.fits'%(cat, selection, masktype)
    troughRA, troughDEC, troughZ, paramlists = utils.import_troughcat(path_troughcat, troughcatname, [])
    troughZ = troughZ[0]
    troughDc = (cosmo.comoving_distance(troughZ).to('pc')).value
    troughDa = troughDc / (1+troughZ)
    
    print('Troughs:', len(troughRA))
    
    # Weights of the troughs
    troughweights = np.ones(len(troughRA))
    
    
    
    # Redshift samples
    if 'lowZ' in selection:
        thetalist = np.array([10.])
    if 'highZ' in selection:
        thetalist = np.array([6.826])
        
    if 'miceZ' in selection:
        thetalist = np.array([20., 12.85, 9.45, 7.44, 6.14]) # Dc
    
    # Select unit (arcmin or Mpc)
    if 'Z' in selection:
        Runit = 'Mpc'
    else:
        Runit = 'arcmin'
    
    if 'arcmin' in Runit:
        Rmin = 2
        Rmax = 100
        Nbins = 20
    
    if 'pc' in Runit:
        Rmin = 0.5
        Rmax = 20
        Nbins = 10
        
    
    # Defining the radial bins
    Rbins = np.logspace(np.log10(Rmin), np.log10(Rmax), Nbins+1)
    Rcenters = Rbins[0:-1] + np.diff(Rbins)/2.
    
    if 'pc' in Runit:
        # Translate radial distances from Mpc to arcmin
        Rarcmin = np.degrees(Rmin/(troughDa/10**6.))*60.
        Rarcmax = np.degrees(Rmax/(troughDa/10**6.))*60.
    else:
        Rarcmin = Rmin
        Rarcmax = Rmax
    
    
    ### Trough selection

    """
    """
    
    ## Percentiles
    if 'arcmin' in Runit:
        
        ymin, ymax = [-1.4e-3,1.5e-3]
        
        dperc = 0.05
        percnames = ['0','0p05','0p1','0p15','0p2','0p25','0p3','0p35','0p4',\
        '0p45','0p5','0p55','0p6','0p65','0p7','0p75','0p8','0p85','0p9','0p95','1']

    if 'pc' in Runit:
        
        ymin, ymax = [-1.4,1.5]
        
        dperc = 0.1
        percnames = ['0','0p1','0p2','0p3','0p4','0p5','0p6','0p7','0p8','0p9','1']

    # Defining the percentile bins
    percmin = 0.
    percmax = 1.

    perclist = np.arange(percmin, percmax, dperc)
    perccenters = perclist+dperc/2.
    perclist = np.append(perclist, percmax)

    Npercs = len(percnames)-1
    
    
    theta = thetalist[thetanum]
    paramnames_tot = np.array([ ['Pmasktheta%g'%theta, 'Ptheta%g'%theta] for p in range(Npercs) ])
    maskvals_tot = np.array([ [[0.8, np.inf], [perclist[p], perclist[p+1]]] for p in range(Npercs) ])
    
    """
    
    ## Fiducial troughs
    
    # Sizes
    
    paramnames_tot = np.array([ ['Pmasktheta%g'%theta, 'Ptheta%g'%theta] for theta in thetalist])
    maskvals_tot = np.array([ [[0.8, np.inf], [0., 0.2]] for theta in thetalist ])
    

    
    # Troughs/ridges
    
    perclist = [ [0., 0.2], [0.8, 1.] ]
    perclist = [ [0.8, 1.] ]
    theta = thetalist[thetanum]
    
    paramnames_tot = np.array([ ['Pmasktheta%g'%theta, 'Ptheta%g'%theta] for p in range(len(perclist)) ])
    maskvals_tot = np.array([ [[0.8, np.inf], perclist[p] ] for p in range(len(perclist)) ])
    

    
    ## Weighted profiles
    
    # Redshifts
    weightcatname = 'amplitude_trough_weights_%s_lowZ_%s.fits'%(cat, masktype)
    
    # Sizes
    #weightcatname = 'amplitude_trough_weights_%s_%s_%s.fits'%(cat, selection, masktype)
    
    perclist = [ [-inf , 0.], [0., inf] ]
    theta = thetalist[thetanum]
    
    paramnames_tot = np.array([ ['Pmasktheta%g'%theta, 'delta%g'%theta] for p in range(len(perclist)) ])
    maskvals_tot = np.array([ [[0.8, np.inf], perclist[p] ] for p in range(len(perclist)) ])

    weightfile = '%s/%s'%(path_troughcat, weightcatname)
    print('Weighted by:', weightfile, theta)
    
    weightcat = pyfits.open(weightfile, memmap=True)[1].data
    troughweights = weightcat['Wtheta%g'%theta]
    weightPthetas = weightcat['Ptheta%g'%theta]
    

    
    """

    # Import source catalog

    # Path to the Mice field
    path_mockcat = '/data2/brouwer/KidsCatalogues'
    mockcatname = 'mice_source_catalog_dc.fits'

    # Importing the Mice sources
    galRA, galDEC, galZ, galDc, rmag, rmag_abs, e1, e2, galmass = \
    utils.import_mockcat(path_mockcat, mockcatname)

    # Boundaries of the field

    if Nruns > 14:
        fieldRAs, fieldDECs = [[i*20.,(i+1)*20.], [j*20.,(j+1)*20.]]
    else:
        fieldRAs, fieldDECs =  [[0.,20.], [0.,20.]]
    
    print(ijnum, fieldRAs, fieldDECs)
    
    # Selecting the galaxies lying within this field
    fieldmask = (fieldRAs[0] < galRA)&(galRA < fieldRAs[1]) & (fieldDECs[0] < galDEC)&(galDEC < fieldDECs[1])

    # Masking the sources
    #srcmask = (0.1 < galZ) & (galZ < 0.9) & (rmag > 20.) & (rmag_abs > -19.3)
    srcmask = (troughZ+0.2 < galZ) & (galZ < 0.9) & (rmag > 20.) & (rmag_abs > -19.3)
    #srcmask = (0.6 < galZ) & (galZ < 0.7) & (rmag > 20.)# & (rmag_abs > -19.3)
    srcmask = srcmask*fieldmask
    
    if 'pc' in Runit:
        # Only use sources behind the lens
        zmask = (galZ > troughZ)
        srcmask = srcmask * zmask

    galRA, galDEC, galZ, galDc, rmag, e1, e2  = \
    galRA[srcmask], galDEC[srcmask], galZ[srcmask], galDc[srcmask], rmag[srcmask], e1[srcmask], e2[srcmask]
    galDc = galDc*1e6 # Convert distances from Mpc to pc
    print('Number of sources:', len(galZ))
    
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
        RA, DEC, weights = troughRA[troughmask], troughDEC[troughmask], troughweights[troughmask]
        
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
        
        troughcat = treecorr.Catalog(ra=RA, dec=DEC, ra_units='deg', dec_units='deg', w=weights)
        
        config = {'min_sep': Rarcmin, 'max_sep': Rarcmax, 'nbins': Nbins, 'sep_units': 'arcmin', 'verbose': 2}
        ng = treecorr.NGCorrelation(config)
        ng.process(troughcat,galcat)   # Compute the cross-correlation.

        output_temp = 'temp_treecor.txt'
        ng.write(output_temp)     # Write out to a file.
        shearfile = np.loadtxt(output_temp).T

        Rbins, gamma_t, gamma_x, gamma_error, Nsrc = \
        [shearfile[0], shearfile[3], shearfile[4], np.sqrt(shearfile[5]), shearfile[7]]

        # Translate to comoving ESD
        Rbins = Rbins*(1+troughZ)
        gamma_t, gamma_x, gamma_error = np.array([gamma_t, gamma_x, gamma_error])/(1+troughZ)**2

        path_output = '/data2/brouwer/shearprofile/trough_results_final/No_bins_%s_%s_%s'%(cat, selection, masktype)
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
