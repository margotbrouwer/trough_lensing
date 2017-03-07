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


# Radii theta of circular regions (in deg)
#thetalist = np.array([5., 10., 15., 20.])/60.
thetalist = np.array([5.])/60.
Ntheta = len(thetalist)

print('Theta:', thetalist*60., 'arcmin')
print()





### 1) Defining the galaxy sample:

## 1a) Importing the galaxy catalogue.

# Select the galaxy catalogue for trough selection (kids/gama)
cat = 'gama'

if cat == 'kids':
    
    # Names of the KiDS fields
    fieldnames = ['G9', 'G12', 'G15', 'G23', 'GS']
    
    # Boundaries of the KiDS fields
    coordsG9 = np.array([[128.0,142.5], [-2.5,3.5]])
    coordsG12 = np.array([[155.0,190.0], [-3.5,3.5]])
    coordsG15 = np.array([[209.5,239.0], [-3.5,3.5]])
    coordsG23 = np.array([[328.0,361.0], [-35.0,-28.0]])
    coordsGS = np.array([[31.0,54.0], [-35.0,-29.0]])
    fieldboundaries = np.array([coordsG9,coordsG12,coordsG15,coordsG23,coordsGS]) # Boundaries of all fields
    
    # Path to the KiDS fieldsc
    path_kidscat = '/data2/brouwer/KidsCatalogues'
    kidscatname = 'KiDS.DR3.zpbpzLess0.6.fits'
    
    # Importing the KiDS coordinates
    galRA, galDEC, galZ, galTB, mag_auto, ODDS, umag, gmag, rmag, imag = \
    utils.import_kidscat(path_kidscat, kidscatname)
    gridmax = 1./60.

if cat == 'gama':
    
    # Names of the GAMA fields
    fieldnames = ['G9', 'G12', 'G15']
    
    # Boundaries of the GAMA fields
    coordsG9 = [[129., 141.], [-2.,3.]]
    coordsG12 = [[174., 186.], [-3.,2.]]
    coordsG15 = [[211.5, 223.5], [-2.,3.]]
    fieldboundaries = np.array([coordsG9,coordsG12,coordsG15])

    # Path to the KiDS fields
    path_gamacat = '/data2/brouwer/MergedCatalogues/'
    gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'
    
    # Importing the GAMA coordinates
    galRA, galDEC, galZ, rmag, rmag_abs = utils.import_gamacat(path_gamacat, gamacatname)
    gridmax = 0.

    # Import GAMA masks for effective area calculation
    gridspace_mask = 0.01 # in degree
    path_gamamasks = ['/data2/brouwer/MergedCatalogues/GamaMasks/%smask08000.fits'%g for g in ['g09', 'g12', 'g15']]
    gamamasks = utils.import_gamamasks(path_gamamasks, gridspace_mask, fieldboundaries)
  
## 1b) Select the galaxy sample to define the troughs

# Name of the pre-defined galaxy selection [all, ell, redseq, redseq4]
selection = 'absmag'
#selection = 'redseq4'

# Defining the selection for the KiDS galaxy sample
if cat=='kids':

    # Redshift cut
    zmin = 0.1
    zmax = 0.4

    galmask = utils.define_galsamp(selection, zmin, zmax, galZ, galTB, gmag, rmag, mag_auto)

# Defining the selection for the GAMA galaxy sample
if cat=='gama':
    if selection == 'absmag':
        galmask = (rmag_abs < -19.7)&(rmag <= 19.8)
    if selection == 'all':
        galmask = (rmag <= 19.8)
    
print( 'Selected: %i/%i galaxies = %g percent'%(np.sum(galmask), len(galmask), float(np.sum(galmask))/float(len(galmask))*100.) )

# These lists will contain the full trough catalog (for all fields)
gridRA_tot = []
gridDEC_tot = []
Ngaltheta_tot = np.array([[]]*Ntheta)
Nmasktheta_tot = np.array([[]]*Ntheta)

# These lists will contain the density information for printing to a text file
field_area = [] # Will contain the masked area of the field (in arcmin)
field_galaxies = [] # Will contain the number of selected LRGs

# Looping over each field
for field in range(len(fieldnames)):
#for field in range(1):

    # Boundaries of the current KiDS field
    fieldRAs = fieldboundaries[field,0]
    fieldDECs = fieldboundaries[field,1]
    print(fieldRAs, fieldDECs)
    
    # Selecting the galaxies lying within this field
    fieldmask = (fieldRAs[0] < galRA)&(galRA < fieldRAs[1]) & (fieldDECs[0] < galDEC)&(galDEC < fieldDECs[1])

    print()
    print('KiDS field %s:'%fieldnames[field])

    # Applying the galaxy mask to the KiDS sample
    fieldmask = galmask*fieldmask

    # Coordinates of the galaxies
    galRA_field = galRA[fieldmask]
    galDEC_field = galDEC[fieldmask]
    galcoords = SkyCoord(ra=galRA_field*u.deg, dec=galDEC_field*u.deg)
    
    print('Selected: %g galaxies'%np.sum(fieldmask))
    field_galaxies = np.append(field_galaxies, np.sum(fieldmask))

    ### 2) Creating the grid:
    
    ## 2a) Creating a cartesian grid of narrowly spaced (1 arcmin) points.
    
    # Spacing of the grid (in deg)
    gridspace = 1./60. # One arcmin
    
    ## 2b) We remove all points that lie within 1 arcmin of masked areas or the edge of the field.
    
    # Define the grid coordinates in this field
    gridRA, gridDEC, gridcoords = utils.define_gridpoints(fieldRAs, fieldDECs, galcoords, gridspace, gridmax)
    Ngrid = len(gridcoords)
    gridID = np.arange(Ngrid)
    
    # Define the mask coordinates in this field
    maskRA, maskDEC, maskcoords = utils.define_gridpoints(fieldRAs, fieldDECs, galcoords, gridspace_mask, gridmax)
    gamamask = gamamasks[field]

    #n, bins, patches = plt.hist(gamamask, 50)
    #plt.show()
    
    maskcoords = maskcoords[gamamask>0]
    gamamask = gamamask[gamamask>0]
    
    
    Nmaskedgrid = np.sum(gamamask)
    
    # Field area information for printing to text file
    field_area = np.append(field_area, Nmaskedgrid)
    
    ### 3) Measuring galaxy density:
    
    # These lists will contain the number of galaxies/mask points counted in the circle around each grid point (for this field)
    Ngaltheta = np.zeros([Ntheta, Ngrid])
    Nmasktheta = np.zeros([Ntheta, Ngrid])
   
    # Creating smaller grid point samples (of 10.000 points) to avoid memory overload
    sampindex = np.array(np.append(np.arange(0., Ngrid, 1e4), Ngrid), dtype=int)
    
    # Looping over different circle radius sizes (theta)
    for theta in range(Ntheta):
        print('Theta =', thetalist[theta]*60., 'arcmin')

        # Looping over the smaller samples of grid points (to avoid memory overload)
        for s in range(len(sampindex)-1):
        #for s in range(1):
            
            print(fieldnames[field], ', Theta=' , thetalist[theta]*60., ', Grid sample: %i of %i'%(s+1, len(sampindex)-1))
            
            # Defining the smaller grid point sample
            gridsamp = gridcoords[sampindex[s]:sampindex[s+1]]
            
            # 3a) For each grid point, we count the number of sample galaxies within a circle of chosen radius \theta.
            galxgrid, gridxgal, d2d, d3d = galcoords.search_around_sky(gridsamp, thetalist[theta]*u.deg)
            galcountlist = np.array(list(Counter(galxgrid).items())).T
            (Ngaltheta[theta, sampindex[s]:sampindex[s+1]])[galcountlist[0]] = galcountlist[1]
            
            # 3b) For each grid point, we count the number of other non-zero mask-points to estimate the effective area.
            maskxgrid, gridxmask, d2d, d3d = maskcoords.search_around_sky(gridsamp, thetalist[theta]*u.deg)
            print('Calculating the mask area...')
            
            """
            # Simple masking: completeness is either 1 or 0.
            maskcountlist = np.array(list(Counter(maskxgrid).items())).T
            (Nmasktheta[theta, sampindex[s]:sampindex[s+1]])[maskcountlist[0]] = maskcountlist[1]
            """
            # Improved masking: completeness varies between 0 and 1.
            maskcountlist = np.array([np.sum(gamamask[gridxmask][maskxgrid==g]) for g in range(len(gridsamp))])
            Nmasktheta[theta, sampindex[s]:sampindex[s+1]] = maskcountlist
            #"""
            print('Maskcountlist:', maskcountlist)

    
    # Combine the field information into one final catalog
    gridRA_tot = np.append(gridRA_tot, gridRA)
    gridDEC_tot = np.append(gridDEC_tot, gridDEC)
    Ngaltheta_tot = np.hstack([Ngaltheta_tot, Ngaltheta])
    Nmasktheta_tot = np.hstack([Nmasktheta_tot, Nmasktheta])
    
# Calculate the percentage of the area that is masked
#area = np.pi*thetalist**2/gridspace_mask**2
Pmasktheta_tot = Nmasktheta_tot/np.amax(Nmasktheta_tot, axis=1)

# Define the total grid (of all fields)
Ngrid_tot = len(gridRA_tot)
gridID_tot = np.arange(Ngrid_tot)
gridcoords_tot = SkyCoord(ra=gridRA_tot*u.deg, dec=gridDEC_tot*u.deg)


#  3c) For each grid point, the galaxy density is defined as the galaxy count within the circle, divided by the effective area of the circle.

# Density
mask_density = 1/(gridspace_mask*60.)**2 # Density of mask gridpoints (in arcmin^-2)
eff_area = Nmasktheta_tot/mask_density # Effective area of each circle
rhotheta_tot = Ngaltheta_tot/eff_area
rho_mean = np.mean(rhotheta_tot, 1)
print('Mean density', rho_mean)


### 4) Flagging overlapping circles

## 4a) We sort the circles by galaxy density, and divide the circles in two samples: underdense and overdense.

# Sorting circles by density (rho) from lowest to highest
sort_rho = np.argsort(rhotheta_tot)
# Sorting overdense circles (second half of the sorted list) from highest to lowest
sort_rho_crab = [np.append( sort_rho[theta, 0:int(Ngrid_tot/2)], \
                    np.flip(sort_rho[theta, int(Ngrid_tot/2):Ngrid_tot], 0) ) for theta in range(Ntheta)]


# This array will contain the selected troughs (for all fields)
Stheta_tot = np.zeros([Ntheta, Ngrid_tot])
Ptheta_tot = np.zeros([Ntheta, Ngrid_tot])

for theta in range(Ntheta):

    ## 4b) Definition: Two circles are "overlapping" when their centers are closer than 0.5 \theta.
    thetamax = 0.5 * thetalist[theta]
    print('theta=', thetalist[theta]*60.)

    # This list will contain the selected (1) and flagged (0) troughs for this radius theta (for all fields)
    selected_theta = np.ones(Ngrid_tot)
    
    # Sort the troughs by density
    selID = gridID_tot[sort_rho[theta]]
    selcoords = gridcoords_tot[sort_rho[theta]]

    # Creating smaller grid point samples to avoid memory overload
    sampindex_tot = np.array(np.append(np.arange(0., Ngrid_tot, 1e4), Ngrid_tot), dtype=int)
    
    # Looping over the smaller samples of grid points (to avoid memory overload)
    for s in range(len(sampindex_tot)-1):
    #for s in range(1):
    
        print('Removing troughs within theta_max =' , thetamax*60., ', Grid sample: %i of %i'%(s+1, len(sampindex_tot)-1))
        
        # Defining the smaller grid point sample
        selsamp = selcoords[sampindex_tot[s]:sampindex_tot[s+1]]
        
        # 4c/d) For the underdense/overdense circles: We start from the lowest/highest density circle and flag all overlapping circles.
        
        # Find grid points within thetamax (= 0.5 theta) of each other
        sampxgrid, gridoverlap, d2d, d3d = selcoords.search_around_sky(selsamp, thetamax*u.deg)
        
        # Looping over grid points
        for g in range(len(selsamp)):
            
            # Index of the grid point in the flag list
            tot_index = sampindex_tot[s]+g
            
            # If the point is not flagged (not flagged = 1), flag all overlapping grid points (flagged = 0).
            if selected_theta[tot_index] == 1:
                removeindex = gridoverlap[(sampxgrid==g)&(gridoverlap!=tot_index)]
                selected_theta[removeindex] = 0.
            
            # If the point is already flagged, do nothing.
            else:
                pass
    
    # Sort the grid points back to their original order
    sort_ID = np.argsort(selID)

    # The list Stheta contains the selected (= not flagged = 1) grid points for the non-overlapping sample.
    Stheta_tot[theta] = selected_theta[sort_ID]


    ### 5) Selecting the trough sample:
    
    ##  5a) For each grid point/circle, we calculate the percentage Ptheta of grid points that has a lower galaxy density.
    Ptheta_tot[theta] = np.argsort(gridID_tot[sort_rho[theta]])/Ngrid_tot
    
    # Following the definition from DES, we usually define the 20% of all circles with the lowest density to be the troughs.
    # The 20% of all circles with the highest density are defined as overdensities ("ridges").


## Write catalog

## 5b) Now we have two trough samples:
# - The overlapping sample, by taking all troughs
# - The non-overlapping sample, by taking only unflagged troughs (not flagged = 1, flagged = 0).

# Writing the combined columns to a fits file
filename = '/data2/brouwer/MergedCatalogues/trough_catalog_%s_%s_masked.fits'%(cat, selection)

# For each grid point/circle, we save the following information to the catalog:
# the location (RA/DEC), number count (Ngaltheta), effective area in arcmin (grid count Ngtheta), galaxy density (rhotheta), 
# percentage of circles below its density (Ptheta), and non-overlapping selection flag (Stheta).

# Defining the names of the columns
outputnames = ['RA', 'DEC']
[outputnames.append('Ngaltheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Nmasktheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Pmasktheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('rhotheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Ptheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Stheta%i'%(theta*60)) for theta in thetalist]

# Defining the output
output = [gridRA_tot, gridDEC_tot]
[output.append(Ngaltheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(Nmasktheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(Pmasktheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(rhotheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(Ptheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(Stheta_tot[theta,:]) for theta in range(Ntheta)]


print('Writing output catalogue...')
utils.write_catalog(filename, gridID_tot, outputnames, output)


# Writing the mean galaxy density information to a text file
field_area = np.append(field_area, np.mean(field_area))
field_galaxies = np.append(field_galaxies, np.mean(field_galaxies))

field_density = field_galaxies/field_area

field_header = 'G9      G12      G15      G23      GS      Average'
field_footer = '1: Area after masking (arcmin^2), 2: Number of red galaxies, 3: Galaxy density (arcmin^-2)'
density_info = np.array([field_galaxies, field_area, field_density])

filename = 'density_info_%s.txt'%selection

np.savetxt(filename, density_info, delimiter='    ',  fmt='%.6g', header = field_header, footer = field_footer)
print('Written:', filename)


# Write mask fits-file
outputnames = ['RA', 'DEC', 'Completeness']
output = [maskRA, maskDEC, gamamask]
utils.write_catalog('/data2/brouwer/MergedCatalogues/gama_mask.fits', np.arange(len(maskRA)), outputnames, output)
