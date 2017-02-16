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

# Names of the KiDS fields
kidsfields = ['G9', 'G12', 'G15', 'G23', 'GS']

# Boundaries of the KiDS fields
coordsG9 = np.array([[128.0,142.5], [-2.5,3.5]])
coordsG12 = np.array([[155.0,190.0], [-3.5,3.5]])
coordsG15 = np.array([[209.5,239.0], [-3.5,3.5]])
coordsG23 = np.array([[328.0,361.0], [-35.0,-28.0]])
coordsGS = np.array([[31.0,54.0], [-35.0,-29.0]])
kidsboundaries = np.array([coordsG9,coordsG12,coordsG15,coordsG23,coordsGS]) # Boundaries of all fields

# Path to the KiDS fields
path_kidscat = '/data2/brouwer/KidsCatalogues'
kidscatname = 'KiDS.DR3.zpbpzLess0.6.fits'

# Importing the KiDS coordinates
galRA, galDEC, galZB, galTB, mag_auto, ODDS, umag, gmag, rmag, imag = \
utils.import_kidscat(path_kidscat, kidscatname)

## 1b) Select the galaxy sample to define the troughs

# Name of the pre-defined galaxy selection [all, ell, redseq, redseq4]
selection = 'all'
#selection = 'redseq4'

# Redshift cut
zmin = 0.1
zmax = 0.4

# Defining the mask for the galaxy sample
galmask = utils.define_galsamp(selection, zmin, zmax, galZB, galTB, gmag, rmag, mag_auto)


# These lists will contain the full trough catalog (for all fields)
gridRA_tot = []
gridDEC_tot = []
Nstheta_tot = np.array([[]]*Ntheta)
Ngtheta_tot = np.array([[]]*Ntheta)

# These lists will contain the density information for printing to a text file
field_area = [] # Will contain the masked area of the field (in arcmin)
field_galaxies = [] # Will contain the number of selected LRGs

# Looping over each KiDS field
for field in range(len(kidsfields)):
#for field in np.arange(1):
    
    print()
    print('KiDS field %s:'%kidsfields[field])

    # Boundaries of the current KiDS field
    fieldRAs = kidsboundaries[field,0]
    fieldDECs = kidsboundaries[field,1]
    print(fieldRAs, fieldDECs)
    
    # Selecting the galaxies lying within this field
    fieldmask = (fieldRAs[0] < galRA)&(galRA < fieldRAs[1]) & (fieldDECs[0] < galDEC)&(galDEC < fieldDECs[1])
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
    
    # Define the masked grid coordinates in this field
    gridRA, gridDEC, gridcoords, Nmaskedgrid = utils.define_gridpoints(fieldRAs, fieldDECs, galcoords, gridspace, 20./60.)
    Ngrid = len(gridcoords)
    gridID = np.arange(Ngrid)

    # Field area information for printing to text file
    field_area = np.append(field_area, Nmaskedgrid)
    
    ### 3) Measuring galaxy density:
    
    # These lists will contain the galaxy/grid point density around each grid point (for this field)
    Nstheta = np.zeros([Ntheta, Ngrid])
    Ngtheta = np.zeros([Ntheta, Ngrid])
   
    # Creating smaller grid point samples (of 10.000 points) to avoid memory overload
    sampindex = np.array(np.append(np.arange(0., Ngrid, 1e4), Ngrid), dtype=int)
    
    # Looping over different circle radius sizes (theta)
    for theta in range(Ntheta):
    #for theta in range(1):
        print('Theta =', thetalist[theta]*60., 'arcmin')

        # Looping over the smaller samples of grid points (to avoid memory overload)
        for s in range(len(sampindex)-1):
            
            print(kidsfields[field], ', Theta=' , thetalist[theta]*60., ', Grid sample: %i of %i'%(s+1, len(sampindex)-1))
            
            # Defining the smaller grid point sample
            gridsamp = gridcoords[sampindex[s]:sampindex[s+1]]
        
            # 3a) For each grid point, we count the number of sample galaxies within a circle of chosen radius \theta.
            galxgrid, idxcatalog, d2d, d3d = galcoords.search_around_sky(gridsamp, thetalist[theta]*u.deg)
            
            # 3b) For each grid point, we count the number of other grid points to estimate the effective area.
            gridxgrid, idxcatalog, d2d, d3d = gridcoords.search_around_sky(gridsamp, thetalist[theta]*u.deg)
            galcountlist = np.array(list(Counter(galxgrid).items())).T
            gridcountlist = np.array(list(Counter(gridxgrid).items())).T
            
            # Add the sample information to the full list for this field
            (Nstheta[theta, sampindex[s]:sampindex[s+1]])[galcountlist[0]] = galcountlist[1]
            (Ngtheta[theta, sampindex[s]:sampindex[s+1]])[gridcountlist[0]] = gridcountlist[1]
    
    # Combine the field information into one final catalog
    gridRA_tot = np.append(gridRA_tot, gridRA)
    gridDEC_tot = np.append(gridDEC_tot, gridDEC)
    Nstheta_tot = np.hstack([Nstheta_tot, Nstheta])
    Ngtheta_tot = np.hstack([Ngtheta_tot, Ngtheta])
    
# Define the total grid (of all fields)
Ngrid_tot = len(gridRA_tot)
gridID_tot = np.arange(Ngrid_tot)
gridcoords_tot = SkyCoord(ra=gridRA_tot*u.deg, dec=gridDEC_tot*u.deg)


#  3c) For each grid point, the galaxy density is defined as the galaxy count within the circle, divided by the effective area of the circle.

# Density
rhotheta_tot = Nstheta_tot/Ngtheta_tot
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
#for theta in range(1):

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
filename = '/data2/brouwer/MergedCatalogues/trough_catalog_test_%s.fits'%(selection)

# For each grid point/circle, we save the following information to the catalog:
# the location (RA/DEC), number count (Nstheta), effective area in arcmin (grid count Ngtheta), galaxy density (rhotheta), 
# percentage of circles below its density (Ptheta), and non-overlapping selection flag (Stheta).

# Defining the names of the columns
outputnames = ['RA', 'DEC']
[outputnames.append('Nstheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Ngtheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('rhotheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Ptheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Stheta%i'%(theta*60)) for theta in thetalist]

# Defining the output
output = [gridRA_tot, gridDEC_tot]
[output.append(Nstheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(Ngtheta_tot[theta,:]) for theta in range(Ntheta)]
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
