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

import trough_modules_all as utils

path_cat = '/data2/brouwer/KidsCatalogues'
catnames = os.listdir('%s/kids_masks'%path_cat)

print(len(catnames))

RAlist_tot = []
DEClist_tot = []
masklist_tot = []

gridspace_mask = 0.02 # in degree
gridspace_orig = 0.002
gapsize = gridspace_mask/gridspace_orig

# Names of the KiDS fields
fieldnames = ['G9', 'G12', 'G15', 'G23', 'GS']

# Boundaries of the KiDS fields
coordsG9 = np.array([[128.0,142.5], [-2.5,3.5]])
coordsG12 = np.array([[155.0,190.0], [-3.5,3.5]])
coordsG15 = np.array([[209.5,239.0], [-3.5,3.5]])
coordsG23 = np.array([[328.0,361.0], [-35.0,-28.0]])
coordsGS = np.array([[31.0,54.0], [-35.0,-29.0]])
fieldboundaries = np.array([coordsG9,coordsG12,coordsG15,coordsG23,coordsGS]) # Boundaries of all fields


for c in range(len(catnames)):
#for c in range(20):
    
    # Full directory & name of the corresponding KiDS catalogue
    catfile = '%s/kids_masks/%s'%(path_cat, catnames[c])
    cat = pyfits.open(catfile, memmap=True)

    kidsmask = 1.-np.array(cat[1].data)
    RAs = (cat[2].data)
    DECs = (cat[3].data)
  
    print('Importing KiDS mask: %s (%i/%i)'%(catnames[c], c+1, len(catnames)))
    #print('Old size:', np.shape(kidsmask))
   
    RAnums = np.arange(int(gapsize/2.), int(len(kidsmask)), int(gapsize))
    DECnums = np.arange(int(gapsize/2.), int(len(kidsmask[0])), int(gapsize))
    
    kidsmask_small = np.zeros([len(RAnums), len(DECnums)])

    for i in range(len(RAnums)):
        for j in range(len(DECnums)):
            maskmean = np.mean(kidsmask[int(i*gapsize):int((i+1)*gapsize), int(j*gapsize):int((j+1)*gapsize)])
            kidsmask_small[i, j] = maskmean
    
    #print('New size:', np.shape(kidsmask_small))
    #print()
    
    kidsmask = kidsmask_small
    RAs = RAs[RAnums]
    DECs = DECs[DECnums]
    
    masklist =  np.array([[ kidsmask[i,j] for i in range(len(RAs)) ] for j in range(len(DECs)) ])
    coordlist = np.array([[ [RAs[i], DECs[j]] for i in range(len(RAs)) ] for j in range(len(DECs)) ])

    RAlist = np.reshape(coordlist[:,:,0], [len(RAs)*len(DECs)])
    DEClist = np.reshape(coordlist[:,:,1], [len(RAs)*len(DECs)])
    masklist = np.reshape(masklist, [len(RAs)*len(DECs)])

    RAlist_tot = np.append(RAlist_tot, RAlist)
    DEClist_tot = np.append(DEClist_tot, DEClist)
    masklist_tot = np.append(masklist_tot, masklist)


for f in range(len(fieldnames)):

    # Boundaries of the current KiDS field
    fieldRAs = fieldboundaries[f,0]
    fieldDECs = fieldboundaries[f,1]
    
    # Selecting the mask points lying within this field
    fieldmask = (fieldRAs[0] < RAlist_tot)&(RAlist_tot < fieldRAs[1]) & (fieldDECs[0] < DEClist_tot)&(DEClist_tot < fieldDECs[1])
    Ngrid_field = np.sum(fieldmask)
    
    RAlist_field, DEClist_field, masklist_field = RAlist_tot[fieldmask], DEClist_tot[fieldmask], masklist_tot[fieldmask]
    
    maskcoords = SkyCoord(ra=RAlist_field*u.deg, dec=DEClist_field*u.deg)

    # This list will contain the selected (1) and flagged (0) mask points (for all fields)
    selected = np.ones(Ngrid_field)

    # Creating smaller grid point samples to avoid memory overload
    sampindex_field = np.array(np.append(np.arange(0., Ngrid_field, 1e4), Ngrid_field), dtype=int)

    # Looping over the smaller samples of grid points (to avoid memory overload)
    for s in range(len(sampindex_field)-1):
    #for s in range(1):

        print('Removing points within 0.5 gridspace of each other. Grid sample: %i of %i'%(s+1, len(sampindex_field)-1))
        
        # Defining the smaller grid point sample
        selsamp = maskcoords[sampindex_field[s]:sampindex_field[s+1]]
        
        # 4c/d) For the underdense/overdense circles: We start from the lowest/highest density circle and flag all overlapping circles.
        
        # Find grid points within 0.5 gridspace of each other
        sampxgrid, gridoverlap, d2d, d3d = maskcoords.search_around_sky(selsamp, (gridspace_mask/2.)*u.deg)
        
        # Looping over grid points
        for g in range(len(selsamp)):
            
            # Index of the grid point in the flag list
            field_index = sampindex_field[s]+g
            
            # If the point is not flagged (not flagged = 1), flag all overlapping grid points (flagged = 0).
            if selected[field_index] == 1:
                removeindex = gridoverlap[(sampxgrid==g)&(gridoverlap!=field_index)]
                selected[removeindex] = 0.
            
            # If the point is already flagged, do nothing.
            else:
                pass

    # Remove overlapping mask points
    selmask = (selected == 1)
    
    # Print output to file
    filename = '%s/KiDS-450_mask_%s.fits'%(path_cat, fieldnames[f])
    outputnames = ['RA', 'DEC', 'mask']
    output = [RAlist_field[selmask], DEClist_field[selmask], masklist_field[selmask]]


    print('Writing output catalogue...')
    utils.write_catalog(filename, np.arange(len(RAlist_field[selmask])), outputnames, output)
