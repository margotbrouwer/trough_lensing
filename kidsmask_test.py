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

cat = 'kids'
gridspace_mask = 0.04 # in degree


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
    
    # Names of the mask files
    kids_path = '/data2/brouwer/KidsCatalogues/kids_masks_new'
    catnames = os.listdir(kids_path)
    #catnames = ['KIDS450_0.0_-31.2_r_sci_masked_pixels.fits']
    
    # Path to the KiDS fields
    path_kidscat = '/data2/brouwer/MergedCatalogues'
    kidscatname = 'KiDS_DR3_GAMA-like_290317.fits'
    
    """
    # Importing the KiDS coordinates
    galRA, galDEC, galZ, galTB, mag_auto, ODDS, umag, gmag, rmag, imag = \
    utils.import_kidscat(path_kidscat, kidscatname)
    galcoords = SkyCoord(ra=galRA*u.deg, dec=galDEC*u.deg)
    """

if cat == 'gama':
    
    # Names of the GAMA fields
    fieldnames = ['G9', 'G12', 'G15']
    
    # Boundaries of the GAMA fields
    coordsG9 = [[129., 141.], [-2.,3.]]
    coordsG12 = [[174., 186.], [-3.,2.]]
    coordsG15 = [[211.5, 223.5], [-2.,3.]]
    
    fieldboundaries = np.array([coordsG9,coordsG12,coordsG15])
    
    # Names of the mask files
    catnames = ['/data2/brouwer/MergedCatalogues/GamaMasks/%smask08000.fits'%g for g in fieldnames]


## Making the big matrix files into a smaller array with [RA, DEC, mask] for every coordinate

RAlist_tot = []
DEClist_tot = []
masklist_tot = []

for c in range(len(catnames)):
#for c in range(1):
    
    print('Importing %s mask: %s (%i/%i)'%(cat, catnames[c], c+1, len(catnames)))
    
    # Import the mask catalogue
    if cat == 'kids':
        maskcat = pyfits.open('%s/%s'%(kids_path, catnames[c]), memmap=True)

        catmask = 1.-np.array(maskcat[1].data)
        maskRA, maskDEC = [maskcat[2].data, maskcat[3].data]

        # Finding the space between the original grid points
        RAdiff, DECdiff = [np.diff(maskRA), np.diff(maskDEC)]
        RAdiff, DECdiff = [RAdiff[RAdiff<1.], DECdiff[DECdiff<1.]]
        
        gridspace_RA = abs(np.mean(RAdiff))
        gridspace_DEC = abs(np.mean(DECdiff))
        gapsize_RA = round(gridspace_mask/gridspace_RA)
        gapsize_DEC = round(gridspace_mask/gridspace_DEC)
        
        #print(gridspace_RA, gapsize_RA)
        #print(gridspace_DEC, gapsize_DEC)
        
        """
        # Highlight fields without galaxies
        
        fieldcoord = SkyCoord(ra=np.mean(maskRA)*u.deg, dec=np.mean(maskDEC)*u.deg)
        idx, d2d, d3d = fieldcoord.match_to_catalog_sky(galcoords)
        d2d = d2d[0].to('degree').value
        
        print(d2d)
        if d2d > 0.1:
            print('No galaxies in:', catnames[c])
        """
        
    if cat == 'gama':
        # Boundaries of the current GAMA field
        fieldRAs = fieldboundaries[c,0]
        fieldDECs = fieldboundaries[c,1]
        
        catmask = np.array(pyfits.open(catnames[c], memmap=True)['PRIMARY'].data).T
        maskRA, maskDEC, maskcoords = utils.define_gridpoints(fieldRAs, fieldDECs, gridspace_mask)
        maskRA, maskDEC = [np.sort(np.unique(maskRA)), np.sort(np.unique(maskDEC))]
        
        # The space between the original grid points
        gridspace_orig = 0.001
        
        # Defining the gap between the points of the big matrix, to create the small matrix
        gapsize_RA = round(gridspace_mask/gridspace_orig)
        gapsize_DEC = gapsize_RA
    
    # Set all negative mask values to 0.
    catmask[catmask < 0.] = 0.

    RAnums = np.arange(int(round(gapsize_RA/2.)), np.shape(catmask)[0], int(gapsize_RA))
    DECnums = np.arange(int(round(gapsize_DEC/2.)), np.shape(catmask)[1], int(gapsize_DEC))
    
    catmask_small = np.zeros([len(RAnums), len(DECnums)])

    for i in range(len(RAnums)):
        for j in range(len(DECnums)):
            maskmean = np.mean(catmask[int(i*gapsize_RA):int((i+1)*gapsize_RA), int(j*gapsize_DEC):int((j+1)*gapsize_DEC)])
            catmask_small[i, j] = maskmean
    
    print('Old size:', np.shape(catmask))
    print('New size:', np.shape(catmask_small))
    
    catmask = catmask_small
    if cat == 'kids':
        maskRA = maskRA[RAnums]
        maskDEC = maskDEC[DECnums]
    
    masklist =  np.array([[ catmask[i,j] for i in range(len(maskRA)) ] for j in range(len(maskDEC)) ])
    coordlist = np.array([[ [maskRA[i], maskDEC[j]] for i in range(len(maskRA)) ] for j in range(len(maskDEC)) ])
    
    RAlist = np.reshape(coordlist[:,:,0], [len(maskRA)*len(maskDEC)])
    DEClist = np.reshape(coordlist[:,:,1], [len(maskRA)*len(maskDEC)])
    masklist = np.reshape(masklist, [len(maskRA)*len(maskDEC)])
    
    RAlist_tot = np.append(RAlist_tot, RAlist)
    DEClist_tot = np.append(DEClist_tot, DEClist)
    masklist_tot = np.append(masklist_tot, masklist)
    
    #print()

## Dividing total masks into different fields, and removing overlapping tiles
for f in range(len(fieldnames)):

    # Boundaries of the current KiDS field
    fieldRAs = fieldboundaries[f,0]
    fieldDECs = fieldboundaries[f,1]
    
    # Selecting the mask points lying within this field
    fieldmask = (fieldRAs[0] <= RAlist_tot)&(RAlist_tot <= fieldRAs[1]) & (fieldDECs[0] <= DEClist_tot)&(DEClist_tot <= fieldDECs[1])
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

        print('Removing points within 0.5*sqrt(2) gridspace of each other. Grid sample: %i of %i'%(s+1, len(sampindex_field)-1))
        
        # Defining the smaller grid point sample
        selsamp = maskcoords[sampindex_field[s]:sampindex_field[s+1]]
        
        # 4c/d) For the underdense/overdense circles: We start from the lowest/highest density circle and flag all overlapping circles.
        
        # Find grid points within 0.8 gridspace of each other
        sampxgrid, gridoverlap, d2d, d3d = maskcoords.search_around_sky(selsamp, (np.sqrt(2.)/2.)*gridspace_mask*u.deg)
        
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
    print('Selected:', np.sum(selmask), '/', len(selmask))
    
    # Print output to file
    
    filename = '/data2/brouwer/MergedCatalogues/Masks/%s_mask_%s_%gdeg.fits'%(cat, fieldnames[f], gridspace_mask)
    outputnames = ['RA', 'DEC', 'mask']
    output = [RAlist_field[selmask], DEClist_field[selmask], masklist_field[selmask]]


    print('Writing output catalogue...')
    utils.write_catalog(filename, np.arange(len(RAlist_field[selmask])), outputnames, output)
