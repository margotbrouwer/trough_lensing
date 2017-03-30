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

import trough_modules_all as utils

path_cat = '/data2/brouwer/KidsCatalogues'
catnames = os.listdir('%s/kids_masks'%path_cat)

print(len(catnames))

RAlist_tot = []
DEClist_tot = []
masklist_tot = []

for catname in catnames:
    # Full directory & name of the corresponding KiDS catalogue
    catfile = '%s/kids_masks/%s'%(path_cat, catname)
    cat = pyfits.open(catfile, memmap=True)

    # List of the observables of all sources in the KiDS catalogue

    matrix = np.array(cat[1].data)
    RAs = (cat[2].data)
    DECs = (cat[3].data)

    masklist =  np.array([[ matrix[i,j] for i in range(len(RAs)) ] for j in range(len(DECs)) ])
    coordlist = np.array([[ [RAs[i], DECs[j]] for i in range(len(RAs)) ] for j in range(len(DECs)) ])

    RAlist = np.reshape(coordlist[:,:,0], [len(RAs)*len(DECs)])
    DEClist = np.reshape(coordlist[:,:,1], [len(RAs)*len(DECs)])
    masklist = np.reshape(masklist, [len(RAs)*len(DECs)])

    print(RAlist)
    print(DEClist)

    RAlist_tot = np.append(RAlist_tot, RAlist)
    DEClist_tot = np.append(DEClist_tot, DEClist)
    masklist_tot = np.append(masklist_tot, masklist)
    

filename = '%s/KiDS-450_mask.fits'%path_cat
outputnames = ['RA', 'DEC', 'mask']
output = [RAlist_tot, DEClist_tot, masklist_tot]


print('Writing output catalogue...')
utils.write_catalog(filename, np.arange(len(RAlist_tot)), outputnames, output)
