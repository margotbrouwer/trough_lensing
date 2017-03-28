#!/usr/bin/python

"Module to determine the local overdensity delta_r within a sphere of radius r."

# Import the necessary libraries
import astropy.io.fits as pyfits
import numpy as np

import trough_modules_all as utils

path_cat = '/data2/brouwer/MergedCatalogues'
catname = 'KIDS450_216.0_-1.5_r_sci_masked_pixels.fits'

# Full directory & name of the corresponding KiDS catalogue
catfile = '%s/%s'%(path_cat, catname)
cat = pyfits.open(catfile, memmap=True)

# List of the observables of all sources in the KiDS catalogue

matrix = np.array(cat[1].data)
RAs = (cat[2].data)
DECs = (cat[3].data)

masklist =  np.array([[ matrix[i,j] for i in range(len(RAs)) ] for j in range(len(DECs)) ])
coordlist = np.array([[ [RAs[i], DECs[j]] for i in range(len(RAs)) ] for j in range(len(DECs)) ])

masklist = np.reshape(masklist, [len(RAs)*len(DECs)])
RAlist = np.reshape(coordlist[:,:,0], [len(RAs)*len(DECs)])
DEClist = np.reshape(coordlist[:,:,1], [len(RAs)*len(DECs)])

print(RAlist)
print(DEClist)

filename = 'KiDS-450_mask.fits'
outputnames = ['RA', 'DEC', 'mask']
output = [RAlist, DEClist, masklist]


print('Writing output catalogue...')
utils.write_catalog(filename, np.arange(len(RAlist)), outputnames, output)
