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

sys.path.append('/data2/brouwer/shearprofile/trough_lensing')
import trough_modules_all as utils

from astropy import constants as const, units as u
from astropy.coordinates import SkyCoord
from collections import Counter

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams


# Radii theta of circular regions (in deg)
#thetalist = np.array([5., 10., 15., 20.])/60.
thetalist = np.array([5.])/60.
Ntheta = len(thetalist)

print('Theta:', thetalist*60., 'arcmin')
print()

### 1) Defining the galaxy sample:

## 1a) Importing the galaxy catalogue.

#x_mock, y_mock, z_mock, galRA, galDEC = np.loadtxt('Mock_upto_z_525_LOS500.txt').T
mockcatfile = '/data2/brouwer/MergedCatalogues/trough_catalogues/mock_galaxies.fits'
mockcat = pyfits.open(mockcatfile, memmap=True)[1].data
mockdatanames = ['x_mock', 'y_mock', 'z_mock', 'galRA', 'galDEC']
x_mock, y_mock, z_mock, galRA, galDEC = np.array([mockcat[name] for name in mockdatanames])

# Galaxy selection
galmask = (1400 < z_mock)
galRA, galDEC = [galRA[galmask], galDEC[galmask]]

# Boundaries of the field
mock_boundary = np.array([[np.amin(galRA), np.amax(galRA)], [np.amin(galDEC), np.amax(galDEC)]])
fieldRAs = mock_boundary[0]
fieldDECs = mock_boundary[1]

# Coordinates of the galaxies
galcoords = SkyCoord(ra=galRA*u.deg, dec=galDEC*u.deg)

# These lists will contain the density information for printing to a text file
field_galaxies = len(galcoords) # Will contain the number of selected galaxies
field_area = (mock_boundary[0,1]-mock_boundary[0,0]) * (mock_boundary[1,1]-mock_boundary[1,0]) * 60.**2.

print('Selected: %g galaxies'%field_galaxies)

### 2) Creating the grid:

## 2a) Creating a cartesian grid of narrowly spaced (1 arcmin) points.

# Spacing of the grid (in deg)
gridspace = 1./60. # One arcmin

# Define the grid coordinates in this field
gridRA, gridDEC, gridcoords = utils.define_gridpoints(fieldRAs, fieldDECs, galcoords, gridspace)
Ngrid = len(gridcoords)
gridID = np.arange(Ngrid)


### 3) Measuring galaxy density:

# These lists will contain the galaxy/grid point density around each grid point (for this field)
Ngaltheta = np.zeros([Ntheta, Ngrid])
Ngridtheta = np.zeros([Ntheta, Ngrid])


# Creating smaller grid point samples (of 10.000 points) to avoid memory overload
sampindex = np.array(np.append(np.arange(0., Ngrid, 1e4), Ngrid), dtype=int)

# Looping over different circle radius sizes (theta)
for theta in range(Ntheta):
#for theta in range(1):
    print('Theta =', thetalist[theta]*60., 'arcmin')

    # Looping over the smaller samples of grid points (to avoid memory overload)
    for s in range(len(sampindex)-1):
    #for s in range(1):
        
        print('Theta =' , thetalist[theta]*60., ', Grid sample: %i of %i'%(s+1, len(sampindex)-1))
        
        # Defining the smaller grid point sample
        gridsamp = gridcoords[sampindex[s]:sampindex[s+1]]
    
        # 3a) For each grid point, we count the number of sample galaxies within a circle of chosen radius \theta.
        galxgrid, idxcatalog, d2d, d3d = galcoords.search_around_sky(gridsamp, thetalist[theta]*u.deg)
        
        # 3b) For each grid point, we count the number of other grid points to estimate the effective area.
        gridxgrid, idxcatalog, d2d, d3d = gridcoords.search_around_sky(gridsamp, thetalist[theta]*u.deg)
        galcountlist = np.array(list(Counter(galxgrid).items())).T
        gridcountlist = np.array(list(Counter(gridxgrid).items())).T
        
        # Add the sample information to the full list for this field
        (Ngaltheta[theta, sampindex[s]:sampindex[s+1]])[galcountlist[0]] = galcountlist[1]
        (Ngridtheta[theta, sampindex[s]:sampindex[s+1]])[gridcountlist[0]] = gridcountlist[1]


#  3c) For each grid point, the galaxy density is defined as the galaxy count within the circle, divided by the effective area of the circle.

# Percentage Pgrid of the full area in each circle (Always remove gridpoints with Pgridtheta<0.8 from your sample)
Pgridtheta = Ngridtheta/np.reshape(np.amax(Ngridtheta, axis=1), [len(thetalist), 1])

# Density
rhotheta = Ngaltheta/Ngridtheta
rho_mean = np.mean(rhotheta, 1)
print('Mean density', rho_mean)

### 5) Selecting the trough sample:
    
##  5a) For each grid point/circle, we calculate the percentage Ptheta of grid points that has a lower galaxy density.
# Percentile
sort_rho = np.array([np.argsort(rhotheta[theta]) for theta in range(Ntheta)])
Ptheta = np.array([np.argsort(gridID[sort_rho[theta]])/Ngrid for theta in range(Ntheta)])

# Following the definition from DES, we usually define the 20% of all circles with the lowest density to be the troughs.
# The 20% of all circles with the highest density are defined as overdensities ("ridges").


## Write catalog


# Writing the combined columns to a fits file

# For each grid point/circle, we save the following information to the catalog:
# the location (RA/DEC), galaxy number count (Ngaltheta), effective area in arcmin (grid count Ngridtheta), galaxy density (rhotheta), 
# percentage of circles below its density (Ptheta), and non-overlapping selection flag (Stheta).
# Get rid of circles that are not 80% complete.  Find number of grid points in each circle because they're the same
# Area of circle should be the number of grid points in the circle.  So edge circles will not have as many

# Defining the names of the columns
outputnames = ['RA', 'DEC']
[outputnames.append('Ngaltheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Ngridtheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Pgridtheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('rhotheta%i'%(theta*60)) for theta in thetalist]
[outputnames.append('Ptheta%i'%(theta*60)) for theta in thetalist]

# Defining the output
output = [gridRA, gridDEC]
[output.append(Ngaltheta[theta,:]) for theta in range(Ntheta)]
[output.append(Ngridtheta[theta,:]) for theta in range(Ntheta)]
[output.append(Pgridtheta[theta,:]) for theta in range(Ntheta)]
[output.append(rhotheta[theta,:]) for theta in range(Ntheta)]
[output.append(Ptheta[theta,:]) for theta in range(Ntheta)]


print('Writing output catalogue...')
filename = '/data2/brouwer/MergedCatalogues/trough_catalogues/mock_trough_catalogue.fits'
utils.write_catalog(filename, gridID, outputnames, output)


# Writing the mean galaxy density information to a text file
field_density = field_galaxies/field_area

field_header = '1: Number of galaxies, 2: Effective area (arcmin^2), 3: Galaxy density (arcmin^-2)'
density_info = np.array([field_galaxies, field_area, field_density])

filename = 'density_info.txt'

np.savetxt(filename, density_info, delimiter='    ', header = field_header)
print('Written:', filename)
