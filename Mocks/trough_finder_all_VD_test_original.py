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

sys.path.append('/disk1/vgd/ZOBOV')
sys.path.append('/disk1/vgd')
from read_vol_zone_void import *

from astropy import constants as const, units as u
from astropy.coordinates import SkyCoord
from collections import Counter

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

import trough_modules_all as utils

args = sys.argv
LOS = [int(val) for val in args[1:]][0]
print "LOS:",LOS
z_upper = 525

# Radii theta of circular regions (in deg)
#thetalist = np.array([5., 10., 15., 20.])/60.
thetalist = np.array([5.])/60.
Ntheta = len(thetalist)

print('Theta:', thetalist*60., 'arcmin')
print()

numpart_0_042, x_gal_0_042,y_gal_0_042,z_gal_0_042, x_halo_0_042,y_halo_0_042,z_halo_0_042, m200c_0_042, r200c_0_042, Rs_0_042, c_0_042 = read_HOD('../../HOD_Mocks/LOS%i/L505Mpc_HOD+0.042_LOS%i.dat' % (LOS,LOS))
numpart_0_130, x_gal_0_130,y_gal_0_130,z_gal_0_130, x_halo_0_130,y_halo_0_130,z_halo_0_130, m200c_0_130, r200c_0_130, Rs_0_130, c_0_130 = read_HOD('../../HOD_Mocks/LOS%i/L505Mpc_HOD+0.130_LOS%i.dat' % (LOS,LOS))
numpart_0_221, x_gal_0_221,y_gal_0_221,z_gal_0_221, x_halo_0_221,y_halo_0_221,z_halo_0_221, m200c_0_221, r200c_0_221, Rs_0_221, c_0_221 = read_HOD('../../HOD_Mocks/LOS%i/L505Mpc_HOD+0.221_LOS%i.dat' % (LOS,LOS))
numpart_0_317, x_gal_0_317,y_gal_0_317,z_gal_0_317, x_halo_0_317,y_halo_0_317,z_halo_0_317, m200c_0_317, r200c_0_317, Rs_0_317, c_0_317 = read_HOD('../../HOD_Mocks/LOS%i/L505Mpc_HOD+0.317_LOS%i.dat' % (LOS,LOS))
numpart_0_418, x_gal_0_418,y_gal_0_418,z_gal_0_418, x_halo_0_418,y_halo_0_418,z_halo_0_418, m200c_0_418, r200c_0_418, Rs_0_418, c_0_418 = read_HOD('../../HOD_Mocks/LOS%i/L505Mpc_HOD+0.418_LOS%i.dat' % (LOS,LOS))
numpart_0_525, x_gal_0_525,y_gal_0_525,z_gal_0_525, x_halo_0_525,y_halo_0_525,z_halo_0_525, m200c_0_525, r200c_0_525, Rs_0_525, c_0_525 = read_HOD('../../HOD_Mocks/LOS%i/L505Mpc_HOD+0.525_LOS%i.dat' % (LOS,LOS))


x_gal_tot = hstack((x_gal_0_042,x_gal_0_130,x_gal_0_221,x_gal_0_317,x_gal_0_418,x_gal_0_525))
y_gal_tot = hstack((y_gal_0_042,y_gal_0_130,y_gal_0_221,y_gal_0_317,y_gal_0_418,y_gal_0_525))
z_gal_tot = hstack((z_gal_0_042,z_gal_0_130,z_gal_0_221,z_gal_0_317,z_gal_0_418,z_gal_0_525))

# Convert physical coordinates to degree coordinates given box size in degrees to be used as an effective RA/DEC
# x_gal_deg is my effective RA
# y_gal_deg is my effective DEC
box_deg = 10.

x_gal_deg = x_gal_tot*(box_deg/(2.*max(z_gal_tot)*tan(radians(box_deg/2.))))
x_gal_deg = x_gal_deg+abs(min(x_gal_deg)) # This is to not have any neg values in my RA
y_gal_deg = y_gal_tot*(box_deg/(2.*max(z_gal_tot)*tan(radians(box_deg/2.))))
print 'Number of RA/DEC is:', len(x_gal_deg),'/',len(y_gal_deg)


galRA = x_gal_deg
galDEC = y_gal_deg

f = open('Mock_upto_z_%i_LOS%i.txt' % (z_upper,LOS), 'w')
f.write('#x [Mpc/h]  \y [Mpc/h]  z [Mpc/h]  RA  DEC\n')
for i in xrange(len(x_gal_tot)):
    f.write('{} {} {} {} {}\n'.format(x_gal_tot[i], y_gal_tot[i], z_gal_tot[i], galRA[i], galDEC[i]))
f.close()
raise()

### 1) Defining the galaxy sample:

## 1a) Importing the galaxy catalogue.

# Names of the KiDS fields
# kidsfields = ['G9', 'G12', 'G15', 'G23', 'GS']

# Boundaries of the KiDS fields
# coordsG9 = np.array([[128.0,142.5], [-2.5,3.5]])
# coordsG12 = np.array([[155.0,190.0], [-3.5,3.5]])
# coordsG15 = np.array([[209.5,239.0], [-3.5,3.5]])
# coordsG23 = np.array([[328.0,361.0], [-35.0,-28.0]])
# coordsGS = np.array([[31.0,54.0], [-35.0,-29.0]])
# kidsboundaries = np.array([coordsG9,coordsG12,coordsG15,coordsG23,coordsGS]) # Boundaries of all fields
mock_boundary = np.array([[min(x_gal_deg),max(x_gal_deg)],[min(y_gal_deg),max(y_gal_deg)]]) # Boundaries of all fields

# Path to the KiDS fields
# path_kidscat = '/data2/brouwer/KidsCatalogues'
# kidscatname = 'KiDS.DR3.zpbpzLess0.6.fits'

# Importing the KiDS coordinates
# galRA, galDEC, galZB, galTB, mag_auto, ODDS, umag, gmag, rmag, imag = \
# utils.import_kidscat(path_kidscat, kidscatname)

## 1b) Select the galaxy sample to define the troughs

# Name of the pre-defined galaxy selection [all, ell, redseq, redseq4]
selection = 'all'
#selection = 'redseq4'

# Redshift cut
# zmin = 0.1
# zmax = 0.4

### This is questionable since I don't know how to go from z-coord in Mpc/h to z (redshift)
zmin = 0.0
zmax = 0.525
print 'Upper z is: %0.3f' % zmax

# Defining the mask for the galaxy sample
# galmask = utils.define_galsamp(selection, zmin, zmax, galZB, galTB, gmag, rmag, mag_auto)


# These lists will contain the full trough catalog (for all fields)
gridRA_tot = []
gridDEC_tot = []
Nstheta_tot = np.array([[]]*Ntheta)
Ngtheta_tot = np.array([[]]*Ntheta)

# These lists will contain the density information for printing to a text file
field_area = [] # Will contain the masked area of the field (in arcmin)
field_galaxies = [] # Will contain the number of selected LRGs

# Looping over each KiDS field
# for field in range(len(kidsfields)):
#for field in np.arange(1):
    
    # print()
    # print('KiDS field %s:'%kidsfields[field])

# Boundaries of the current MOCK field
fieldRAs = mock_boundary[0]
fieldDECs = mock_boundary[1]
print(fieldRAs, fieldDECs)

# Selecting the galaxies lying within this field
# fieldmask = (fieldRAs[0] < galRA)&(galRA < fieldRAs[1]) & (fieldDECs[0] < galDEC)&(galDEC < fieldDECs[1])
fieldmask = len(galRA) #galmask

# Coordinates of the galaxies
galRA_field = galRA
galDEC_field = galDEC
galcoords = SkyCoord(ra=galRA_field*u.deg, dec=galDEC_field*u.deg)
    
print('Selected: %g galaxies'%np.sum(fieldmask))
field_galaxies = np.append(field_galaxies, np.sum(fieldmask))

### 2) Creating the grid:

## 2a) Creating a cartesian grid of narrowly spaced (1 arcmin) points.

# Spacing of the grid (in deg)
gridspace = 1./60. # One arcmin

## 2b) We remove all points that lie within 1 arcmin of masked areas or the edge of the field.

# Define the masked grid coordinates in this field
gridRA, gridDEC, gridcoords = utils.define_gridpoints(fieldRAs, fieldDECs, galcoords, gridspace, 20./60.) #, Nmaskedgrid
Ngrid = len(gridcoords)
gridID = np.arange(Ngrid)

# Field area information for printing to text file
field_area = ((max(galRA)-min(galRA))*(max(galDEC)-min(galDEC))*60.**2.)# The area of the field on the sky that's not masked #np.append(field_area, Nmaskedgrid)
    
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
        
        # print(kidsfields[field], ', Theta=' , thetalist[theta]*60., ', Grid sample: %i of %i'%(s+1, len(sampindex)-1))
        
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
# sort_rho_crab = [np.append( sort_rho[theta, 0:int(Ngrid_tot/2)], \
#                     np.flip(sort_rho[theta, int(Ngrid_tot/2):Ngrid_tot], 0) ) for theta in range(Ntheta)]


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
filename = 'trough_catalog_test_upto_z_%i_LOS%i.fits'%(z_upper,LOS)

# For each grid point/circle, we save the following information to the catalog:
# the location (RA/DEC), number count (Nstheta), effective area in arcmin (grid count Ngtheta), galaxy density (rhotheta), 
# percentage of circles below its density (Ptheta), and non-overlapping selection flag (Stheta).
# Get rid of circles that are not 80% complete.  Find number of grid points in each circle because they're the same
# Area of circle should be the number of grid points in the circle.  So edge circles will not have as many

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

field_header = 'Mock      Average'
field_footer = '1: Number of red galaxies, 2: Area after masking (arcmin^2), 3: Galaxy density (arcmin^-2)'
density_info = np.array([field_galaxies, field_area, field_density])

filename = 'density_info_%s.txt'%selection

np.savetxt(filename, density_info, delimiter='    ',  fmt='%.6g', header = field_header, footer = field_footer)
print('Written:', filename)
