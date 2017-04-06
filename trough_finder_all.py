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
from astropy.cosmology import LambdaCDM

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

import trough_modules_all as utils


# Radii theta of circular regions (in deg)
#thetalist = np.array([5., 10., 15., 20.])/60.
thetalist = np.array([2.764, 5., 5.527, 8.291, 10., 11.054, 15., 20.])/60.
#thetalist = np.array([5.])/60.

Ntheta = len(thetalist)


### 1) Defining the galaxy sample:

## 1a) Importing the galaxy catalogue.

# Select the galaxy catalogue for trough selection (kids/gama)
cat = 'gama'
masktype = 'complex'

# Spacing of the trough and mask grids (in degree)
gridspace = 0.04
mask_density = 1/(gridspace*60.)**2 # Density of mask gridpoints (in arcmin^-2)

# Import maskfile if present
maskfilename = '/data2/brouwer/MergedCatalogues/Masks/mask_catalogue_%s_%gdeg_%s.fits'%(cat, gridspace, masktype)
masktextname = 'area_info_%s.txt'%masktype
if os.path.isfile(maskfilename):
    nomaskfile = False
else:
    nomaskfile = True


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
    
    """
    # Names of the GAMA fields
    fieldnames = ['G9', 'G12', 'G15']
    
    # Boundaries of the GAMA fields
    coordsG9 = [[129., 141.], [-2.,3.]]
    coordsG12 = [[174., 186.], [-3.,2.]]
    coordsG15 = [[211.5, 223.5], [-2.,3.]]
    fieldboundaries = np.array([coordsG9,coordsG12,coordsG15])
    """
    
    # Path to the KiDS fieldsc
    path_kidscat = '/data2/brouwer/MergedCatalogues'
    kidscatname = 'KiDS_DR3_GAMA-like_290317.fits'
    
    # Importing the KiDS coordinates
    galRA, galDEC, galZ, galTB, mag_auto, ODDS, umag, gmag, rmag, imag = \
    utils.import_kidscat(path_kidscat, kidscatname)


if cat == 'gama':
    
    # Names of the GAMA fields
    fieldnames = ['G9', 'G12', 'G15']
    
    # Boundaries of the GAMA fields
    coordsG9 = [[129., 141.], [-2.,3.]]
    coordsG12 = [[174., 186.], [-3.,2.]]
    coordsG15 = [[211.5, 223.5], [-2.,3.]]
    fieldboundaries = np.array([coordsG9,coordsG12,coordsG15])

    # Path to the GAMA fields
    path_gamacat = '/data2/brouwer/MergedCatalogues/'
    gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'
    
    # Importing the GAMA coordinates
    galRA, galDEC, galZ, rmag, rmag_abs = utils.import_gamacat(path_gamacat, gamacatname)
    

## 1b) Select the galaxy sample to define the troughs

# Name of the pre-defined galaxy selection [all, ell, redseq, redseq4]
#selection = 'all'
selection = 'absmag'
#selection = 'redseq4'
#selection = 'lowZ'
#selection = 'highZ'

# Defining the selection for the KiDS galaxy sample
if cat=='kids':

    if selection == 'all':
        galmask = (galZ >= 0.)
        
    if selection == 'absmag':
        cosmo = LambdaCDM(H0=70., Om0=0.315, Ode0=0.685)
        galDc = (cosmo.comoving_distance(galZ).to('pc')).value
        rmag_abs = rmag - 5.*np.log10(galDc) + 5.
        selection = (rmag_abs < -19.7)

    if ('redseq' in selection) or (selection=='ell'):
        galmask = utils.define_galsamp(selection, zmin, zmax, galZ, galTB, gmag, rmag, mag_auto)
    
# Defining the selection for the GAMA galaxy sample
if cat=='gama':

    if selection == 'all':
        galmask = (rmag <= 19.8)
    if selection == 'absmag':
        galmask = (rmag_abs < -19.7) & (rmag <= 19.8)

    # Redshift samples
    zmin = 0.05
    zlim = 0.17075622469594484
    zmax = 0.3
    
    if selection == 'lowZ':
        galmask = (zmin < galZ)&(galZ < zlim) & (rmag_abs < -21.) & (rmag <= 19.8)
        thetalist = np.array([10.])/60.
        Ntheta = len(thetalist)

    if selection == 'highZ':
        galmask = (zlim < galZ)&(galZ < zmax) & (rmag_abs < -21.) & (rmag <= 19.8)
        thetalist = np.array([5.527])/60.
        Ntheta = len(thetalist)


print('Theta:', thetalist*60., 'arcmin')
print()
print( 'Total galaxy selection: %i/%i = %g percent'%(np.sum(galmask), len(galmask), float(np.sum(galmask))/float(len(galmask))*100.) )

# These lists will contain the full trough catalog (for all fields)
gridRA_tot = []
gridDEC_tot = []
Ngaltheta_tot = np.array([[]]*Ntheta)
Nmasktheta_tot = np.array([[]]*Ntheta)

# These lists will contain the density information for printing to a text file
field_area = [] # Will contain the effective area of the field (in arcmin)
field_galaxies = [] # Will contain the number of selected LRGs

# Looping over each field
for field in range(len(fieldnames)):

    # Boundaries of the current KiDS field
    fieldRAs = fieldboundaries[field,0]
    fieldDECs = fieldboundaries[field,1]
    
    # Selecting the galaxies lying within this field
    fieldmask = (fieldRAs[0] < galRA)&(galRA < fieldRAs[1]) & (fieldDECs[0] < galDEC)&(galDEC < fieldDECs[1])
    
    print()
    print('%s field %s:'%(cat, fieldnames[field]))

    # Applying the galaxy mask to the KiDS sample
    fieldmask = galmask*fieldmask

    # Coordinates of the galaxies
    galRA_field = galRA[fieldmask]
    galDEC_field = galDEC[fieldmask]
    galcoords = SkyCoord(ra=galRA_field*u.deg, dec=galDEC_field*u.deg)
    
    print('Selected: %g galaxies'%np.sum(fieldmask))
    field_galaxies = np.append(field_galaxies, np.sum(fieldmask))

    ### 2) Creating the grid:
    
    ## 2a) Creating a cartesian grid of narrowly spaced (2 arcmin) points.
   
    # Define the grid coordinates in this field

    # Import the mask coordinates of this field
    path_catmask = '/data2/brouwer/MergedCatalogues/Masks/%s_mask_%s_%gdeg.fits'%(cat, fieldnames[field], gridspace)
    
    maskcat = pyfits.open(path_catmask, memmap=True)[1].data
    gridRA, gridDEC, catmask = [maskcat['RA'], maskcat['DEC'], maskcat['mask']]
    gridcoords = SkyCoord(ra=gridRA*u.deg, dec=gridDEC*u.deg)
    Ngrid = len(gridcoords)
    
    #n, bins, patches = plt.hist(catmask, 50)
    #plt.show()
    
    # Count only the mask coordinates above 0.
    #gridcoords = gridcoords[catmask>0]
    #catmask = catmask[catmask>0]

    # Field area information for printing to text file
    Nmaskedgrid = np.sum(catmask)
    field_area = np.append(field_area, Nmaskedgrid/mask_density)

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
            
            if nomaskfile:
                # 3b) For each grid point, we count the number of other non-zero mask-points to estimate the effective area.
                maskxgrid, gridxmask, d2d, d3d = gridcoords.search_around_sky(gridsamp, thetalist[theta]*u.deg)
                print('Calculating the mask area...')
                
                if masktype == 'simple':
                    # Simple masking: completeness is either 1 or 0.
                    maskcountlist = np.array(list(Counter(maskxgrid).items())).T
                    (Nmasktheta[theta, sampindex[s]:sampindex[s+1]])[maskcountlist[0]] = maskcountlist[1]
                else:
                    # Improved masking: completeness varies between 0 and 1.
                    maskcountlist = np.array([np.sum(catmask[gridxmask][maskxgrid==g]) for g in range(len(gridsamp))])
                    Nmasktheta[theta, sampindex[s]:sampindex[s+1]] = maskcountlist

    
    # Combine the field information into one final catalog
    gridRA_tot = np.append(gridRA_tot, gridRA)
    gridDEC_tot = np.append(gridDEC_tot, gridDEC)
    Ngaltheta_tot = np.hstack([Ngaltheta_tot, Ngaltheta])
    Nmasktheta_tot = np.hstack([Nmasktheta_tot, Nmasktheta])

Ngrid_tot = len(gridRA_tot)
gridID_tot = np.arange(Ngrid_tot)

# Calculate the percentage of the circle areas that is masked

# If there is no file containing the masked percentage ...
if nomaskfile:
    # Define the percentage of effective area
    Pmasktheta_tot = Nmasktheta_tot/np.reshape(np.amax(Nmasktheta_tot, axis=1), [len(thetalist), 1])
    
    # Print this to a catalogue for later use
    outputnames = ['RA', 'DEC']
    output = [gridRA_tot, gridDEC_tot]

    [outputnames.append('Nmasktheta%i'%(theta*60)) for theta in thetalist]
    [output.append(Nmasktheta_tot[theta,:]) for theta in range(Ntheta)]


    [outputnames.append('Pmasktheta%i'%(theta*60)) for theta in thetalist]
    [output.append(Pmasktheta_tot[theta,:]) for theta in range(Ntheta)]

    utils.write_catalog(maskfilename, gridID_tot, outputnames, output)
    
    # Save the effective survey area into a text file for later use
    np.savetxt(masktextname, field_area, header = 'Total effective field area (in arcmin^2)')
    print('Written:', masktextname)

else:
    # Import the masked percentage
    maskcat = pyfits.open(maskfilename, memmap=True)[1].data
    
    Nmasktheta_tot= np.array([maskcat['Nmasktheta%i'%(theta*60)] for theta in thetalist])
    Pmasktheta_tot= np.array([maskcat['Pmasktheta%i'%(theta*60)] for theta in thetalist])
    
    field_area = np.loadtxt(masktextname)

# Effective area and density of the fields, for printing to text file
field_area = np.append(field_area, np.mean(field_area))
field_galaxies = np.append(field_galaxies, np.mean(field_galaxies))
field_density = field_galaxies/field_area

# Define the total grid (of all fields)
Ngrid_tot = len(gridRA_tot)
gridID_tot = np.arange(Ngrid_tot)
gridcoords_tot = SkyCoord(ra=gridRA_tot*u.deg, dec=gridDEC_tot*u.deg)

#  3c) For each grid point, the galaxy density is defined as the galaxy count within the circle, divided by the effective area of the circle.

### 5) Selecting the trough sample:

##  5a) For each grid point/circle, we calculate the percentage Ptheta of grid points that has a lower galaxy density.

# Density
Aefftheta_tot = Nmasktheta_tot/mask_density # Effective area of each circle
rhotheta_tot = Ngaltheta_tot/Aefftheta_tot

# Percentile and delta = (rho - av(rho))/av(rho)
sort_rho = np.array([np.argsort(rhotheta_tot[theta]) for theta in range(Ntheta)])
Ptheta_tot = np.array([np.argsort(gridID_tot[sort_rho[theta]])/Ngrid_tot for theta in range(Ntheta)])
delta_tot = (rhotheta_tot - field_density[-1]) / field_density[-1]

# Average redshift of the trough sample
redshift_av = np.mean(galZ[galmask])
redshift_tot = np.array([redshift_av]*Ngrid_tot)


# Following the definition from DES, we usually define the 20% of all circles with the lowest density to be the troughs.
# The 20% of all circles with the highest density are defined as overdensities ("ridges").


## Write catalog

## 5b) Now we have two trough samples:
# - The overlapping sample, by taking all troughs
# - The non-overlapping sample, by taking only unflagged troughs (not flagged = 1, flagged = 0).

# Writing the combined columns to a fits file
filename = '/data2/brouwer/MergedCatalogues/trough_catalogues/trough_catalog_%s_%s_%gdeg_%s.fits'%(cat, selection, gridspace, masktype)

# For each grid point/circle, we save the following information to the catalog:
# the location (RA/DEC), number count (Ngaltheta), effective area in arcmin (grid count Ngtheta), galaxy density (rhotheta), 
# percentage of circles below its density (Ptheta), and non-overlapping selection flag (Stheta).

# Defining the names of the columns
outputnames = ['RA', 'DEC', 'Z']
[outputnames.append('Ngaltheta%g'%(theta*60)) for theta in thetalist]
[outputnames.append('Aefftheta%g'%(theta*60)) for theta in thetalist]
[outputnames.append('Pmasktheta%g'%(theta*60)) for theta in thetalist]
[outputnames.append('rhotheta%g'%(theta*60)) for theta in thetalist]
[outputnames.append('Ptheta%g'%(theta*60)) for theta in thetalist]
[outputnames.append('delta%g'%(theta*60)) for theta in thetalist]
#[outputnames.append('Stheta%g'%(theta*60)) for theta in thetalist]

# Defining the output
output = [gridRA_tot, gridDEC_tot, redshift_tot]
[output.append(Ngaltheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(Aefftheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(Pmasktheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(rhotheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(Ptheta_tot[theta,:]) for theta in range(Ntheta)]
[output.append(delta_tot[theta,:]) for theta in range(Ntheta)]
#[output.append(Stheta_tot[theta,:]) for theta in range(Ntheta)]


print('Writing output catalogue...')
utils.write_catalog(filename, gridID_tot, outputnames, output)

# Writing the mean galaxy density information to a text file

field_header = '     '.join(np.append(fieldnames, 'Average'))
field_footer = '1: Number of galaxies, 2: Effective area (arcmin^2), 3: Galaxy density (arcmin^-2)'
density_info = np.array([field_galaxies, field_area, field_density])

filename = 'density_info_%s.txt'%selection

np.savetxt(filename, density_info, delimiter='    ', header = field_header, footer = field_footer)
print('Written:', filename)
