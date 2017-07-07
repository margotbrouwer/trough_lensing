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

h, O_matter, O_lambda = [0.7, 0.25, 0.75]
cosmo = LambdaCDM(H0=h*100, Om0=O_matter, Ode0=O_lambda)
micecor = 5*np.log10(h) # Correction on MICE absmag_r (= -0.7745)


ijlist = np.array([ [ [i, j] for i in range(4) ] for j in range(4) ])
ijlist = np.reshape(ijlist, [16,2])

#for ij in np.arange(0, len(ijlist)):
for ij in range(1):

    ### 1) Defining the galaxy sample:

    cat = 'mice'

    ## 1a) Importing the galaxy catalogue.

    # Name of the pre-defined galaxy selection
    selection = 'miceZ'
    
    
    # Select mask type (nomask or complex)
    i, j = ijlist[ij]
    ijnum = i + 4.*j + 1.
    print(i, j, ijnum)
    
    #masktype = 'nomask-%g'%(ijnum)
    masktype = 'nomask-new'
    
    # Spacing of the trough and mask grids (in degree)
    gridspace = 0.04
    mask_density = 1/(gridspace*60.)**2 # Density of mask gridpoints (in arcmin^-2)
    
    
    # Import galaxy catalog
    
    # Names of the MICE fields
    fieldnames = ['M']
    
    # Boundaries of the MICE field
    coordsM = [[i*20.,(i+1)*20.], [j*20.,(j+1)*20.]]
    fieldboundaries = np.array([coordsM]) # Boundaries of all fields
    
    print(coordsM)
    
    # Path to the Mice field
    path_mockcat = '/data2/brouwer/MergedCatalogues'
    mockcatname = 'mice_gama_catalog.fits'
    
    # Importing the Mice galaxies
    galRA, galDEC, galZ, galDc, rmag, rmag_abs, e1, e2 = \
    utils.import_mockcat(path_mockcat, mockcatname)
    rmag_abs = rmag_abs + micecor
    
    gama_rlim = 20.2


    # Import maskfile if present
    maskfilename = '/data2/brouwer/MergedCatalogues/Masks/mask_catalog_%s_%gdeg_%s.fits'%(cat, gridspace, masktype)
    masktextname = 'density_info/area_info_%s_%s.txt'%(cat, masktype)
    if os.path.isfile(maskfilename):
        nomaskfile = False
        print('Importing mask file:', maskfilename)
    else:
        nomaskfile = True
        print('No mask file present')
    

    ## 1b) Select the galaxy sample to define the troughs

    # Defining the selection for the galaxy sample
        
    if 'miceZ' in selection:
        
        #zlims = np.array([0.1, 0.198, 0.3])
        #thetalist = np.array([10., 6.826])/60.
        
        zlims = np.array([ 0.1, 0.193, 0.290, 0.392, 0.5])
        thetalist = np.array([10., 6.892, 5.447, 4.616])/60.
        Ntheta = len(thetalist)
        
        galmask = [ (zlims[theta] < galZ)&(galZ < zlims[theta+1]) & (rmag_abs < -21.) & (rmag < gama_rlim) \
                        for theta in range(Ntheta) ]


    Ngalsel = np.array([ float(sum(galmask[theta])) for theta in range(Ntheta) ])
    Ngaltot = np.array([ float(len(galmask[theta])) for theta in range(Ntheta) ])
    
    print('Theta:', thetalist*60., 'arcmin')
    print()
    print('Total galaxy selection:', Ngalsel, Ngaltot, Ngalsel/Ngaltot*100., 'percent')

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
        field_galaxies = np.append(field_galaxies, np.sum(fieldmask))
        
        print()
        print('%s field %s: contains %g galaxies'%(cat, fieldnames[field], np.sum(fieldmask)))
        
        ### 2) Creating the grid:
        
        ## 2a) Creating a cartesian grid of narrowly spaced (2 arcmin) points.
       
        # Define the grid coordinates in this field
        
        # Creating the grid for each field
        gridspace_mask = 0.04 # in degree
        gridRA, gridDEC, gridcoords = utils.define_gridpoints(fieldRAs, fieldDECs, gridspace_mask, True)

        Ngrid = len(gridcoords)    
        catmask = np.ones(Ngrid)
        
        
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
            
            # Applying the galaxy mask to the KiDS sample
            sampmask = galmask[theta] * fieldmask
            
            # Coordinates of the galaxies
            galRA_field = galRA[sampmask]
            galDEC_field = galDEC[sampmask]
            galcoords = SkyCoord(ra=galRA_field*u.deg, dec=galDEC_field*u.deg)
            
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
                    
                    if 'nomask' in masktype:
                        # No masking: completeness is always 1.
                        maskcountlist = np.array(list(Counter(maskxgrid).items())).T
                        (Nmasktheta[theta, sampindex[s]:sampindex[s+1]])[maskcountlist[0]] = maskcountlist[1]
                    if 'complex' in masktype:
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
        Nmaxtheta = np.pi * (thetalist * 60.)**2. * mask_density
        #Nmaxtheta = np.amax(Nmasktheta_tot, axis=1)
        
        Pmasktheta_tot = Nmasktheta_tot/np.reshape(Nmaxtheta, [len(thetalist), 1])
        
        # Print this to a catalogue for later use
        outputnames = ['RA', 'DEC']
        output = [gridRA_tot, gridDEC_tot]

        [outputnames.append('Nmasktheta%g'%(theta*60)) for theta in thetalist]
        [output.append(Nmasktheta_tot[theta,:]) for theta in range(Ntheta)]


        [outputnames.append('Pmasktheta%g'%(theta*60)) for theta in thetalist]
        [output.append(Pmasktheta_tot[theta,:]) for theta in range(Ntheta)]

        utils.write_catalog(maskfilename, gridID_tot, outputnames, output)
        
        # Save the effective survey area into a text file for later use
        np.savetxt(masktextname, field_area, header = 'Total effective field area (in arcmin^2)')
        print('Written:', masktextname)
        
    else:
        # Import the masked percentage
        print('Importing mask from:', maskfilename)
        maskcat = pyfits.open(maskfilename, memmap=True)[1].data
        
        Nmasktheta_tot= np.array([(maskcat['Nmasktheta%g'%(theta*60)])[0:Ngrid_tot] for theta in thetalist])
        Pmasktheta_tot= np.array([(maskcat['Pmasktheta%g'%(theta*60)])[0:Ngrid_tot] for theta in thetalist])
        
        field_area = np.array(np.loadtxt(masktextname))
        try:
            field_area = field_area[0:len(fieldnames)]
        except:
            field_area = np.array([field_area])
        
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
    
    # Average redshift of the trough sample

    troughZ = np.array([ np.mean(galZ[(galmask[theta])]) for theta in range(Ntheta) ])
    troughZ_tot = np.array([ [troughZ[theta]]*Ngrid_tot for theta in range(Ntheta) ])
    
    troughDc = np.array([ np.mean(galDc[(galmask[theta])]) for theta in range(Ntheta) ])
    troughDa = np.array([ np.mean( galDc[(galmask[theta])] / (1 + galZ[(galmask[theta])]) ) for theta in range(Ntheta) ])
    
    am_to_rad = np.pi/(60.*180.)
    Aefftheta_tot = (Nmasktheta_tot/mask_density) * (am_to_rad * troughDa)**2. # Effective area of each circle (in Mpc^2)
    field_density = field_density / (am_to_rad * troughDa)**2.


    rhotheta_tot = Ngaltheta_tot/Aefftheta_tot

    # Percentile and delta = (rho - av(rho))/av(rho)
    sort_rho = np.array([np.argsort(rhotheta_tot[theta]) for theta in range(Ntheta)])
    Ptheta_tot = np.array([np.argsort(gridID_tot[sort_rho[theta]])/float(Ngrid_tot) for theta in range(Ntheta)])
    delta_tot = (rhotheta_tot - field_density[-1]) / field_density[-1]

    # Following the definition from DES, we usually define the 20% of all circles with the lowest density to be the troughs.
    # The 20% of all circles with the highest density are defined as overdensities ("ridges").

    ## Write catalog

    ## 5b) Now we have two trough samples:
    # - The overlapping sample, by taking all troughs
    # - The non-overlapping sample, by taking only unflagged troughs (not flagged = 1, flagged = 0).

    # Writing the combined columns to a fits file
    filename = '/data2/brouwer/MergedCatalogues/trough_catalogs/trough_catalog_%s_%s_%s.fits'%(cat, selection, masktype)

    # For each grid point/circle, we save the following information to the catalog:
    # the location (RA/DEC), number count (Ngaltheta), effective area in arcmin (grid count Ngtheta), galaxy density (rhotheta), 
    # percentage of circles below its density (Ptheta), and non-overlapping selection flag (Stheta).

    # Defining the names of the columns
    outputnames = ['RA', 'DEC']
    
    [outputnames.append('Ngaltheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('Aefftheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('Pmasktheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('rhotheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('Ptheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('delta%g'%(theta*60)) for theta in thetalist]
    #[outputnames.append('Stheta%g'%(theta*60)) for theta in thetalist]

    # Defining the output
    output = [gridRA_tot, gridDEC_tot]
    
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

    filename = 'density_info/density_info_%s_%s_%s.txt'%(cat, selection, masktype)

    np.savetxt(filename, density_info, delimiter='    ', header = field_header, footer = field_footer)
    print('Written:', filename)
