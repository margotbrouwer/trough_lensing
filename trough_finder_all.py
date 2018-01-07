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

# Radii theta of circular regions (in deg)
#thetalist = np.array([5., 10., 15., 20.])/60.
#thetalist = np.array([5., 10., 15., 20., 6.303, 6.288])/60.
thetalist = np.array([5., 10., 6.303, 6.288])/60.
#thetalist = np.array([20., 12.85, 9.45, 7.44, 6.14])
#thetalist = np.array([5.])/60.


ijlist = np.array([ [ [i, j] for i in range(4) ] for j in range(4) ])
ijlist = np.reshape(ijlist, [16,2])

Nruns = len(ijlist) # Mock patches
#Nruns = len(thetalist) # Theta
#Nruns = 1
#Nruns = 5

for ij in np.arange(0, Nruns):

    # Calculate i and j
    i, j = ijlist[ij]
    ijnum = ij+1
    print(i, j, ijnum)

    ### 1) Defining the galaxy sample:

    ## 1a) Importing the galaxy catalogue.

    # Select the galaxy catalogue for trough selection (kids, gama or mice)
    #cat = 'gama'
    #cat = 'kids'
    cat = 'mice'

    # Name of the pre-defined galaxy selection
    #selection = 'all'
    #selection = 'absmag'
    #selection = 'mice'
    #selection = 'lowZ'
    selection = 'highZ'
    #selection = 'miceZ-%g'%(ijnum)
    
    #selection = 'gama_all'
    #selection = 'kids_all'
    
    
    # Select mask type (nomask or complex)
    if 'mice' in cat:

        if 'miceZ' in selection:
            masktype = 'nomask-Z'
        else:
            masktype = 'nomask-1'

        if Nruns > 14:
            masktype = 'nomask-%g'%(ijnum)
            coordsM = [[i*20.,(i+1)*20.], [j*20.,(j+1)*20.]]
        else:
            coordsM = [[0.,20.], [0.,20.]]
            
            #Wlist = np.array([83.539, 40., 25.015, 17.948, 14.086])
            #highW = Wlist[ij]
            #coordsM = [[0.,highW], [0.,20.]]
    else:
        masktype = 'complex'

    # Spacing of the trough and mask grids (in degree)
    gridspace = 0.04
    mask_density = 1/(gridspace*60.)**2 # Density of mask gridpoints (in arcmin^-2)


    # Redshift samples
    zmin = 0.1
    zlim = 0.198
    zmax = 0.3


    # Import galaxy catalog
    if 'kids' in cat:
        
        if 'gama' in selection:
            # Names of the GAMA fields
            fieldnames = ['G9', 'G12', 'G15']
            
            # Boundaries of the KiDS fields
            coordsG9 = [[128.0,142.5], [-2.5,3.5]]
            coordsG12 = [[155.0,190.0], [-3.5,3.5]]
            coordsG15 = [[209.5,239.0], [-3.5,3.5]]
            fieldboundaries = np.array([coordsG9,coordsG12,coordsG15]) # Boundaries of all fields
        else:
            # Names of the KiDS fields
            fieldnames = ['G9', 'G12', 'G15', 'G23', 'GS']
            
            # Boundaries of the KiDS fields
            coordsG9 = [[128.0,142.5], [-2.5,3.5]]
            coordsG12 = [[155.0,190.0], [-3.5,3.5]]
            coordsG15 = [[209.5,239.0], [-3.5,3.5]]
            coordsG23 = [[328.0,361.0], [-35.0,-28.0]]
            coordsGS = [[31.0,54.0], [-35.0,-29.0]]
            fieldboundaries = np.array([coordsG9,coordsG12,coordsG15,coordsG23,coordsGS]) # Boundaries of all fields
        
        # Path to the KiDS fields
        path_kidscat = '/data2/brouwer/KidsCatalogues'
        if 'nomask' in masktype:
            kidscatname = '/KiDS_DR3_GAMA-like_Maciek_NOMASKING_01.06.17-withNEWzANNz2.fits'
        else:
            kidscatname = '/KiDS_DR3_GAMA-like_Maciek_revised_1905.fits'
        
        # Importing the KiDS galaxies
        galRA, galDEC, galZB, galZ, galTB, mag_auto, ODDS, umag, gmag, rmag, imag = \
        utils.import_kidscat(path_kidscat, kidscatname)
        
        # Calculating the absolute magnitudes
        h, O_matter, O_lambda = [0.7, 0.25, 0.75]
        rmag_abs = utils.calc_absmag(rmag, galZ, gmag, imag, h, O_matter, O_lambda)
        
        gama_rlim = 20.2
        
        """
        # Write test absmag_r catalogue
        outputnames = ['RA', 'DEC', 'Z', 'mag_auto', 'mag_r', 'absmag_r']
        output = [galRA, galDEC, galZ, mag_auto, rmag, rmag_abs]

        testfilename = '/data2/brouwer/KidsCatalogues/absmag_test.fits'
        print('Writing output catalogue...')
        utils.write_catalog(testfilename, np.arange(len(galRA))+1, outputnames, output)
        """

    if 'gama' in cat:
        
        # Names of the GAMA fields
        fieldnames = ['G9', 'G12', 'G15']
        
        # Boundaries of the GAMA fields
        coordsG9 = [[129., 141.], [-2.,3.]]
        coordsG12 = [[174., 186.], [-3.,2.]]
        coordsG15 = [[211.5, 223.5], [-2.,3.]]
        fieldboundaries = np.array([coordsG9,coordsG12,coordsG15])

        # Path to the GAMA fields
        path_gamacat = '/data2/brouwer/MergedCatalogues'
        gamacatname = 'ShearMergedCatalogueAll_sv0.8.fits'
        
        # Importing the GAMA galaxies
        galRA, galDEC, galZ, rmag, rmag_abs = utils.import_gamacat(path_gamacat, gamacatname)
        
        gama_rlim = 19.8

    # Import galaxy catalog
    if 'mice' in cat:
        
        # Names of the GAMA fields
        fieldnames = ['M']
        
        # Boundaries of the MICE field
        fieldboundaries = np.array([coordsM]) # Boundaries of all fields
        print(coordsM)
        
        # Path to the Mice field
        path_mockcat = '/data2/brouwer/MergedCatalogues'
        if 'miceZ' in selection:
            mockcatname = 'mice_gama_highZ_catalog.fits'
        else:
            mockcatname = 'mice_gama_catalog.fits'
        
        # Importing the Mice galaxies
        galRA, galDEC, galZ, galDc, rmag, rmag_abs, e1, e2, galmass = \
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

    print(maskfilename)



    ## 1b) Select the galaxy sample to define the troughs

    # Defining the selection for the galaxy sample

    if selection == 'all':
        galmask = (galZ < 0.5) & (rmag < gama_rlim)
    if selection == 'absmag':
        galmask = (galZ < 0.5) & (rmag_abs < -19.7) & (rmag < gama_rlim)
    if 'mice' in selection:
        galmask = (galZ < 0.5) & (rmag_abs < -18.9 + micecor) & (rmag < gama_rlim)
        print('rmag_abs(MICE) =', -18.9 + micecor)

    if 'lowZ' in selection:
        galmask = (zmin < galZ)&(galZ < zlim) & (rmag_abs < -21.) & (rmag < gama_rlim)
        #thetalist = np.array([5., 10., 15., 20.])/60.
        thetalist = np.array([10.])/60.
    
    if 'highZ' in selection:
        galmask = (zlim < galZ)&(galZ < zmax) & (rmag_abs < -21.) & (rmag < gama_rlim)
        #thetalist = np.array([3.41290231, 6.82580462, 10.23870693, 13.65160924])/60.
        thetalist = np.array([6.303, 6.288])/60.

    
    if 'miceZ' in selection:
        zlims = np.array([0.1, 0.191, 0.286, 0.385, 0.489, 0.6])
        zmask = (zlims[ij] < galZ)&(galZ < zlims[ij+1])
        
        # AbsMag cut
        galmask = zmask * (rmag_abs < -21.)
        
        # DM halo mass
        #galmask = zmask * (galmass > 12.5)
        
        # AbsMag cut with Luminosity evolution
        #Z = np.mean(galZ[zmask])
        #dM = -2.5*np.log10(1+Z)
        #galmask = zmask * (rmag_abs < (-21.+dM))
        #print('Abs. Magnitude limit:', -21.+dM)
        
        miceZthetalist = np.array([20., 12.85, 9.45, 7.44, 6.14])/60.
        thetalist = np.array([miceZthetalist[ij]])
        
        print('Zlims:', zlims[ij], zlims[ij+1], ', theta:', thetalist[0]*60.)

    
    if 'kids' in selection:
        # Path to the KiDS fields
        path_kidscat = '/data2/brouwer/KidsCatalogues'
        kidscatname = '/KiDS_DR3_GAMA-like_Maciek_NOMASKING_01.06.17-withNEWzANNz2.fits'
        
        # Importing the KiDS coordinates
        galRA, galDEC, galZ, galTB, mag_auto, ODDS, umag, gmag, rmag, imag = \
        utils.import_kidscat(path_kidscat, kidscatname)
        
        if 'all' in selection:
            galmask = (galZ >= 0.)
        if 'absmag' in selection:
            rmag_abs = calc_absmag()
            galmask = (rmag_abs < -19.7)

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
        
        if 'mice' in cat:
            # Creating the grid for each field
            gridspace_mask = 0.04 # in degree
            gridRA, gridDEC, gridcoords = utils.define_gridpoints(fieldRAs, fieldDECs, gridspace_mask, True)

            Ngrid = len(gridcoords)    
            catmask = np.ones(Ngrid)
        else:
            # Import the mask coordinates of this field
            path_catmask = '/data2/brouwer/MergedCatalogues/Masks/%s_mask_%s_%gdeg.fits'%(cat, fieldnames[field], gridspace)
            
            maskcat = pyfits.open(path_catmask, memmap=True)[1].data
            gridRA, gridDEC, catmask = [maskcat['RA'], maskcat['DEC'], maskcat['mask']]
            gridcoords = SkyCoord(ra=gridRA*u.deg, dec=gridDEC*u.deg)
            Ngrid = len(gridcoords)
        
        
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
    
    troughZ = galZ[galmask]
    troughZ_tot = np.array([np.mean(troughZ)] * Ngrid_tot)
    
    if 'Z' in selection:
        if 'mice' in cat:
            troughDc = galDc[galmask]
        else:
            troughDc = (cosmo.comoving_distance(troughZ).to('Mpc')).value

        troughDa = np.mean( troughDc/(1+troughZ) )
        
        am_to_rad = np.pi/(60.*180.)
        Aefftheta_tot = (Nmasktheta_tot/mask_density) * (am_to_rad * troughDa)**2. # Effective area of each circle (in Mpc^2)
        field_density = field_density / (am_to_rad * troughDa)**2.
    
    else:
        Aefftheta_tot = Nmasktheta_tot/mask_density # Effective area of each circle (in arcmin^2)

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
    outputnames = ['RA', 'DEC', 'Z']
    [outputnames.append('Ngaltheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('Aefftheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('Pmasktheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('rhotheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('Ptheta%g'%(theta*60)) for theta in thetalist]
    [outputnames.append('delta%g'%(theta*60)) for theta in thetalist]
    #[outputnames.append('Stheta%g'%(theta*60)) for theta in thetalist]

    # Defining the output
    output = [gridRA_tot, gridDEC_tot, troughZ_tot]
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
