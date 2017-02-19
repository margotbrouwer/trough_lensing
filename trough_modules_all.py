#!/usr/bin/python
"""
# This contains all the modules for the trough project.
"""
import astropy.io.fits as pyfits
import gc
import numpy as np
import sys
import os
import time
from glob import glob

from astropy import constants as const, units as u
from astropy.coordinates import SkyCoord

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

# Import and mask all used data from the sources in this KiDS field
def import_kidscat(path_kidscat, kidscatname):
    
    # Full directory & name of the corresponding KiDS catalogue
    kidscatfile = '%s/%s'%(path_kidscat, kidscatname)
    kidscat = pyfits.open(kidscatfile, memmap=True)[1].data
    
    # List of the observables of all sources in the KiDS catalogue
    srcRA = kidscat['RAJ2000']
    srcDEC = kidscat['DECJ2000']
    srcZB = kidscat['Z_B_BPZ']
    srcTB = kidscat['T_B_BPZ']
    ODDS = kidscat['ODDS_BPZ'] # Quality of the photometric redshift
    mag_auto = kidscat['MAG_AUTO_R']
    
    umag = kidscat['MAG_GAAP_u']
    gmag = kidscat['MAG_GAAP_g']
    rmag = kidscat['MAG_GAAP_r']
    imag = kidscat['MAG_GAAP_i']

    # Adding homogenization ZPT offset and subtracting Galactic foreground extinction following Schlegel et al. maps.
    umag = umag + kidscat['ZPT_OFFSET_U'] - kidscat['EXT_SFD_U']
    gmag = gmag + kidscat['ZPT_OFFSET_G'] - kidscat['EXT_SFD_G']
    rmag = rmag + kidscat['ZPT_OFFSET_R'] - kidscat['EXT_SFD_R']
    imag = imag + kidscat['ZPT_OFFSET_I'] - kidscat['EXT_SFD_I']

    umagerr = kidscat['MAGERR_GAAP_u']
    gmagerr = kidscat['MAGERR_GAAP_g']
    rmagerr = kidscat['MAGERR_GAAP_r']
    imagerr = kidscat['MAGERR_GAAP_i']
    
    uflag = kidscat['IMAFLAGS_ISO_U']
    gflag = kidscat['IMAFLAGS_ISO_G']
    rflag = kidscat['IMAFLAGS_ISO_R']
    iflag = kidscat['IMAFLAGS_ISO_I']
    
    class_star = kidscat['CLASS_STAR']
    
    # Masks
    srcmask = (np.abs(umag) < 90.) & (np.abs(gmag) < 90.) & (np.abs(rmag) < 90.) & (np.abs(imag) < 90.) & \
    (uflag==0)&(gflag==0)&(rflag==0)&(iflag==0) & (class_star < 0.8)# & \
    #(0.<umagerr)&(umagerr<1.) & (0.<gmagerr)&(gmagerr<1.) & (0.<rmagerr)&(rmagerr<1.) & (0.<imagerr)&(imagerr<1.)
    
    print('Imported: %i of %i sources \n from %s'%(np.sum(srcmask), len(srcRA), kidscatfile))

    srcRA = srcRA[srcmask]
    srcDEC = srcDEC[srcmask]
    srcZB = srcZB[srcmask]
    srcTB = srcTB[srcmask]
    
    umag = umag[srcmask]
    gmag = gmag[srcmask]
    rmag = rmag[srcmask]
    imag = imag[srcmask]
    mag_auto = mag_auto[srcmask]
    ODDS = ODDS[srcmask]

    return srcRA, srcDEC, srcZB, srcTB, mag_auto, ODDS, umag, gmag, rmag, imag


def import_gamacat(path_gamacat, gamacatname):
    
    # Full directory & name of the corresponding KiDS catalogue
    gamacatfile = '%s/%s'%(path_gamacat, gamacatname)
    gamacat = pyfits.open(gamacatfile, memmap=True)[1].data
    
    # List of the observables of all sources in the KiDS catalogue
    galRA = gamacat['RA']
    galDEC = gamacat['DEC']
    galZ = gamacat['Z_1']
    
    return galRA, galDEC, galZ
    

# Define grid points for trough selection
def define_gridpoints(fieldRAs, fieldDECs, srccoords, gridspace, gridmax):

    ## Grid
    
    # Creating a grid to measure the galaxy density
    gridRAlist = np.arange(fieldRAs[0], fieldRAs[1], gridspace)
    gridDEClist = np.arange(fieldDECs[0], fieldDECs[1], gridspace)
    gridlist = np.vstack([np.array([np.array([RA, DEC]) for RA in gridRAlist]) for DEC in gridDEClist])
    
    gridRA, gridDEC = gridlist[:,0], gridlist[:,1]
    gridcoords = SkyCoord(ra=gridRA*u.deg, dec=gridDEC*u.deg) # All grid coordinates    

    print('Number of grid coordinates:', len(gridRAlist), 'x', len(gridDEClist), '=', len(gridcoords))

    ## Masking grid points
    
    # Distances of grind points to the nearest source
    idx, d2dsrc, d3d = gridcoords.match_to_catalog_sky(srccoords)
    
    # Find grid points that are outside the field
    outmask = (d2dsrc > gridmax*u.deg) # Points that lie outside the source field
    outcoords = gridcoords[outmask]
    
    # Remove points that lie close to the edge of the galaxy field
    idx, d2dout, d3d = gridcoords.match_to_catalog_sky(outcoords)
    inmask = (d2dout > 2.*gridmax*u.deg)
    
    
    # Remove grid points that are in masked areas or outside the field
    starmask = (d2dsrc < gridmax*u.deg)
    Nmaskedgrid = np.sum(starmask)    
    gridmask = starmask*inmask

    # Define new grid coordinates
    gridRA, gridDEC, gridcoords = gridRA[gridmask], gridDEC[gridmask], gridcoords[gridmask]
    print('New gridcoords:', len(gridcoords))
    
    return gridRA, gridDEC, gridcoords, Nmaskedgrid


# Defining the mask for the galaxy sample
def define_galsamp(selection, zmin, zmax, srcZB, srcTB, gmag, rmag, mag_auto):
    
    gminr = gmag - rmag
    
    colormask = (0.5 < gminr)&(gminr < 2.)
    magmask = (0. < mag_auto)&(mag_auto<23.)
    zmask = (zmin < srcZB)&(srcZB < zmax)
    Tmask = (srcTB < 1.5)

    elmask = Tmask*colormask*magmask

    if selection == 'ell':
        print('Selection: Ellipticals')
        redmask = elmask*zmask
    
    if 'redseq' in selection:
        print('Selection: Red Sequence')
        
        # Find the red sequence
        if '4' in selection:        
            redseq = np.loadtxt('redseq_Ascaso.txt').T
            problim = 0.75
        else:
            redseq = np.loadtxt('redseq_tot.txt').T
            problim = 0.5
        
        # Determine the probability of galaxies to be on the red sequence
        zlist, dz, redseq_dist, redseq_prob, dist_dev0 = calc_redseq_prob(redseq, srcZB, gminr, mag_auto, elmask)
        
        probmask = (redseq_prob > problim)
        redmask = probmask*colormask*magmask
    
    if selection == 'all':
        print('Selection: All galaxies')
        redmask = zmask*magmask
    
    return redmask
    


# Define the red sequence lines and probability of being an LRG
def calc_redseq_prob(redseq, srcZB, gminr, mag_auto, elmask):

    Alist, Blist, Clist = redseq[1::]
    
    zlist = redseq[0]
    dz = zlist[1]-zlist[0]

    redseq_dist = np.zeros(len(srcZB))
    redseq_prob = np.zeros(len(srcZB))
    dist_dev0 = []
        
    print('The A, B and C of the red sequence lines, where: A*x + B*y + C = 0.')
    print('Z_B          A                 B                  C               Sigma(T_B<1.5)')
    
    for r in range(len(zlist)):
        
        zmask = ((zlist[r]-dz/2.) < srcZB)&(srcZB < (zlist[r]+dz/2.))

        # Calculating the distance between the galaxies at redshift Z and their red sequence.
        point_dist = np.abs(Alist[r]*mag_auto[zmask] + Blist[r]*gminr[zmask] + Clist[r])/np.sqrt(Alist[r]**2 + Blist[r]**2)
        redseq_dist[zmask] = point_dist
        
        # Defining the elliptical sample
        elmask_z = zmask*elmask
        
        # For each red sequence, calculate the standard deviation of ellipticals
        dist_dev0 = np.append(dist_dev0, np.sqrt(np.mean(np.abs(redseq_dist[elmask_z] - 0.)**2)))
        
        # For each galaxy, calculate the probability of lying on the red sequence
        sigma = dist_dev0[r]
        redseq_prob[zmask] = np.exp(-(redseq_dist[zmask])**2./(2.*sigma**2.))

        print(zlist[r], '    ', Alist[r], '    ', Blist[r], '    ', Clist[r], '    ', dist_dev0[r])
        
    return zlist, dz, redseq_dist, redseq_prob, dist_dev0


# Write the results to a fits catalogue
def write_catalog(filename, galIDlist, outputnames, output):

    fitscols = []

    # Adding the lens IDs
    fitscols.append(pyfits.Column(name = 'ID', format='J', array = galIDlist))

    # Adding the output
    [fitscols.append(pyfits.Column(name = outputnames[c], format = '1D', array = output[c])) \
        for c in range(len(outputnames))]

    cols = pyfits.ColDefs(fitscols)
    tbhdu = pyfits.new_table(cols)

    #	print
    if os.path.isfile(filename):
        os.remove(filename)
        print('Old catalog overwritten:', filename)
    else:
        print('New catalog written:', filename)
    print()

    tbhdu.writeto(filename)


# Importing the ESD profiles
def read_esdfiles(esdfiles):
    
    data = np.loadtxt(esdfiles[0]).T
    data_x = data[0]

    data_x = []
    data_y = []
    error_h = []
    error_l = []
    
    print('Imported ESD profiles: %i'%len(esdfiles))
    
    for f in range(len(esdfiles)):
        # Load the text file containing the stacked profile
        data = np.loadtxt(esdfiles[f]).T
    
        bias = data[4]
        bias[bias==-999] = 1
    
        datax = data[0]
        datay = data[1]/bias
        datay[datay==-999] = np.nan
    
        errorh = (data[3])/bias # covariance error
        errorl = (data[3])/bias # covariance error
        errorh[errorh==-999] = np.nan
        errorl[errorl==-999] = np.nan
        
        data_x.append(datax)     
        data_y.append(datay) 
        error_h.append(errorh) 
        error_l.append(errorl) 
    
    return data_x, data_y, error_h, error_l
