#!/usr/bin/python

import astropy.io.fits as pyfits
import gc
import numpy as np
import sys
import os
import time
from glob import glob

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

from matplotlib import gridspec
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

import scipy.optimize as optimization
import trough_modules_all as utils
from astropy.cosmology import LambdaCDM

h, O_matter, O_lambda = [0.7, 0.25, 0.75]
cosmo = LambdaCDM(H0=h*100, Om0=O_matter, Ode0=O_lambda)
am_to_rad = np.pi/(60.*180.)

# Make use of TeXm
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})

colors = ['red', 'orange', 'cyan', 'blue']
#colors = ['#d7191c', '#fdae61', '#92c5de', '#0571b0']

# Model to fit the troughs
def trough_model(x, A):
    model_y = A*x**-0.5
    return model_y

show = True

ijlist = np.array([ [ [i, j] for i in range(4) ] for j in range(4) ])
ijlist = np.reshape(ijlist, [16,2])
thetalist = np.array([5., 10., 15., 20.]) # in arcmin

#Nruns = len(ijlist) # Mock patches
#Nruns = len(thetalist) # Theta
#Nruns = 5
Nruns = 1

# Configuration

for ij in np.arange(0, Nruns):
    
    # Number of the current run
    ijnum = ij+1
    
    # Defining the paths to the data
    blind = 'A'

    #selection = 'kids_mice_complex'
    #selection = 'kids_absmag_complex'
    selection = 'kids_highZ_complex'
    #selection = 'mice_all_nomask-%g'%ijnum
    #selection = 'mice_highZ_nomask-%g'%ijnum
    #selection = 'mice_miceZ-%g_nomask-Z'%ijnum
    #selection = 'mice_miceZa_nomask-Za-%g'%ijnum
    #selection = 'slics_mocks_nomask'
    #selection = 'slics_mockZ_nomask'
    
    mocksel = 'mice_highZ_nomask-%g'%ijnum
    randomsel = 'kids_highZ_complex'

    # Select unit (arcmin or Mpc)
    if 'Z' in selection:
        Runit = 'Mpc'
    else:
        Runit = 'arcmin'
    
    
    if ('all' in selection) or ('absmag' in selection) or ('mice' in selection) or ('slics' in selection):
        thetalist = np.array([5., 10., 15., 20.]) # in arcmin
        thetanum = ij
        
    if 'lowZ' in selection:
        thetalist = np.array([10.])
        thetanum = 0
        
    if 'highZ' in selection:
        thetalist = np.array([6.288])
        thetanum = 0

    if 'miceZ' in selection:
        thetalist = np.array([20., 12.85, 9.45, 7.44, 6.14]) # Dc
        thetanum = ij

    if 'mockZ' in selection:
        thetalist = np.array([15., 9.554, 7.283, 5.770]) # Dc
        thetanum = ij
        h, O_matter, O_lambda = [0.7, 0.29, 0.71]
        cosmo = LambdaCDM(H0=h*100, Om0=O_matter, Ode0=O_lambda)
        
        
    theta = thetalist[thetanum]
    
    if Runit == 'arcmin':
        Rmin = 2
        Rmax = 100
        Nbins = 20
        
        dperc = 0.05
        percnames = ['0','0p05','0p1','0p15','0p2','0p25','0p3','0p35','0p4',\
        '0p45','0p5','0p55','0p6','0p65','0p7','0p75','0p8','0p85','0p9','0p95','1']

    if Runit == 'Mpc':
        Rmin = 0.5
        Rmax = 20
        Nbins = 10
        
        dperc = 0.1
        percnames = ['0','0p1','0p2','0p3','0p4','0p5','0p6','0p7','0p8','0p9','1']

    # Defining the percentile bins
    percmin = 0.
    percmax = 1.
    
    perclist = np.arange(percmin, percmax, dperc)
    perccenters = perclist+dperc/2.
    perclist = np.append(perclist, percmax)

    Npercs = len(percnames)-1

    # Import trough catalog
    if 'slics' in selection:
        mean_zs = np.array([0.1539674, 0.24719192, 0.33112174, 0.42836386])
        troughZ = mean_zs[ij]
    else:
        path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'
        troughcatname = 'trough_catalog_%s.fits'%(selection)
        troughRA, troughDEC, troughZ, paramlists = utils.import_troughcat(path_troughcat, troughcatname, [])
        troughZ = troughZ[0]

    # Import lensing profiles

    if ('kids' in selection) or ('gama' in selection):
      
        # Observed lensing profiles
        path_sheardata = 'data2/brouwer/shearprofile/trough_results_final'
        path_lenssel = ['No_bins_%s/Pmasktheta%s_0p8_inf-Ptheta%s_%s_%s'\
            %(selection, ('%g'%theta).replace('.','p'), ('%g'%theta).replace('.','p'), percnames[p], percnames[p+1]) for p in range(Npercs)]
        path_cosmo = 'ZB_0p1_0p9-Om_0p25-Ol_0p75-Ok_0-h_0p7/Rbins%i_%s_%s_%s/shearcovariance'%(Nbins, ('%g'%Rmin).replace('.','p'), ('%g'%Rmax).replace('.','p'), Runit)
        path_filename = 'No_bins_%s.txt'%(blind)

        esdfiles = np.array([('/%s/%s/%s/%s'%(path_sheardata, path_lenssel[p], path_cosmo, path_filename)) \
                   for p in range(len(path_lenssel))])

        covfiles = np.array([e.replace('bins_%s.txt'%blind, 'matrix_%s.txt'%blind) for e in esdfiles])
        covariance_tot = np.array([ np.loadtxt(covfiles[c]).T for c in range(len(covfiles)) ])

        # Mock lensing profiles
        path_mockdata = 'data2/brouwer/shearprofile/trough_results_final'
        path_mocksel = ['No_bins_%s/Pmasktheta%s_0p8_inf-Ptheta%s_%s_%s'\
            %(mocksel, ('%g'%theta).replace('.','p'), ('%g'%theta).replace('.','p'), percnames[p], percnames[p+1]) for p in range(Npercs)]
        mockfiles = np.array([('/%s/%s.txt'%(path_mockdata, path_mocksel[p])) \
                   for p in range(len(path_mocksel))])
        
    else:
        if 'mice' in selection:
            # Mock lensing profiles
            path_sheardata = 'data2/brouwer/shearprofile/trough_results_final'
            path_lenssel = ['No_bins_%s/Pmasktheta%s_0p8_inf-Ptheta%s_%s_%s'\
                %(selection, ('%g'%theta).replace('.','p'), ('%g'%theta).replace('.','p'), percnames[p], percnames[p+1]) for p in range(Npercs)]
            esdfiles = np.array([('/%s/%s.txt'%(path_sheardata, path_lenssel[p])) \
                       for p in range(len(path_lenssel))])
                       
        if 'slics' in selection:

            path_sheardata = 'data2/brouwer/shearprofile/trough_results_final'
            path_lenssel = '%s'%(selection)

            if 'mockZ' in selection:
                esdfiles = 'Redshift_bins_covariance.npy'
                data=np.load('/%s/%s/%s'%(path_sheardata, path_lenssel, esdfiles))
                
                error_factor = 100./15000.
                covariance_tot = 0.7 * error_factor * np.array((data[2][ijnum]).values())
                
                theta_x = np.array(data[0][ijnum][0])
                data_x = [theta_x*am_to_rad * (cosmo.angular_diameter_distance(troughZ).to('Mpc')).value] # Physical radial distance
                
                data_y = 0.7 * np.array(data[1][ijnum].values()) # Convert h100 to h70

            else:
                esdfiles = 'Ptheta%s.npy'%(('%g'%theta).replace('.','p'))
                data=np.load('/%s/%s/%s'%(path_sheardata, path_lenssel, esdfiles))

                error_factor = 100./360.3                
                covfiles = 'Ptheta%s_cov.npy'%(('%g'%theta).replace('.','p'))
                covariance_tot = error_factor *\
                    np.array((np.load('/%s/%s/%s'%(path_sheardata, path_lenssel, covfiles))[0]).values())
                
                data_x = data[0]
                data_y = np.array(data[1].values())

                
    path_plots = '/%s/Plots/%s'%(path_sheardata, selection)
    
    # Importing the shearprofiles and lens IDs
    if 'slics' in selection:
        
        errors = [np.sqrt(np.diag(covariance_tot[x])) for x in range(len(covariance_tot))]
        error_h = np.array(errors)
        error_l = np.array(errors)
        
    else:
        print('Import shear signal:', esdfiles[0])
        data_x, data_y, error_h, error_l = utils.read_esdfiles(esdfiles)

    data_x = data_x[0]
    
    if 'pc' in Runit:
        # Translate to comoving ESD
        data_x = data_x*(1+troughZ)
        data_y, error_h, error_l = np.array([data_y, error_h, error_l])/(1+troughZ)**2
        covariance_tot = covariance_tot/(1+troughZ)**4
    
    try:
        print('Import mock signal:')
        mock_x, mock_y, mock_error_h, mock_error_l = utils.read_esdfiles(mockfiles)
        mock_x = mock_x[0]

        if 'pc' in Runit:
            # Translate to comoving ESD
            mock_x = mock_x*(1+troughZ)
            mock_y, mock_error_h, mock_error_l = np.array([mock_y, mock_error_h, mock_error_l])/(1+troughZ)**2

    except:
        pass

    if ('kids' in selection) or ('gama' in selection):
        # Import random signal

        print('Import random signal:')
        path_randoms = ['No_bins_%s/Pmasktheta%s_0p8_inf'%(randomsel, ('%g'%theta).replace('.','p'))]
        random_esdfile = np.array(['/%s/%s/%s/%s'%(path_sheardata, path_random, path_cosmo, path_filename) for path_random in path_randoms])
        random_data_x, random_data_y, random_error_h, random_error_l = utils.read_esdfiles(random_esdfile)
        random_data_x, random_data_y, random_error_h, random_error_l = random_data_x[0], random_data_y[0], random_error_h[0], random_error_l[0]

        # Subtract random signal
        data_y = data_y-random_data_y
        error_h = np.sqrt(error_h**2. + random_error_h**2)
        error_l = np.sqrt(error_l**2. + random_error_l**2)


    # Define the part of the trough profile that contributes to the fit
    if Runit == 'arcmin':
        xmin = theta*1.2
        xmax = 70.
    if Runit == 'Mpc':
        Rlist = theta*am_to_rad * (cosmo.comoving_distance(troughZ).to('Mpc')).value #2.77797224336
        print(Rlist)
        xmin = Rlist*1.2
        xmax = 14.
    
    xmask = (xmin < data_x) & (data_x < xmax)
    xwhere = np.where(xmask)[0]

    # Plotting the ueber matrix
    Nbins = Npercs
    Nrows = Nbins/5
    Ncolumns = int(Nbins/Nrows)
    
    fig = plt.figure(figsize=(12,8))
    canvas = FigureCanvas(fig)

    gs_full = gridspec.GridSpec(1,1)
    gs = gridspec.GridSpecFromSubplotSpec(Nrows, Ncolumns, wspace=0, hspace=0, subplot_spec=gs_full[0,0])

    ax = fig.add_subplot(gs_full[0,0])

    # Fit the model of every trough, and plot the result
    Alist = []
    Alist_error = []
    for N1 in range(Nrows):
        for N2 in range(Ncolumns):
        
            N = np.int(N1*Ncolumns + N2)
            print(N)
            ax_sub = fig.add_subplot(gs[N1, N2])
          
            if ('kids' in selection) or ('gama' in selection):
                # With covariance
                covariance = covariance_tot[N]
                ind = np.lexsort((covariance[3,:], covariance[1,:], covariance[2,:], covariance[0,:]))
                covmatrix = np.reshape(covariance[4][ind], [len(data_x), len(data_x)])
                covmatrix = covmatrix[int(xwhere[0]):int(xwhere[-1]+1), int(xwhere[0]):int(xwhere[-1]+1)]
                
                A, Acov = optimization.curve_fit(f=trough_model, xdata=data_x[xmask], ydata=(data_y[N])[xmask], p0=[0.], \
                sigma=covmatrix, absolute_sigma=True)
            else:
                if 'slics' in selection:
                    covariance = np.array(covariance_tot[N])
                    covmatrix = covariance[int(xwhere[0]):int(xwhere[-1]+1), int(xwhere[0]):int(xwhere[-1]+1)]
                    
                    A, Acov = optimization.curve_fit(f=trough_model, xdata=data_x[xmask], ydata=(data_y[N])[xmask], p0=[0.], \
                    sigma=covmatrix, absolute_sigma=True)
                    
                    #A, Acov = optimization.curve_fit(f=trough_model, xdata=data_x[xmask], ydata=(data_y[N])[xmask], p0=[0.], \
                    #sigma=(error_h[N])[xmask], absolute_sigma=True)
                    
                    #A, Acov = optimization.curve_fit(f=trough_model, xdata=data_x[xmask], ydata=(data_y[N])[xmask], p0=[0.])
                    
                else:
                    # Without covariance
                    A, Acov = optimization.curve_fit(f=trough_model, xdata=data_x[xmask], ydata=(data_y[N])[xmask], p0=[0.])
                    #sigma=(error_h[N])[xmask], absolute_sigma=True)
            A = A[0]
                
            print('%g < P(x) < %g: Amplitude = %g'%(perclist[N], perclist[N+1], A))
            
            Alist = np.append(Alist, A)
            Alist_error = np.append(Alist_error, np.sqrt(Acov[0,0]))

            model_x = np.linspace(xmin, xmax, 10)
            model_y = trough_model(model_x, A)

            # Plot mock profile
            try:
                ax_sub.plot(mock_x, mock_y[N], ls='-', color='black', zorder=1)
            except:
                pass
            
            # Plot fitted model
            ax_sub.plot(model_x, model_y, ls='-', color='red', label=r'$%g < P \leq %g$'%(perclist[N], perclist[N+1]), zorder=6)
            
            # Plot observed profile
            ax_sub.errorbar(data_x, data_y[N], yerr=error_h[N], ls='', marker='.', color=colors[3], zorder=4)
            
            # Plot lines
            ax_sub.axvline(x=xmin, ls=':', color='black', zorder=1)
            ax_sub.axvline(x=xmax, ls=':', color='black', zorder=2)
            ax_sub.axhline(y=0., ls=':', color='black', zorder=3)

            ax_sub.xaxis.set_label_position('bottom')
            ax_sub.yaxis.set_label_position('left')

            ax.tick_params(labelleft='off', labelbottom='off', top='off', bottom='off', left='off', right='off')
            
            plt.legend(loc='best', handlelength=0, handletextpad=0, numpoints=1)
            
            # Last row
            if N1 != Nrows-1:
                ax_sub.tick_params(axis='x', labelbottom='off')
                
            # First column
            if N2 != 0:
                ax_sub.tick_params(axis='y', labelleft='off')

            if Runit == 'arcmin':
                plt.ylim(-0.002, 0.0039)
            if Runit == 'Mpc':
                plt.ylim(-4., 7.9)

            plt.xlim(Rmin, Rmax)
            
            plt.xscale('log')

    if Runit == 'arcmin':
        ax.set_xlabel(r'Radial separation $\theta$ (arcmin)', fontsize=14)
        ax.set_ylabel(r'Shear $\gamma$', fontsize=14)
    if Runit == 'Mpc':
        ax.set_xlabel(r'Radial distance $R ({\rm Mpc} \, {h_{70}}^{-1})$', fontsize=14)
        ax.set_ylabel(r'Shear $\gamma$', fontsize=14)
        
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.07, 0.5)

    for ext in ['png', 'pdf']:

        plotfilename = '%s_trough_amplitude_fits_%g%s'%(path_plots, theta, Runit)
        plotname = '%s.%s'%(plotfilename, ext)

        plt.savefig(plotname, format=ext, bbox_inches='tight')
        
    print('Written plot:', plotname)
    if show:
        plt.show()


    # Import trough catalog
    path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'
    troughcatname = 'trough_catalog_%s_%s_%s.fits'%(selection.split('_')[0], selection.split('_')[1], selection.split('_')[2])
    print(troughcatname)

    # Full directory & name of the trough catalogue
    troughcatfile = '%s/%s'%(path_troughcat, troughcatname)
    troughcat = pyfits.open(troughcatfile, memmap=True)[1].data

    Ptheta = troughcat['Ptheta%g'%theta]
    delta = troughcat['delta%g'%theta]
    Pmasktheta = troughcat['Pmasktheta%g'%theta]


    # Calculate the mean delta for each percentile bin
    deltacenters = np.zeros(len(perccenters))
    for p in range(len(perccenters)):
        percmask = (perclist[p] < Ptheta) & (Ptheta <= perclist[p+1]) & (0.8 < Pmasktheta)
        deltacenters[p] = np.mean(delta[percmask])


    # Write amplitude text-file
    Afilename = '/%s/Plots/trough_amplitudes_%s_%g%s.txt'%(path_sheardata, selection, theta, Runit)
    np.savetxt(Afilename, np.array([perccenters, deltacenters, Alist, Alist_error]).T, header = 'Trough Percentile     Delta     Amplitude     Error(Amplitude)')
    print('Written textfile:', Afilename)
