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
import matplotlib.style
matplotlib.style.use('classic')

from matplotlib import gridspec
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid.inset_locator import inset_axes

from scipy.stats import chi2, norm
import scipy.optimize as optimization
import trough_modules_all as utils

# Make use of TeX
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})

#colors = ['red', 'orange', 'cyan', 'blue', 'green']
colors = ['#d7191c', '#fdae61', '#92c5de', '#0571b0', '#0ba52d']

Runit = 'arcmin'
show = True

theta = 20

percnames = ['0','0p05','0p1','0p15','0p2','0p25','0p3','0p35','0p4',\
'0p45','0p5','0p55','0p6','0p65','0p7','0p75','0p8','0p85','0p9','0p95','1']
percvals = np.arange(0.0, 1.05, 0.05)
#percnames = ['0','0p2']

#for Pbin in np.arange(0,20):
for Pbin in [0, 10, 19]:
#for Pbin in [0]: 
    
    print('Percentile bin %i: %s < P < %s'%(Pbin+1, percnames[Pbin], percnames[Pbin+1]))
    
    # Import KiDS errors
    path_filename_kids = '/data2/brouwer/shearprofile/trough_results_final/No_bins_kids_mice_complex/Pmasktheta%i_0p8_inf-Ptheta%i_%s_%s/ZB_0p1_0p9-Om_0p25-Ol_0p75-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'%(theta,theta,percnames[Pbin],percnames[Pbin+1])
    filename_kids = 'No_bins_A'

    data_kids = np.loadtxt('%s/%s.txt'%(path_filename_kids, filename_kids)).T
    data_x_kids = data_kids[0]
    data_y_kids = data_kids[3]
    print(data_y_kids)
    
    # Import GAMA errors
    path_filename_gama = '/data2/brouwer/shearprofile/trough_results_final/No_bins_gama_mice_complex/Pmasktheta5_0p8_inf-Ptheta5_0_0p2/ZB_0p1_0p9-Om_0p25-Ol_0p75-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
    filename_gama = 'No_bins_A'

    data_gama = np.loadtxt('%s/%s.txt'%(path_filename_gama, filename_gama)).T
    data_x_gama = data_gama[0]
    data_y_gama = data_gama[3]

    # Import mock errors
    path_filename_mock = '/data2/brouwer/shearprofile/trough_results_final/slics_mocks_nomask'

    """
    # Fiducial troughs (P<0.2)
    filename_mock = 'Ptheta5-0_0p2'
    data_mock = np.load('%s/%s.npy'%(path_filename_mock, filename_mock))
    data_x_mock, data_y_mock = [data_mock[0], data_mock[1]]
    """
    # All trough bins
    filename_mock = 'Ptheta%i'%(theta)
    data_mock = np.load('%s/%s.npy'%(path_filename_mock, filename_mock))
    data_x_mock = data_mock[0][Pbin]
    
    error_factor = 100./360.3
    covfiles = '%s_cov.npy'%filename_mock
    covariance_tot = error_factor * np.array((np.load('/%s/%s'%(path_filename_mock, covfiles))[0]).values())
    data_y_mock = np.sqrt(np.diag(covariance_tot[Pbin]))
    #"""
    
    error_diff = (data_y_kids-data_y_mock)/data_y_kids
    print(error_diff)
    print('Min: %g'%np.amin(error_diff))
    
    
    #print('GAMA/KiDS errors:', np.mean(data_y_gama[data_x_gama<30.])/np.mean(data_y_kids[data_x_kids<30.]))
    #print(np.sqrt(360.3/180))
    #print()
    
    # Plot error comparison
    fig = plt.figure(figsize=(5,4))

    plt.plot(data_x_kids, data_y_kids, label='KiDS (analytical covariance)', marker='.', ls='-', color=colors[0], zorder=3)
    #plt.plot(data_x_gama, data_y_gama, label='GAMA (analytical covariance)', marker='.', ls='-', color=colors[1], zorder=3)
    plt.plot(data_x_mock, data_y_mock, label='SLICS (mock covariance)', marker='.', ls='-', color=colors[3], zorder=3)


    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='best')

    xmin = 5
    xmax = 70
    plt.axvline(x=xmin, ls=':', color='black', zorder=2)
    plt.axvline(x=xmax, ls=':', color='black', zorder=2)
    
    # Define the labels for the plot
    plt.title(r'Radius: $\theta_{\rm A} = %g$ arcmin, Percentile bin: $%g < P < %g$'%(theta, percvals[Pbin], percvals[Pbin+1]))
    
    if 'pc' in Runit:
        plt.xlabel(r'Radius $R ({\rm Mpc} \, {h_{70}}^{-1})$', fontsize=14)
        plt.ylabel(r'ESD Amplitude [h$_{%g}$ M$_{\odot}$/pc$^2$]'%(h*100), fontsize=14)

    if Runit == 'arcmin':
        plt.xlabel(r'Angular separation $\theta$ (arcmin)', fontsize=14)
        plt.ylabel(r'Shear error $\sigma_{\gamma}$', fontsize=14)
        plt.xlim([2., 100])
        plt.axis([2.,100.,6e-5,1e-3])

    for ext in ['pdf']:

        plotfilename = '%s/Plots/error_comparison_Ptheta%s-%s_%s'%(path_filename_mock, theta, percnames[Pbin], percnames[Pbin+1])
        plotname = '%s.%s'%(plotfilename, ext)

        plt.savefig(plotname, format=ext, bbox_inches='tight')
        
    print('Written plot:', plotname)
    if show:
        plt.show()
    plt.close()
