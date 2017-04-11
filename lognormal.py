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

# Make use of TeX
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})


# Model to fit the troughs
def trough_model(x, A):
    
    model_y = A/x
    
    return model_y

# Defining the percentile bins
percmin = 0.
percmax = 1.

#dperc = 0.05
dperc = 0.1

#percnames = ['0','0p05','0p1','0p15','0p2','0p25','0p3','0p35','0p4',\
#'0p45','0p5','0p55','0p6','0p65','0p7','0p75','0p8','0p85','0p9','0p95','1']
percnames = ['0','0p1','0p2','0p3','0p4','0p5','0p6','0p7','0p8','0p9','1']


perclist = np.arange(percmin, percmax, dperc)
perccenters = perclist+dperc/2.
perclist = np.append(perclist, percmax)

# Delta
deltalimits = (np.loadtxt('/data2/brouwer/shearprofile/KiDS-GAMA/brouwer/configs_margot/deltalimits.txt').T)[0]
deltacenters = (np.loadtxt('/data2/brouwer/shearprofile/KiDS-GAMA/brouwer/configs_margot/deltacenters.txt').T)[0]

deltanames = ['m1', 'm0p553201', 'm0p487318', 'm0p404268', 'm0p354624', 'm0p30498', 'm0p255335', 'm0p216329', 'm0p179709', 'm0p128441', \
            'm0p0771725', 'm0p0304349', '0p0205948', '0p0766321', '0p1279', '0p191464', '0p281705', '0p363309', '0p48933', '0p699793']#, '20p9428']
       
Npercs = len(percnames)-1
deltacenters = deltacenters[0:Npercs]

# Defining the paths to the data
blind = 'A'
sample = 'gama_lowZ_complex'
Runit = 'Mpc' # arcmin or Mpc

thetalist = np.array([5., 10., 15., 20.]) # in arcmin
if 'highZ' in sample:
    thetalist = np.array([3.163, 6.326, 9.490, 12.653])
Rlist = [0.94126266, 1.8825293, 2.8238039, 3.76509045] # Physical size of the troughs (in Mpc)

thetanum = 1
theta = thetalist[thetanum]

if Runit == 'arcmin':
    Rmin = 2
    Rmax = 100

if Runit == 'Mpc':
    Rmin = 0.5
    Rmax = 20

"""
path_sheardata = 'data2/brouwer/shearprofile/trough_results_Feb'
path_lenssel = ['No_bins_gama_absmag/Pmasktheta5_0p8_1-Ptheta5_%s_%s'%(percnames[i], percnames[i+1]) for i in range(Npercs)]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins25_1_300_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)
"""

# Redshift
path_sheardata = 'data2/brouwer/shearprofile/trough_results_Apr'
path_lenssel = ['No_bins_%s/Pmasktheta%s_0p8_1-Ptheta%s_%s_%s'\
    %(sample, ('%g'%theta).replace('.','p'), ('%g'%theta).replace('.','p'), percnames[i], percnames[i+1]) for i in range(Npercs)]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins10_%s_%s_%s/shearcovariance'%(('%g'%Rmin).replace('.','p'), ('%g'%Rmax).replace('.','p'), Runit)
path_filename = 'No_bins_%s.txt'%(blind)

path_plots = '/data2/brouwer/shearprofile/trough_results_Apr/Plots/%s'%sample

esdfiles = np.array([('/%s/%s/%s/%s'%(path_sheardata, path_lenssel[i], path_cosmo, path_filename)) \
           for i in range(len(path_lenssel))])

covfiles = np.array([e.replace('bins_%s.txt'%blind, 'matrix_%s.txt'%blind) for e in esdfiles])

# Importing the shearprofiles and lens IDs
data_x, data_y, error_h, error_l = utils.read_esdfiles(esdfiles)

covariance_tot = np.array([ np.loadtxt(covfiles[c]).T for c in range(len(covfiles)) ])

# Import random signal
path_randoms = 'No_bins_randoms/Pmasktheta5_0_inf'
random_esdfile = np.array(['/%s/%s/%s/%s'%(path_sheardata, path_randoms, path_cosmo, path_filename)])

random_data_x, random_data_y, random_error_h, random_error_l = utils.read_esdfiles(random_esdfile)
random_data_x, random_data_y, random_error_h, random_error_l = random_data_x[0], random_data_y[0], random_error_h[0], random_error_l[0]

# Subtract random signal
data_x = data_x[0]
data_y = data_y-random_data_y
error_h = np.sqrt(error_h**2. + random_error_h**2)
error_l = np.sqrt(error_l**2. + random_error_l**2)


# Define the part of the trough profile that contributes to the fit
if Runit == 'arcmin':
    xmin = theta*1.2
    xmax = 70.
if Runit == 'Mpc':
    xmin = Rlist[thetanum]*1.2
    xmax = 10

xmask = (xmin < data_x) & (data_x < xmax)
xwhere = np.where(xmask)[0]

# Plotting the ueber matrix
Nbins = Npercs
Nrows = 2
Ncolumns = 5#int(Nbins/Nrows)

fig = plt.figure(figsize=(10,6))
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
        ax_sub = fig.add_subplot(gs[N1, N2])
      
        # Without covariance
        A, Acov = optimization.curve_fit(f=trough_model, xdata=data_x[xmask], ydata=(data_y[N])[xmask], p0=[0.], \
        sigma=(error_h[N])[xmask], absolute_sigma=True)
        
        """
        # With covariance
        covariance = covariance_tot[N]
        ind = np.lexsort((covariance[3,:], covariance[1,:], covariance[2,:], covariance[0,:]))
        covmatrix = np.reshape(covariance[4][ind], [len(data_x), len(data_x)])
        covmatrix = covmatrix[int(xwhere[0]):int(xwhere[-1]+1), int(xwhere[0]):int(xwhere[-1]+1)]
        covmatrix = np.array([i for i in covmatrix], dtype=float)
        
        A, Acov = optimization.curve_fit(f=trough_model, xdata=data_x[xmask], ydata=(data_y[N])[xmask], p0=[0.], \
        sigma=covmatrix, absolute_sigma=True)
        """
        
        print('%g < P(x) < %g: Amplitude = %g'%(perclist[N], perclist[N+1], A))
        
        Alist = np.append(Alist, A)
        Alist_error = np.append(Alist_error, np.sqrt(Acov[0,0]))

        model_x = np.linspace(xmin, xmax, 10)
        model_y = trough_model(model_x, A)

        ax_sub.plot(model_x, model_y, ls='-', color='red', label=r'$%g < P \leq %g$'%(perclist[N], perclist[N+1]), zorder=4)
        
        ax_sub.errorbar(data_x, data_y[N], yerr=error_h[N], ls='', marker='.', color='blue', zorder=3)
        
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

ax.set_xlabel(r'Radial distance $R$ (arcmin)', fontsize=14)
ax.set_ylabel(r'Shear $\gamma$', fontsize=14)

ax.xaxis.set_label_coords(0.5, -0.07)
ax.yaxis.set_label_coords(-0.07, 0.5)

for ext in ['png', 'pdf']:

    plotfilename = '%s_trough_amplitude_fits'%path_plots
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()


## AMPLITUDE (Percentile)
fig = plt.figure(figsize=(4,3))

# Fitting a polynomial curve to the Amplitude a.f.o. percentile
poly_param_amplitude = np.polyfit(perccenters, Alist, 5, w=1/Alist_error)
poly_func_amplitude = np.poly1d(poly_param_amplitude)
model_y = poly_func_amplitude(perccenters)

# Plot best-fit polynomial
#plt.plot(perccenters, model_y, label=r'Best-fit $5^{\rm th}$ degree polynomial', color='red')

plt.errorbar(perccenters, Alist, yerr=Alist_error, marker='.', ls='', color='blue')
plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0.5, ls=':', color='black')

Afilename = '/%s/Plots/trough_amplitudes_%s_%g%s.txt'%(path_sheardata, sample, theta, Runit)
np.savetxt(Afilename, np.array([perccenters, Alist, Alist_error]).T, header = 'Trough percentile     Amplitude     Error(Amplitude)')
print('Written textfile:', Afilename)

# Save plot

# Define the labels for the plot
plt.ylabel(r'Lensing detection Amplitude')
plt.xlabel(r"Percentile $P(\theta=%g')$"%theta)

plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '%s_trough_amplitudes'%path_plots
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()



## AMPLITUDE (Delta)
fig = plt.figure(figsize=(4,3))

"""
# Fitting a polynomial curve to the Amplitude a.f.o. percentile
poly_param_amplitude = np.polyfit(deltacenters, Alist, 2, w=1/Alist_error)
poly_func_amplitude = np.poly1d(poly_param_amplitude)
model_y = poly_func_amplitude(deltacenters)

# Plot best-fit polynomial
plt.plot(deltacenters, model_y, label=r'Best-fit $2^{\rm nd}$ degree polynomial', color='red')
"""

plt.errorbar(deltacenters, Alist, yerr=Alist_error, marker='.', ls='', color='blue')
plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0., ls=':', color='black')


# Save plot

# Define the labels for the plot
plt.ylabel(r'Lensing detection Amplitude')
plt.xlabel(r"Overdensity $\delta(\theta=%g')$"%theta)

plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '%s_trough_amplitudes_delta'%path_plots
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()


## WEIGHT
fig = plt.figure(figsize=(4,3))

# Calculating the weight corresponding to Amplitude/error
weightlist = Alist/Alist_error
weightlist_error = abs(1/Alist_error) * (Alist_error)


# Fitting a polynomial curve to the weight a.f.o. percentile
poly_param_weights = np.polyfit(perccenters, weightlist, 4, w=1/weightlist_error)
poly_func_weights = np.poly1d(poly_param_weights)
model_x = np.linspace(0., 1., 100)
model_y = poly_func_weights(model_x)

# Plot best-fit polynomial
plt.plot(model_x, model_y, label=r'Best-fit $4^{\rm th}$ degree polynomial', color='red')

# Save plot
plt.errorbar(perccenters, weightlist, yerr=weightlist_error, marker='.', ls='', color='blue')

plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0.5, ls=':', color='black')

# Define the labels for the plot
plt.ylabel(r'Lensing ``signal to noise" ($A/\delta A$)')

plt.xlabel(r"Percentile $P(\theta=%g')$"%theta)
#plt.xlabel(r'Overdensity $\delta$')

plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '%s_trough weights'%path_plots
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()



# Write weight fits-file

# Import trough catalog
path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'
troughcatname = 'trough_catalog_%s_%s_0.04deg_%s.fits'%(sample.split('_')[0], sample.split('_')[1], sample.split('_')[2])

# Full directory & name of the trough catalogue
troughcatfile = '%s/%s'%(path_troughcat, troughcatname)
troughcat = pyfits.open(troughcatfile, memmap=True)[1].data

# List of the observables of all sources in the KiDS catalogue
Ptheta = troughcat['Ptheta%g'%theta]
Atheta = poly_func_amplitude(Ptheta)
Wtheta = poly_func_weights(Ptheta)

outputnames = ['Ptheta%g'%theta, 'Atheta%g'%theta, 'Wtheta%g'%theta]
output = [Ptheta, Atheta, Wtheta]

weightcatname = '%s/amplitude_trough_weights_%s_%g%s.fits'%(path_troughcat, sample, theta, Runit)
utils.write_catalog(weightcatname, np.arange(len(Wtheta)), outputnames, output)
