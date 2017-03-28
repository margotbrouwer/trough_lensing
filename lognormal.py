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
dperc = 0.05

# Troughs
percmin = 0.
percmax = 1.
percnames = ['0','0p05','0p1','0p15','0p2','0p25','0p3','0p35','0p4',\
'0p45','0p5','0p55','0p6','0p65','0p7','0p75','0p8','0p85','0p9','0p95','1']

perclist = np.arange(percmin, percmax, dperc)
perccenters = perclist+dperc/2.
perclist = np.append(perclist, percmax)

# Delta
deltalimits = (np.loadtxt('/data2/brouwer/shearprofile/KiDS-GAMA/brouwer/configs_margot/deltalimits.txt').T)[0]
deltacenters = (np.loadtxt('/data2/brouwer/shearprofile/KiDS-GAMA/brouwer/configs_margot/deltacenters.txt').T)[0]

deltanames = ['m1', 'm0p553201', 'm0p487318', 'm0p404268', 'm0p354624', 'm0p30498', 'm0p255335', 'm0p216329', 'm0p179709', 'm0p128441', \
            'm0p0771725', 'm0p0304349', '0p0205948', '0p0766321', '0p1279', '0p191464', '0p281705', '0p363309', '0p48933', '0p699793']#, '20p9428']
       
Npercs = 20
deltacenters = deltacenters[0:Npercs]

# Defining the paths to the data
blind = 'A'

path_sheardata = 'data2/brouwer/shearprofile/trough_results_Feb'
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins25_1_300_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

# Percentiles
path_lenssel = ['No_bins_gama_absmag/Pmasktheta5_0p8_1-Ptheta5_%s_%s'%(percnames[i], percnames[i+1]) for i in range(Npercs)]

# Delta
#path_lenssel = ['No_bins_gama_simplemask/Pmasktheta5_0p8_1-delta5_%s_%s'%(deltanames[i], deltanames[i+1]) for i in range(Npercs)]


esdfiles = np.array([('/%s/%s/%s/%s'%(path_sheardata, path_lenssel[i], path_cosmo, path_filename)) \
           for i in range(len(path_lenssel))])

covfiles = np.array([e.replace('bins_%s.txt'%blind, 'matrix_%s.txt'%blind) for e in esdfiles])

# Importing the shearprofiles and lens IDs
data_x, data_y, error_h, error_l = utils.read_esdfiles(esdfiles)
data_x = data_x[0]

# Define the part of the trough profile that contributes to the fit
xmin = 6.
xmax = 60.
xmask = (xmin < data_x) & (data_x < xmax)


# Plotting the ueber matrix
Nbins = Npercs
Nrows = 4
Ncolumns = 5#int(Nbins/Nrows)

fig = plt.figure(figsize=(15,9))
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

        A, Acov = optimization.curve_fit(trough_model, data_x[xmask], (data_y[N])[xmask], [0.], (error_h[N])[xmask])
        print('%g < P(x) < %g: Amplitude = %g'%(perclist[N], perclist[N+1], A))
        Alist = np.append(Alist, A)
        Alist_error = np.append(Alist_error, np.sqrt(Acov[0,0]))

        model_y = trough_model(data_x[xmask], A)

        ax_sub.errorbar(data_x, data_y[N], yerr=error_h[N], \
        ls='', marker='.', color='blue', label=r'$%g < P \leq %g$'%(perclist[N], perclist[N+1]))
        ax_sub.plot(data_x[xmask], model_y, ls='-', color='red')
        
        ax_sub.axvline(x=xmin, ls='--', color='black')
        ax_sub.axvline(x=xmax, ls='--', color='black')
        
        ax_sub.xaxis.set_label_position('bottom')
        ax_sub.yaxis.set_label_position('left')

        ax.tick_params(labelleft='off', labelbottom='off', top='off', bottom='off', left='off', right='off')
        
        # Last row
        if N1 == Nrows-1:
            plt.xlabel(r'Radial distance $R$ (arcmin)')
            plt.legend(loc='best')
        else:
            ax_sub.tick_params(axis='x', labelbottom='off')
            plt.legend(loc='best')

        # First column
        if N2 == 0:
            plt.ylabel(r'Shear $\gamma$')
        else:
            ax_sub.tick_params(axis='y', labelleft='off')
        
        plt.ylim(-0.0019, 0.0039)
        plt.xlim(1, 300)
        plt.xscale('log')

    
for ext in ['png', 'pdf']:

    plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/trough_amplitude_fits_delta'
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()


## AMPLITUDE
fig = plt.figure(figsize=(5,4))

# Fitting a polynomial curve to the Amplitude a.f.o. percentile
poly_param_amplitude = np.polyfit(deltacenters, Alist, 2, w=1/Alist_error)
poly_func_amplitude = np.poly1d(poly_param_amplitude)
model_y = poly_func_amplitude(deltacenters)

# Plot best-fit polynomial
plt.plot(deltacenters, model_y, label=r'Best-fit $2^{\rm th}$ degree polynomial', color='red')

plt.errorbar(deltacenters, Alist, yerr=Alist_error, marker='.', ls='', color='blue')
plt.axhline(y=0., ls='--', color='black')
plt.axvline(x=0., ls='--', color='black')


# Save plot

# Define the labels for the plot
plt.ylabel(r'Lensing detection Amplitude')

#plt.xlabel(r'Percentile $P(\theta)$')
plt.xlabel(r'Overdensity $\delta$')

plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/trough_amplitudes_delta'
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()

## WEIGHT
fig = plt.figure(figsize=(5,4))

# Calculating the weight corresponding to Amplitude/error
weightlist = Alist/Alist_error
weightlist_error = abs(1/Alist_error) * (Alist_error)


# Fitting a polynomial curve to the weight a.f.o. percentile
poly_param_weights = np.polyfit(deltacenters, weightlist, 6, w=1/weightlist_error)
poly_func_weights = np.poly1d(poly_param_weights)
model_y = poly_func_weights(deltacenters)

# Plot best-fit polynomial
plt.plot(deltacenters, model_y, label=r'Best-fit $6^{\rm th}$ degree polynomial', color='red')

print(Alist_error)

# Save plot
plt.errorbar(deltacenters, weightlist, yerr=weightlist_error, marker='.', ls='', color='blue')

plt.axhline(y=0., ls='--', color='black')
plt.axvline(x=0., ls='--', color='black')

# Define the labels for the plot
plt.ylabel(r'Lensing detection weight ($A/\delta A$)')
#plt.xlabel(r'Percentile $P(\theta)$')
plt.xlabel(r'Overdensity $\delta$')


# Save plot
plt.axhline(y=0., ls='--', color='black')
plt.axvline(x=0., ls='--', color='black')

# Define the labels for the plot
plt.ylabel(r'Lensing detection weight ($A/\delta A$)')


plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/trough_weights_delta'
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()



# Write weight fits-file

# Import trough catalog
path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogues'
troughcatname = '/trough_catalog_gama_absmag_masked.fits'

# Full directory & name of the trough catalogue
troughcatfile = '%s/%s'%(path_troughcat, troughcatname)
troughcat = pyfits.open(troughcatfile, memmap=True)[1].data

# List of the observables of all sources in the KiDS catalogue
Ptheta = troughcat['Ptheta5']
Atheta = poly_func_amplitude(Ptheta)
Wtheta = poly_func_weights(Ptheta)

outputnames = ['Ptheta5', 'Atheta5', 'Wtheta5']
output = [Ptheta, Atheta, Wtheta]
utils.write_catalog('%s/amplitude_trough_weights.fits'%path_troughcat, np.arange(len(Wtheta)), outputnames, output)
