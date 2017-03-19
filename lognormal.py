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

# Model to fit the troughs
def trough_model(x, A):
    
    model_y = A/x
    
    return model_y

# Defining the percentile bins
dperc = 0.05

# Troughs
percmin = 0.
percmax = 1.
percnames = ['0','0p05','0p1','0p15','0p2','0p25','0p3','0p35','0p4','0p45','0p5','0p55','0p6','0p65','0p7','0p75','0p8','0p85','0p9','0p95','1']

perclist = np.arange(percmin, percmax, dperc)
perccenters = perclist+dperc/2.
perclist = np.append(perclist, percmax)
Npercs = len(perccenters)

# Defining the paths to the data
blind = 'A'

path_sheardata = 'data2/brouwer/shearprofile/trough_results_Feb'
path_lenssel = ['No_bins_gama_absmag/Pmasktheta5_0p8_1-Ptheta5_%s_%s'%(percnames[i], percnames[i+1]) for i in range(Npercs)]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins25_1_300_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

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

fig = plt.figure(figsize=(12,6))
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
        ls='', marker='.', color='blue', label=r'$%g < P(x) < %g$'%(perclist[N], perclist[N+1]))
        ax_sub.plot(data_x[xmask], model_y, ls='-', color='red')
        
        ax_sub.axvline(x=xmin, ls='--')
        ax_sub.axvline(x=xmax, ls='--')
        
        ax_sub.xaxis.set_label_position('top')
        ax_sub.yaxis.set_label_position('right')

        ax.tick_params(labelleft='off', labelbottom='off', top='off', bottom='off', left='off', right='off')
                
        if N2 != 0:
            ax_sub.tick_params(axis='y', labelleft='off')
        
        plt.ylim(-0.003, 0.003)
        plt.xscale('log')


for ext in ['png', 'pdf']:

    plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/trough_amplitude_fits'
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()




# Fitting a polynomial curve to the Amplitude a.f.o. percentile
poly_param_amplitude = np.polyfit(perccenters, Alist, 6, w=1/Alist_error)
poly_func_amplitude = np.poly1d(poly_param_amplitude)
model_y = poly_func_amplitude(perccenters)

# Plot best-fit polynomial
plt.plot(perccenters, model_y, label=r'Best-fit $6^{\rm th}$ degree polynomial')

print(Alist_error)

plt.errorbar(perccenters, Alist, yerr=Alist_error, marker='.', ls='')
plt.axhline(y=0., ls='--')
plt.axvline(x=0.5, ls='--')

# Save plot

# Define the labels for the plot
plt.xlabel(r'Trough percentile $P(\theta)$')
plt.ylabel(r'Lensing detection Amplitude')

plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/trough_amplitudes'
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()


# Calculating the weight corresponding to Amplitude/error
weightlist = Alist/Alist_error
weightlist_error = abs(1/Alist_error) * (Alist_error)
print(weightlist_error)

# Fitting a polynomial curve to the weight a.f.o. percentile
poly_param_weights = np.polyfit(perccenters, weightlist, 3, w=1/weightlist_error)
poly_func_weights = np.poly1d(poly_param_weights)
model_y = poly_func_weights(perccenters)

# Plot best-fit polynomial
plt.plot(perccenters, model_y, label=r'Best-fit $3^{\rm rd}$ degree polynomial')

print(Alist_error)

# Save plot
plt.errorbar(perccenters, weightlist, yerr=weightlist_error, marker='.', ls='')
plt.axhline(y=0., ls='--')
plt.axvline(x=0.5, ls='--')

# Define the labels for the plot
plt.xlabel(r'Trough percentile $P(\theta)$')
plt.ylabel(r'Lensing detection weight ($A/\delta A$)')

plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/trough_weights'
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
