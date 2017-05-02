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

#colors = ['blue', 'cyan', 'orange', 'red']
colors = ['#d7191c', '#fdae61', '#92c5de', '#0571b0']

#"""

# Redshifts
thetalist = np.array([10., 6.326]) # in arcmin
Runit = 'Mpc'

samples = ['gama_lowZ_complex', 'gama_highZ_complex']
samplelist = ['gama_lowZ_complex_10Mpc', 'gama_highZ_complex_6.326Mpc']
labels = [r'$0.1<z<0.197$', r'$0.197<z<0.3$']


"""

# Sizes
thetalist = np.array([5., 10., 15., 20.]) # in arcmin
Runit = 'arcmin'

samples = ['gama_absmag_complex']
samplelist = ['%s_%g%s'%(samples[0], theta, Runit) for theta in thetalist]

labels = np.array([r"$\theta_{\rm A} = %g'$"%theta for theta in thetalist])

"""



# Import amplitudes
path_filename = '/data2/brouwer/shearprofile/trough_results_Apr/Plots'
filenames = ['%s/trough_amplitudes_%s.txt'%(path_filename, sample) for sample in samplelist]
sample = samples[0]

amplitude_data = np.array([np.loadtxt(filename).T for filename in filenames])

perccenters = amplitude_data[:,0]
deltacenters = amplitude_data[:, 1]
Alist = amplitude_data[:,2]
Alist_error = amplitude_data[:,3]

#Alist = np.array([Alist[i]/thetalist[i] for i in range(len(thetalist))])
#Alist_error = np.array([Alist_error[i]/thetalist[i] for i in range(len(thetalist))])

## AMPLITUDE (Percentile)
fig = plt.figure(figsize=(5,4))

[plt.errorbar(perccenters[i], Alist[i], yerr=Alist_error[i],\
label=labels[i], marker='.', ls=':', color=colors[i]) for i in range(len(filenames))]

plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0.5, ls=':', color='black')

# Define the labels for the plot
plt.ylabel(r'Lensing detection Amplitude')
plt.xlabel(r"Percentile $P(\theta_{\rm A})$")

plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '%s/trough_amplitude_%s'%(path_filename, sample)
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()
plt.close()

## WEIGHT (Percentile)

# Calculating the weight corresponding to Amplitude/error
weightlist = Alist/Alist_error
weightlist_error = abs(1/Alist_error) * (Alist_error)

fig = plt.figure(figsize=(5,4))

model_x = np.linspace(0., 1., 100)
model_y = np.zeros([len(filenames), len(model_x)])
for i in range(len(filenames)):

    # Fitting a polynomial curve to the weight a.f.o. percentile
    poly_param_weights = np.polyfit(perccenters[i], weightlist[i], 5, w=1/weightlist_error[i])
    poly_func_weights = np.poly1d(poly_param_weights)

    model_y[i] = poly_func_weights(model_x)
    
    plt.plot(model_x, model_y[i], color=colors[i])
    plt.errorbar(perccenters[i], weightlist[i], yerr=weightlist_error[i],\
    label=labels[i], marker='.', ls='', color=colors[i])

plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0.5, ls=':', color='black')

# Define the labels for the plot
plt.ylabel(r'Lensing ``signal to noise" ($A/\delta A$)')
plt.xlabel(r"Percentile $P(\theta_{\rm A})$")

plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '%s/trough_weight_%s'%(path_filename, sample)
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()
plt.close()

## AMPLITUDE (delta)
fig = plt.figure(figsize=(5,4))

[plt.errorbar(deltacenters[i], Alist[i], yerr=Alist_error[i], \
label=labels[i], marker='.', ls=':', color=colors[i]) for i in range(len(filenames))]

plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0., ls=':', color='black')

# Save plot

# Define the labels for the plot
plt.ylabel(r'Lensing detection Amplitude')
plt.xlabel(r"Overdensity $\delta(\theta_{\rm A})$")

plt.legend(loc='upper left')

for ext in ['png', 'pdf']:

    plotfilename = '%s/trough_amplitude_delta_%s'%(path_filename, sample)
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()
plt.close()

# Import trough catalog
path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'

Ptheta = np.array([])

for s in range(len(samples)):
    
    troughcatname = 'trough_catalog_%s.fits'%samples[s]

    # Full directory & name of the trough catalogue
    troughcatfile = '%s/%s'%(path_troughcat, troughcatname)
    troughcat = pyfits.open(troughcatfile, memmap=True)[1].data
    print(troughcatfile)

    # Write weight fits-file
    if 'Z' in sample:
        Ptheta_z = troughcat['Ptheta%g'%thetalist[s]]
        Ptheta = np.vstack([Ptheta, Ptheta_z]) if Ptheta.size else Ptheta_z
    else:
        Ptheta = [troughcat['Ptheta%g'%theta] for theta in thetalist]

Wtheta = [poly_func_weights(Ptheta[t]) for t in range(len(thetalist))]

outputnames = ['Wtheta%g'%theta for theta in thetalist]
output = [np.abs(Wtheta[t]) for t in range(len(thetalist))]

sample = samples[0]
weightcatname = '%s/amplitude_trough_weights_%s.fits'%(path_troughcat, sample)
utils.write_catalog(weightcatname, np.arange(len(Wtheta[0])), outputnames, output)

