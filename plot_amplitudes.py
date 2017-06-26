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
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import scipy.optimize as optimization
import trough_modules_all as utils

# Make use of TeX
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})

#colors = ['red', 'orange', 'cyan', 'blue']
colors = ['#d7191c', '#fdae61', '#92c5de', '#0571b0']

h=0.7

#"""

# Redshifts
thetalist = np.array([10., 6.326]) # in arcmin
Runit = 'Mpc'

selection = ['kids_lowZ_complex', 'kids_highZ_complex']
#selection = ['mice_lowZ_nomask', 'mice_highZ_nomask']
mocksel = ['mice_lowZ_nomask', 'mice_highZ_nomask']

labels = [r'$0.1<z<0.197$', r'$0.197<z<0.3$']


"""

# Sizes
thetalist = np.array([5., 10., 15., 20.]) # in arcmin
#thetalist = np.array([5.]) # in arcmin
Runit = 'arcmin'

selection = ['kids_mice_complex' for t in range(len(thetalist))]
mocksel = ['mice_all_nomask' for t in range(len(thetalist))]

labels = np.array([r"$\theta_{\rm A} = %g'$"%theta for theta in thetalist])

"""


# Import observed amplitudes
path_filename = '/data2/brouwer/shearprofile/trough_results_June/Plots'
filenames = ['%s/trough_amplitudes_%s_%g%s.txt'%(path_filename, selection[s], thetalist[s], Runit) \
                                                for s in range(len(selection))]
selection_name = selection[0]

amplitude_data = np.array([np.loadtxt(filename).T for filename in filenames])

perccenters = amplitude_data[:,0]
deltacenters = amplitude_data[:, 1]
Alist = amplitude_data[:,2]
Alist_error = amplitude_data[:,3]

#Alist = np.array([Alist[i]/thetalist[i] for i in range(len(thetalist))])
#Alist_error = np.array([Alist_error[i]/thetalist[i] for i in range(len(thetalist))])


# Import mock amplitudes
path_filename_mock = path_filename
filenames_mock = ['%s/trough_amplitudes_%s_%g%s.txt'%(path_filename_mock, mocksel[s], thetalist[s], Runit) \
                                                for s in range(len(mocksel))]

amplitude_data_mock = np.array([np.loadtxt(filename).T for filename in filenames_mock])

perccenters_mock = amplitude_data_mock[:,0]
deltacenters_mock = amplitude_data_mock[:, 1]
Alist_mock = amplitude_data_mock[:,2]
Alist_error_mock = amplitude_data_mock[:,3]


## AMPLITUDE (Percentile)
fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot(111)

model_x = np.linspace(0., 1., 100)
model_y = np.zeros([len(filenames), len(model_x)])
for i in range(len(filenames)):
    
    # Plot observed amplitudes
    plt.errorbar(perccenters[i], Alist[i], yerr=Alist_error[i],\
    label=labels[i], marker='.', ls='', color=colors[i])
    
    if ('kids' in selection_name) or ('gama' in selection_name):
        # Plot mock amplitudes
        plt.plot(perccenters_mock[i], Alist_mock[i], ls='--', color=colors[i])
        
    # Fitting a polynomial curve to the amplitude a.f.o. percentile
    poly_param_amps = np.polyfit(perccenters[i], Alist[i], 5, w=1/Alist_error[i])
    poly_func_amps = np.poly1d(poly_param_amps)
    
    #model_y[i] = poly_func_amps(model_x)
    #plt.plot(model_x, model_y[i], color=colors[i])

plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0.5, ls=':', color='black')

# Define the labels for the plot
if 'pc' in Runit:
    plt.ylabel(r'ESD Amplitude [h$_{%g}$ M$_{\odot}$/pc$^2$]'%(h*100))
    plt.axis([0.,1.,-6.,15.])
if 'arcmin' in Runit:
    plt.ylabel(r'Shear Amplitude')
    plt.axis([0.,1.,-0.008,0.012])
    
plt.xlabel(r"Density percentile $P(\theta_{\rm A})$")

plt.legend(loc='best')

if 'arcmin' in Runit:
    # Include zoom-in of center
    
    left, bottom, width, height = [0.555, 0.16, 0.31, 0.23]
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    ax2 = fig.add_axes([left, bottom, width, height])

    # Plot observed amplitudes
    [plt.errorbar(perccenters[i], Alist[i], yerr=Alist_error[i], \
    marker='.', ls='', color=colors[i]) for i in range(len(filenames))]
    
    if ('kids' in selection_name) or ('gama' in selection_name):
        # Plot mock amplitudes
        [plt.plot(perccenters_mock[i], Alist_mock[i], ls='--', color=colors[i]) \
        for i in range(len(filenames))]

    plt.axhline(y=0., ls=':', color='black')
    plt.axvline(x=0.5, ls=':', color='black')

    ax2.axis([0.45,0.65,-0.002,0.002])
    ax2.tick_params(labelleft='off', top='off', left='off', right='off')

for ext in ['png', 'pdf']:

    plotfilename = '%s/trough_amplitude_%s'%(path_filename, selection_name)
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
plt.ylabel(r'Shear signal to noise ($A/\delta A$)')
plt.xlabel(r"Density percentile $P(\theta_{\rm A})$")

plt.legend(loc='best')

for ext in ['png', 'pdf']:

    plotfilename = '%s/trough_weight_%s'%(path_filename, selection_name)
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()
plt.close()


## AMPLITUDE (delta)
fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot(111)

# Plot observed amplitudes
[plt.errorbar(deltacenters[i], Alist[i], yerr=Alist_error[i], \
label=labels[i], marker='.', ls='', color=colors[i]) for i in range(len(filenames))]

if ('kids' in selection_name) or ('gama' in selection_name):
    # Plot mock amplitudes
    [plt.plot(deltacenters_mock[i], Alist_mock[i], ls='--', color=colors[i]) \
    for i in range(len(filenames))]

plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0., ls=':', color='black')

# Save plot

# Define the labels for the plot
if 'pc' in Runit:
    plt.ylabel(r'ESD Amplitude [h$_{%g}$ M$_{\odot}$/pc$^2$]'%(h*100))
    plt.axis([-1.,2.,-10.,15.])
if 'arcmin' in Runit:
    plt.ylabel(r'Shear Amplitude')
    plt.axis([-0.8,1.5,-0.007,0.012])

plt.xlabel(r"Overdensity $\delta(\theta_{\rm A})$")

plt.legend(loc='upper left')

# Include zoom-in of center
# These are in unitless percentages of the figure size. (0,0 is bottom left)
left, bottom, width, height = [0.48, 0.16, 0.39, 0.25]
ax2 = fig.add_axes([left, bottom, width, height])

# Plot observed amplitudes
[plt.errorbar(deltacenters[i], Alist[i], yerr=Alist_error[i], \
label=labels[i], marker='.', ls='', color=colors[i]) for i in range(len(filenames))]

if ('kids' in selection_name) or ('gama' in selection_name):
    # Plot mock amplitudes
    [plt.plot(deltacenters_mock[i], Alist_mock[i], ls='--', color=colors[i]) \
    for i in range(len(filenames))]

plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0., ls=':', color='black')


if 'pc' in Runit:
    ax2.axis([-0.2,0.2,-2.,2.])
if 'arcmin' in Runit:
    ax2.axis([-0.2,0.2,-0.002,0.002])

ax2.tick_params(labelleft='off', top='off', left='off', right='off')
plt.xticks(np.arange(-0.2,0.3,0.1))

for ext in ['png', 'pdf']:

    plotfilename = '%s/trough_amplitude_delta_%s'%(path_filename, selection_name)
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()
plt.close()

# Import trough catalog
path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'

Ptheta = np.array([])

for s in range(len(selection)):
    
    troughcatname = 'trough_catalog_%s.fits'%selection[s]

    # Full directory & name of the trough catalogue
    troughcatfile = '%s/%s'%(path_troughcat, troughcatname)
    troughcat = pyfits.open(troughcatfile, memmap=True)[1].data
    print()
    print(troughcatfile)

    # Write weight fits-file
    if 'Z' in selection_name:
        Ptheta_z = troughcat['Ptheta%g'%thetalist[s]]
        Ptheta = np.vstack([Ptheta, Ptheta_z]) if Ptheta.size else Ptheta_z
    else:
        Ptheta = [troughcat['Ptheta%g'%theta] for theta in thetalist]

Wtheta = [poly_func_weights(Ptheta[t]) for t in range(len(thetalist))]
Atheta = [poly_func_amps(Ptheta[t]) for t in range(len(thetalist))]


outputnames = []
[outputnames.append('Wtheta%g'%theta) for theta in thetalist]
[outputnames.append('Atheta%g'%theta) for theta in thetalist]

output = []
[output.append(np.abs(Wtheta[t])) for t in range(len(thetalist))]
[output.append(np.abs(Atheta[t])) for t in range(len(thetalist))]

weightcatname = '%s/amplitude_trough_weights_%s.fits'%(path_troughcat, selection_name)
utils.write_catalog(weightcatname, np.arange(len(Wtheta[0])), outputnames, output)

