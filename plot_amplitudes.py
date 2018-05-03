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

import matplotlib.style
matplotlib.style.use('classic')

from scipy.stats import chi2, norm
import scipy.optimize as optimization
import trough_modules_all as utils

# Make use of TeX
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})

#colors = ['red', 'orange', 'cyan', 'blue', 'green']
colors = ['#d7191c', '#fdae61', '#92c5de', '#0571b0', '#0ba52d']

h=0.7

show = True

# Trough selection

"""
"""

# Mice redshifts

thetalist = np.array([20., 12.85, 9.45, 7.44, 6.14]) # Dc
zlims = np.array([0.1, 0.191, 0.286, 0.385, 0.489, 0.6])      


Runit = 'Mpc'
valpha = 0.3
linestyle = '-'

selection = ['mice_miceZ-%i_nomask-Z'%(t+1) for t in range(len(thetalist))]
#selection = ['mice_miceZa_nomask-Za-%i'%(t+1) for t in range(len(thetalist))]
selection_name = selection[0]

labels = [r'$%g<z<%g$'%(zlims[t], zlims[t+1]) for t in range(len(thetalist))]

"""

# lowZ/highZ

thetalist = np.array([10., 6.303]) # in arcmin

Runit = 'Mpc'
valpha = 0.3
linestyle = '-'

selection = ['kids_lowZ_complex', 'kids_highZ_complex']
#selection = ['mice_lowZ_nomask-1', 'mice_highZ_nomask-1']
selection_name = selection[0]

zlims = np.array([ 0.1, 0.198, 0.3])
labels = [r'$%g<z<%g$'%(zlims[t], zlims[t+1]) for t in range(len(thetalist))]

#mocksel = ['mice_lowZ_nomask-1', 'mice_highZ_nomask-1']
#mockthetalist = thetalist
#mockcolors = colors

mocksel = np.append(['mice_lowZ_nomask-%g'%ij for ij in np.arange(16)+1.], ['mice_highZ_nomask-%g'%ij for ij in np.arange(16)+1.])
mockthetalist = np.append( [10.]*(len(mocksel)/2), [6.303]*(len(mocksel)/2) ) # in arcmin
mockcolors = np.append( ['#d7191c']*(len(mocksel)/2), ['#fdae61']*(len(mocksel)/2) )



# Sizes
thetalist = np.array([5., 10., 15., 20.]) # in arcmin

Runit = 'arcmin'
valpha = 1.
linestyle = '--'

selection = ['kids_mice_complex' for t in range(len(thetalist))]
selection_name = selection[0]

labels = np.array([r"$\theta_{\rm A} = %g'$"%theta for theta in thetalist])

#mocksel = ['mice_all_nomask' for t in range(len(thetalist))]
mocksel = ['slics_mocks_nomask' for t in range(len(thetalist))]

mockthetalist = thetalist
mockcolors = colors



# SLICS mock redshifts

thetalist = np.array([15., 9.554, 7.283, 5.770]) # Dc
zlims = np.array( [0.1, 0.192, 0.289, 0.391, 0.5])

Runit = 'Mpc'
valpha = 0.3
linestyle = '-'
labels = np.array([r"$\theta_{\rm A} = %g'$"%theta for theta in thetalist])

selection = ['slics_mockZ_nomask' for t in range(len(thetalist))]
selection_name = selection[0]

labels = [r'$%g<z<%g$'%(zlims[t], zlims[t+1]) for t in range(len(thetalist))]



"""


# Import observed amplitudes
path_filename = '/data2/brouwer/shearprofile/trough_results_final/Plots'
filenames = ['%s/trough_amplitudes_%s_%g%s.txt'%(path_filename, selection[s], thetalist[s], Runit) \
                                                for s in range(len(selection))]
print('Imported:', filenames)

amplitude_data = np.array([np.loadtxt(filename).T for filename in filenames])

perccenters = amplitude_data[:,0]
deltacenters = amplitude_data[:, 1]
Alist = amplitude_data[:,2]
Alist_error = amplitude_data[:,3]

#Alist = np.array([Alist[i]/thetalist[i] for i in range(len(thetalist))])
#Alist_error = np.array([Alist_error[i]/thetalist[i] for i in range(len(thetalist))])

if ('kids' in selection_name) or ('gama' in selection_name):
    # Import mock amplitudes
    path_filename_mock = '/data2/brouwer/shearprofile/trough_results_final/Plots'
    filenames_mock = ['%s/trough_amplitudes_%s_%g%s.txt'%(path_filename_mock, mocksel[s], mockthetalist[s], Runit) \
                                                    for s in range(len(mocksel))]

    amplitude_data_mock = np.array([np.loadtxt(filename).T for filename in filenames_mock])

    perccenters_mock = amplitude_data_mock[:,0]
    deltacenters_mock = amplitude_data_mock[:, 1]
    Alist_mock = amplitude_data_mock[:,2]
    Alist_error_mock = amplitude_data_mock[:,3]

if 'lowZ' in selection_name:

    # Compute significance of the difference

    diffAlist = abs(Alist[0] - Alist[1])
    diffAlist_error = np.sqrt(Alist_error[0]**2. + Alist_error[1]**2.)
    diffmodel = np.zeros(len(Alist[0]))

    diffchi2 = np.sum( (diffAlist - diffmodel)**2. / diffAlist_error**2. )
    diffprob = chi2.cdf(diffchi2, len(diffAlist)-1)
    diffsigma = norm.ppf( diffprob + (1.-diffprob)/2. )

    print('chi^2:', 'Prob:', 'Sigma')
    print(diffchi2, diffprob, diffsigma)

if 'mockZ' in selection_name:

    # Compute significance of the difference
    for x in range(len(Alist)-1):
        diffAlist = abs(Alist[x] - Alist[x+1])
        diffAlist_error = np.sqrt(Alist_error[x]**2. + Alist_error[x+1]**2.)
        diffmodel = np.zeros(len(Alist[0]))

        diffchi2 = np.sum( (diffAlist - diffmodel)**2. / diffAlist_error**2. )
        diffprob = chi2.cdf(diffchi2, len(diffAlist)-1)
        diffsigma = norm.ppf( diffprob + (1.-diffprob)/2. )

        print('chi^2:', 'Prob:', 'Sigma')
        print(diffchi2, diffprob, diffsigma)

    diffAlist = abs(Alist[0] - Alist[len(thetalist)-1])
    diffAlist_error = np.sqrt(Alist_error[0]**2. + Alist_error[len(thetalist)-1]**2.)
    diffmodel = np.zeros(len(Alist[0]))

    diffchi2 = np.sum( (diffAlist - diffmodel)**2. / diffAlist_error**2. )
    diffprob = chi2.cdf(diffchi2, len(diffAlist)-1)
    diffsigma = norm.ppf( diffprob + (1.-diffprob)/2. )

    print('chi^2:', 'Prob:', 'Sigma')
    print(diffchi2, diffprob, diffsigma)

## AMPLITUDE (Percentile)
fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot(111)

Npoly = 5
model_x = np.linspace(0., 1., 100)
poly_param_amps = []

for i in range(len(selection)):
    # Fitting a polynomial curve to the amplitude a.f.o. percentile
    poly_param_amps.append( np.polyfit(perccenters[i], Alist[i], Npoly, w=1/Alist_error[i]) )
    poly_func_amps = np.poly1d(poly_param_amps[i])
  
if ('kids' in selection_name) or ('gama' in selection_name):
    # Plot mock amplitudes
    [plt.plot(perccenters_mock[i], Alist_mock[i], ls=linestyle, color=mockcolors[i], alpha=valpha, zorder=1) \
                for i in range(len(mocksel))]

# Plot observed amplitudes
if ('kids' in selection_name) or ('gama' in selection_name):
    
    if 'Z' in selection_name:
        [plt.errorbar(perccenters[i], Alist[i], yerr=Alist_error[i],\
        label=labels[i], marker='o', ms=5, ls='', color=colors[i], ecolor='black', zorder=3) \
        for i in range(len(selection))]
    else:
        [plt.errorbar(perccenters[i], Alist[i], yerr=Alist_error[i],\
        label=labels[i], marker='.', ls='', color=colors[i], zorder=3)
        for i in range(len(selection))]
else:
    if 'miceZ' in selection_name:
        [plt.plot(perccenters[i], Alist[i], label=labels[i],  color=colors[i], \
            marker='', ls='-', alpha=1., zorder=3) for i in range(len(selection))]
    if 'slics' in selection_name:
        [plt.errorbar(perccenters[i], Alist[i], yerr=Alist_error[i], label=labels[i], color=colors[i], \
            marker='.', ls='--', alpha=1., zorder=3) for i in range(len(selection))]
        #[plt.fill_between(perccenters[i], Alist[i]-Alist_error[i]/2., Alist[i]+Alist_error[i]/2., \
        #color=colors[i], alpha=valpha)  for i in range(len(selection))]

plt.axhline(y=0., ls=':', color='black', zorder=2)
plt.axvline(x=0.5, ls=':', color='black', zorder=2)

# Define the labels for the plot
if 'pc' in Runit:
    plt.ylabel(r'Comoving ESD Amplitude [h$_{%g}^{1/2}$ M$_{\odot}$/pc$^{3/2}$]'%(h*100))
    plt.axis([0.,1.,-3.,9.])
    if 'miceZ' in selection_name:
        plt.axis([0.,1.,-5.,9.])
    if 'slics' in selection_name:
        plt.axis([0.,1.,-3.,5.])
if 'arcmin' in Runit:
    plt.ylabel(r'Shear Amplitude [${\rm arcmin}^{1/2}$]')
    plt.axis([0.,1.,-0.008,0.012])
    
plt.xlabel(r"Density percentile rank $P(\theta_{\rm A})$")

plt.legend(loc='best')

if 'arcmin' in Runit:
    # Include zoom-in of center
    
    left, bottom, width, height = [0.555, 0.16, 0.31, 0.23]
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    ax2 = fig.add_axes([left, bottom, width, height])

    # Plot observed amplitudes
    if 'Z' in selection_name:
        [plt.errorbar(perccenters[i], Alist[i], yerr=Alist_error[i], \
        marker='o', ms=5, ls='', color=colors[i], ecolor='black', zorder=3) for i in range(len(selection))]
    else:
        [plt.errorbar(perccenters[i], Alist[i], yerr=Alist_error[i], \
        marker='.', ls='', color=colors[i], zorder=3) for i in range(len(selection))]
    
    if ('kids' in selection_name) or ('gama' in selection_name):
        # Plot mock amplitudes
        [plt.plot(perccenters_mock[i], Alist_mock[i], ls=linestyle, color=colors[i], zorder=1) \
        for i in range(len(mocksel))]

    plt.axhline(y=0., ls=':', color='black', zorder=2)
    plt.axvline(x=0.5, ls=':', color='black', zorder=2)

    ax2.axis([0.45,0.65,-0.002,0.002])
    ax2.tick_params(labelleft='off', top='off', left='off', right='off')

for ext in ['png', 'pdf']:

    plotfilename = '%s/trough_amplitude_%s'%(path_filename, selection_name)
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
if show:
    plt.show()
plt.close()


if ('kids' in selection_name) or ('gama' in selection_name):

    ## WEIGHT (Percentile)

    fig = plt.figure(figsize=(5,4))

    # Calculating the weight corresponding to Amplitude/error
    weightlist = Alist/Alist_error
    model_y = np.zeros([len(selection), len(model_x)])
    poly_param_weights = []

    # Calculating the weight corresponding to the mock Amplitude/error
    weightlist_mock = Alist_mock/Alist_error_mock
    model_y_mock = np.zeros([len(selection), len(model_x)])
    poly_param_weights_mock = []

    for i in range(len(selection)):
        
        # Fitting a polynomial curve to the data weight a.f.o. percentile
        poly_param_weights.append( np.polyfit(perccenters[i], weightlist[i], Npoly) )
        poly_func_weights = np.poly1d(poly_param_weights[i])
        model_y[i] = poly_func_weights(model_x)
        
        # Fitting a polynomial curve to the mock weight a.f.o. percentile
        poly_param_weights_mock.append( np.polyfit(perccenters[i], weightlist_mock[i], Npoly) )
        poly_func_weights_mock = np.poly1d(poly_param_weights_mock[i])
        model_y_mock[i] = poly_func_weights_mock(model_x)

        if 'slics' in mocksel[0]:
            plt.plot(perccenters[i], weightlist_mock[i],\
            label=labels[i], marker='o', ms='5', ls='', color=colors[i], zorder=3)

            #plt.plot(perccenters[i], weightlist[i], color=colors[i], zorder=1)
            
            plt.plot(model_x, model_y_mock[i], ls='--', color=colors[i], zorder=1)
        else:
            plt.plot(perccenters[i], weightlist[i],\
            label=labels[i], marker='o', ms='5', ls='', color=colors[i], zorder=3)

            plt.plot(model_x, model_y[i], color=colors[i], zorder=1)
 
    plt.axhline(y=0., ls=':', color='black', zorder=2)
    plt.axvline(x=0.5, ls=':', color='black', zorder=2)

    # Define the labels for the plot
    plt.ylabel(r'Shear signal to noise ($A/\sigma_{\rm A}$)')
    plt.xlabel(r"Density percentile rank $P(\theta_{\rm A})$")

    plt.legend(loc='best')

    for ext in ['png', 'pdf']:

        plotfilename = '%s/trough_weight_%s'%(path_filename, selection_name)
        plotname = '%s.%s'%(plotfilename, ext)

        plt.savefig(plotname, format=ext, bbox_inches='tight')
        
    print('Written plot:', plotname)
    if show:
        plt.show()
    plt.close()


## AMPLITUDE (delta)
fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot(111)

# Plot observed amplitudes

if ('kids' in selection_name) or ('gama' in selection_name):
    # Plot observed amplitudes
    if 'Z' in selection_name:
        [plt.errorbar(deltacenters[i], Alist[i], yerr=Alist_error[i], \
        label=labels[i], marker='o', ms=5, ls='', color=colors[i], ecolor='black', zorder=3) for i in range(len(selection))]
    else:
        [plt.errorbar(deltacenters[i], Alist[i], yerr=Alist_error[i], \
        label=labels[i], marker='.', ls='', color=colors[i], zorder=3) for i in range(len(selection))]

    # Plot mock amplitudes
    [plt.plot(deltacenters_mock[i], Alist_mock[i], ls=linestyle, color=mockcolors[i], alpha=valpha, zorder=1) \
    for i in range(len(mocksel))]
else:
    if 'miceZ' in selection_name:
        [plt.plot(deltacenters[i], Alist[i], label=labels[i],  color=colors[i], \
        marker='', ls='-', alpha=1., zorder=3) for i in range(len(selection))]
    else:
        [plt.plot(deltacenters[i], Alist[i], marker='', ls='-', alpha=valpha, zorder=3) \
                                        for i in range(len(selection))]

plt.axhline(y=0., ls=':', color='black', zorder=2)
plt.axvline(x=0., ls=':', color='black', zorder=2)

# Save plot

# Define the labels for the plot
if 'pc' in Runit:
    plt.ylabel(r'Comoving ESD Amplitude [h$_{%g}^{1/2}$ M$_{\odot}$/pc$^{3/2}$]'%(h*100))
    plt.axis([-1.,1.7,-3.,9.])
    if 'miceZ' in selection_name:
        plt.axis([-1.2,1.7,-3.,7.])
if 'arcmin' in Runit:
    plt.ylabel(r'Shear Amplitude [${\rm arcmin}^{1/2}]$')
    plt.axis([-0.8,1.5,-0.007,0.012])

plt.xlabel(r"Overdensity $\delta(\theta_{\rm A})$")

plt.legend(loc='upper left')

if 'arcmin' in Runit:
    # Include zoom-in of center
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left, bottom, width, height = [0.48, 0.16, 0.39, 0.25]
    ax2 = fig.add_axes([left, bottom, width, height])

    # Plot observed amplitudes
    if 'Z' in selection_name:
        [plt.errorbar(deltacenters[i], Alist[i], yerr=Alist_error[i], \
        label=labels[i], marker='o', ms=5, ls='', color=colors[i], ecolor='black', zorder=3) for i in range(len(selection))]
    else:
        [plt.errorbar(deltacenters[i], Alist[i], yerr=Alist_error[i], \
        label=labels[i], marker='.', ls='', color=colors[i], zorder=3) for i in range(len(selection))]

    if ('kids' in selection_name) or ('gama' in selection_name):
        # Plot mock amplitudes
        [plt.plot(deltacenters_mock[i], Alist_mock[i], ls=linestyle, color=colors[i], zorder=1) \
        for i in range(len(mocksel))]

    plt.axhline(y=0., ls=':', color='black', zorder=2)
    plt.axvline(x=0., ls=':', color='black', zorder=2)

    ax2.axis([-0.2,0.2,-0.002,0.002])
    ax2.tick_params(labelleft='off', top='off', left='off', right='off')
    plt.xticks(np.arange(-0.2,0.3,0.1))

for ext in ['png', 'pdf']:

    plotfilename = '%s/trough_amplitude_delta_%s'%(path_filename, selection_name)
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
if show:
    plt.show()
plt.close()


# Write weight catalog for KiDS/GAMA and MICE
if ('kids' in selection_name) or ('gama' in selection_name):
    
    selcat = np.array([selection, mocksel])
    poly_params = np.array([poly_param_weights, poly_param_weights_mock])
    
    for s in range(len(selcat)):

        # Import observed trough catalog
        path_troughcat = '/data2/brouwer/MergedCatalogues/trough_catalogs'
        
        # Select catalogue and data/mock fitted polynomial
        sel = selcat[s]
        poly_param_weights = poly_params[s]
       
        # Write weight catalog
        outputnames = []
        output = []

        for theta in range(len(thetalist)):
            
            troughcatname = 'trough_catalog_%s.fits'%sel[theta]
            
            # Full directory & name of the trough catalogue
            troughcatfile = '%s/%s'%(path_troughcat, troughcatname)
            troughcat = pyfits.open(troughcatfile, memmap=True)[1].data
            
            print('Creating weight catalog for:', troughcatfile)
            
#            Atheta = poly_func_amps(Ptheta)
#            poly_func_amps = np.poly1d(poly_param_amps[theta])
#            outputnames.append('Atheta%g'%thetalist[theta])
#            output.append(np.abs(Atheta))

            # Write weight fits-file
            Ptheta = troughcat['Ptheta%g'%thetalist[theta]]
            poly_func_weights = np.poly1d(poly_param_weights[theta])
            Wtheta = poly_func_weights(Ptheta)
            
            outputnames.append('Ptheta%g'%thetalist[theta])
            outputnames.append('Wtheta%g'%thetalist[theta])
            
            output.append(np.abs(Ptheta))
            output.append(np.abs(Wtheta))

        weightcatname = '%s/amplitude_trough_weights_%s.fits'%(path_troughcat, sel[0])
        utils.write_catalog(weightcatname, np.arange(len(Ptheta)), outputnames, output)

