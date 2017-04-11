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

samples = ['gama_lowZ_complex_10Mpc', 'gama_highZ_complex_6.326Mpc']
labels = [r'$0.1<Z<0.197$', r'$0.197<Z<0.3$']

# Import amplitudes
path_filename = 'data2/brouwer/shearprofile/trough_results_Apr/Plots/'
filenames = ['/%s/trough_amplitudes_%s.txt'%(path_filename, sample) for sample in samples]

amplitude_data = np.array([np.loadtxt(filename).T for filename in filenames])

perccenters = amplitude_data[:,0]
Alist = amplitude_data[:,1]
Alist_error = amplitude_data[:,2]

## AMPLITUDE (Percentile)
fig = plt.figure(figsize=(5,4))

[plt.errorbar(perccenters[i], Alist[i], yerr=Alist_error[i], label=labels[i], marker='.', ls=':') for i in range(len(filenames))]
plt.legend(loc='upper left')

plt.axhline(y=0., ls=':', color='black')
plt.axvline(x=0.5, ls=':', color='black')

# Define the labels for the plot
plt.ylabel(r'Lensing detection Amplitude')
plt.xlabel(r"Percentile $P(\theta)$")

for ext in ['png', 'pdf']:

    plotfilename = '/%s/trough_amplitude_redshifts'%path_filename
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written plot:', plotname)
plt.show()
