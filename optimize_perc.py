
# coding: utf-8

# In[17]:

#!/usr/bin/python

# Import the necessary libraries
import sys

import numpy as np
import pyfits
import os
import trough_modules_all as utils

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

from matplotlib import gridspec
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# Make use of TeX
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})



# Colours
# Blue, green, turquoise, cyan
blues = ['#332288', '#44AA99', '#117733', '#88CCEE']

# Light red, Red, light pink, pink
reds = ['#CC6677', '#882255', '#CC99BB', '#AA4499']
colors = np.array([reds,blues])


# Defining the paths to the data
blind = 'A'

Nbins = 1
Nrows = 1

perclist = np.arange(0.0, 0.4, 0.05)
percnames = ['0p05', '0p1', '0p15', '0p2', '0p25', '0p3', '0p35']
Npercs = len(percnames)

path_sheardata = 'data2/brouwer/shearprofile/trough_results_Feb'
path_lenssel = ['No_bins_gama_absmag/Ptheta5_0_%s'%percnames[perc] for perc in np.arange(Npercs)]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins25_1_300_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datalabels = [r'Troughs, $\theta=5$ arcmin, $M_r<-19.7$']
plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/optimize_perc'


esdfiles = np.array([('/%s/%s/%s/%s'%(path_sheardata, path_lenssel[i], path_cosmo, path_filename)) \
           for i in range(len(path_lenssel))])

covfiles = np.array([e.replace('bins_%s.txt'%blind, 'matrix_%s.txt'%blind) for e in esdfiles])

# Importing the shearprofiles and lens IDs
data_x, data_y, error_h, error_l = utils.read_esdfiles(esdfiles)

model = np.zeros(len(data_x[0]))

chi2list = np.zeros(Npercs)
chi2covlist = np.zeros(Npercs)

for perc in range(Npercs):

    data = data_y[perc]
    error = error_h[perc]
    covariance = (np.loadtxt(covfiles[perc]).T)

    chi2list[perc] = np.sum(data**2/error**2)
    chi2covlist[perc] = utils.calc_chi2(data, model, covariance, Nbins)

    chilist, chicovlist = np.sqrt(chi2list), np.sqrt(chi2covlist)


print('Chi (without covariance):', chilist)
print('Chi (with covariance):', chicovlist)
print('Ratio:', chilist/chicovlist)

fig = plt.figure(figsize=(8,6))
plt.plot(perclist[1::], chicovlist, 'o', ls='--', label='With covariance')
plt.plot(perclist[1::], chilist, 'o', ls='--', label='Without covariance')

plt.xlabel(r'Percentage of troughs $P(\theta)$')
plt.ylabel(r'Detection ($\chi$ value)')
plt.legend(loc='best')

# Save plot
plotname = '%s.png'%plotfilename

plt.savefig(plotname, format='png', bbox_inches='tight')
print('Written plot:', plotname)

plt.show()
plt.clf

# Import trough catalog

path_troughcat = '/data2/brouwer/MergedCatalogues/'
troughcatname = '/trough_catalog_gama_absmag_masked.fits'

# Full directory & name of the trough catalogue
troughcatfile = '%s/%s'%(path_troughcat, troughcatname)
troughcat = pyfits.open(troughcatfile, memmap=True)[1].data

# List of the observables of all sources in the KiDS catalogue
Ptheta = troughcat['Ptheta5']

Wtheta = np.zeros(len(Ptheta))
for perc in range(Npercs):
    
    weight = chi2covlist[perc]/np.amax(chi2covlist)
    
    percmask = (perclist[perc] < Ptheta) & (Ptheta < perclist[perc+1])
    Wtheta[percmask] = weight

    print('Percentile: %g - %g, weight = %g'%(perclist[perc], perclist[perc+1], weight))    

# Write weight fits-file
outputnames = ['Ptheta5', 'Wtheta5']
output = [Ptheta, Wtheta]
utils.write_catalog('/data2/brouwer/MergedCatalogues/trough_weights.fits', np.arange(len(Wtheta)), outputnames, output)
