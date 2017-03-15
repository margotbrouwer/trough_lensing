
# coding: utf-8

# In[17]:

#!/usr/bin/python

# Import the necessary libraries
import sys

import numpy as np
import pyfits
import os
import trough_modules_all as utils
import scipy.optimize as optimization

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

Nbins = 1
dperc = 0.05

# Troughs
percmin = 0.
percmax = 0.5
percnames = ['0','0p05','0p1','0p15','0p2','0p25','0p3','0p35','0p4','0p45','0p5']

# Ridges
#percmin = 0.5
#percmax = 1.0
#percnames = ['0p5','0p55','0p6','0p65','0p7','0p75','0p8','0p85','0p9','0p95','0p1']

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

datalabels = [r'Troughs, $\theta=5$ arcmin, $M_r<-19.7$']
plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/perc_weights'


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

fig = plt.figure(figsize=(6,6))
plt.plot(perccenters, chi2covlist, 'o', ls='', label=r'$\chi^2(P)$ per percentile bin')
#plt.plot(perccenters, chi2list, 'o', ls='--', label='$\chi^2(P)$ without covariance')


# Fitting the optimal curve to the chi2 a.f.o. percentile
def calc_model(x, A, B, C):
    
    model_y = A*x**2 + B*x + C
    
    return model_y

A, B, C = optimization.curve_fit(calc_model, perccenters, chi2covlist, [-300., 20., 120.])[0]
model_y = calc_model(perccenters, A, B, C)

# Plot best-fit polynomial
plt.plot(perccenters, model_y, label=r'Best-fit $2^{\rm nd}$ degree polynomial')


# Write weight fits-file

# Import trough catalog
path_troughcat = '/data2/brouwer/MergedCatalogues/'
troughcatname = '/trough_catalog_gama_absmag_masked.fits'

# Full directory & name of the trough catalogue
troughcatfile = '%s/%s'%(path_troughcat, troughcatname)
troughcat = pyfits.open(troughcatfile, memmap=True)[1].data

# List of the observables of all sources in the KiDS catalogue
Ptheta = troughcat['Ptheta5']
Wtheta = calc_model(Ptheta, A, B, C)

outputnames = ['Ptheta5', 'Wtheta5']
output = [Ptheta, Wtheta]
utils.write_catalog('/data2/brouwer/MergedCatalogues/trough_weights.fits', np.arange(len(Wtheta)), outputnames, output)


# Calculating covariance of the weighted trough signal
path_lenssel_weighted = 'No_bins_gama_absmag/Pmasktheta5_0p8_1-Ptheta5_0_0p5_lw-Wtheta5/'
esdfiles_weighted = np.array(['/%s/%s/%s/%s'%(path_sheardata, path_lenssel_weighted, path_cosmo, path_filename)])
data_x, data_y_weighted, error_h_weighted, error_l_weighted = utils.read_esdfiles(esdfiles_weighted)
covfiles_weighted = np.array([e.replace('bins_%s.txt'%blind, 'matrix_%s.txt'%blind) for e in esdfiles_weighted])

data_weighted = data_y_weighted[0]
error_weighted = error_h_weighted[0]
covariance_weighted = (np.loadtxt(covfiles_weighted[0]).T)

chi2_weighted = np.sum(data_weighted**2/error_weighted**2)
chi2cov_weighted = utils.calc_chi2(data_weighted, model, covariance_weighted, Nbins)

plt.axhline(y=chi2cov_weighted, ls='--', label=r'$\chi^2$ of the weighted troughs')


# Save plot

# Define the labels for the plot
plt.xlabel(r'Trough percentile $P(\theta)$')
plt.ylabel(r'Lensing detection $\chi^2$')

plt.ylim(50,400)

for ext in ['png', 'pdf']:
    
    plt.legend(loc='best')
    plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/chi2_trough_weights'
    plotname = '%s.%s'%(plotfilename, ext)

    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written: Eplot:', plotname)
plt.show()



