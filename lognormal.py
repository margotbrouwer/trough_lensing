#!/usr/bin/python

import numpy as np
import sys
import trough_modules_all as utils

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

import scipy.optimize as optimization

A = 0.05
sigmalist = np.array([1.0, 0.5, 0.25])*100.
mulist = np.array([0., 0.5, 1.0])*100.

"""
def calc_model(x, A, mu, sigma):
    
    model_y = -A * 1/(x * sigma * np.sqrt(2*np.pi)) * \
    np.exp( -(np.log(x)-mu)**2./(2*sigma**2.) )
    
    return model_y
"""

def calc_model(x, A, eps, alpha, kappa, sigma):
    
    y = -(1/kappa) * np.log(1 - kappa*(x-eps)/alpha)
    
    normal = 1/np.sqrt(2.*np.pi*sigma**2.) * \
            np.exp( -(y-eps)**2. / (2.*sigma**2.) )
    
    model_y = A * normal / (alpha - kappa*(x-eps))
    
    return model_y

# Defining the paths to the data
blind = 'A'

perclist = np.arange(0.05, 0.4, 0.05)
percnames = ['0p05', '0p1', '0p15', '0p2', '0p25', '0p3', '0p35']

path_sheardata = 'data2/brouwer/shearprofile/trough_results_Feb'
path_lenssel = ['No_bins_gama_absmag/Ptheta5_0_%s'%pn for pn in percnames]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins25_1_300_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datalabels = [r'Troughs, $\theta=5$ arcmin, $M_r<-19.7$']
plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/optimize_perc'


esdfiles = np.array([('/%s/%s/%s/%s'%(path_sheardata, path_lenssel[i], path_cosmo, path_filename)) \
           for i in range(len(path_lenssel))])

covfiles = np.array([e.replace('bins_%s.txt'%blind, 'matrix_%s.txt'%blind) for e in esdfiles])

# Importing the shearprofiles and lens IDs
data_x, data_y, error_h, error_l = utils.read_esdfiles(esdfiles)

A, eps, alpha, kappa, sigma = [-0.01, 2., 1., -0.5, 1.]
model_y = calc_model(data_x[0], A, eps, alpha, kappa, sigma)
plt.plot(data_x[0], model_y)

plt.errorbar(data_x[3], data_y[3], yerr=error_h[3], ls='', marker='.')
plt.xscale('log')
plt.show()


for i in range(len(percnames)):

    A, eps, alpha, kappa, sigma = optimization.curve_fit(calc_model, data_x[i], data_y[i], [-0.01, 2., 1., -1., 1.], error_h[i])[0]
    print(A, eps, alpha, kappa, sigma)

    model_y = calc_model(data_x[i], A, eps, alpha, kappa, sigma)

    plt.plot(data_x[i], model_y)
    plt.errorbar(data_x[i], data_y[i], yerr=error_h[i], ls='', marker='.')

    plt.xscale('log')
    
    plt.show()

