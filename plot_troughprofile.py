
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

# In[18]:

# Defining the paths to the data
blind = 'A'
Nbins = 1
Nrows = 1

path_sheardata = 'data2/brouwer/shearprofile/trough_results_Feb'

"""
percnames = ['0','0p05','0p1','0p15','0p2','0p25','0p3','0p35','0p4','0p45','0p5']
perclist = np.arange(0., 0.55, 0.05)

path_lenssel = ['No_bins_gama_absmag/Pmasktheta5_0p8_1-Ptheta5_%s_%s'%(percnames[i], percnames[i+1]) for i in range(Nbins)]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins25_1_300_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datalabels = [r'$%g<P<%g$'%(perclist[i], perclist[i+1]) for i in range(Nbins)]
"""

path_lenssel = ['No_bins_gama_absmag/Pmasktheta5_0p8_1-Ptheta5_0_0p5_lw-Wtheta5/']
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins25_1_300_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datalabels = [r'Troughs: $\theta = 5$ arcmin, weighted by $\chi^2$']
#"""

plotfilename = '/data2/brouwer/shearprofile/trough_results_Feb/Plots/chi2_weighted_troughs'


esdfiles = np.array([('/%s/%s/%s/%s'%(path_sheardata, path_lenssel[i], path_cosmo, path_filename)) \
           for i in range(len(path_lenssel))])

lensIDfiles = np.array([e.replace('_%s.txt'%blind, '_lensIDs.txt').replace('randomsub_','') for e in esdfiles])

# Importing the shearprofiles and lens IDs
data_x, data_y, error_h, error_l = utils.read_esdfiles(esdfiles)
lensIDs = np.array([np.loadtxt(x) for x in lensIDfiles])


# Create a plot of the halo model fit

# Define the labels for the plot
xlabel = r'Angle $\theta$ [arcmin]'
ylabel = r'Shear $\gamma$'

Ncolumns = int(Nbins/Nrows)

# Plotting the ueber matrix
fig = plt.figure(figsize=(Ncolumns*6,6))
canvas = FigureCanvas(fig)

gs_full = gridspec.GridSpec(1,1)
gs = gridspec.GridSpecFromSubplotSpec(Nrows, Ncolumns, wspace=0, hspace=0, subplot_spec=gs_full[0,0])

ax = fig.add_subplot(gs_full[0,0])

for N1 in range(Nrows):
    for N2 in range(Ncolumns):

        N = np.int(N1*Ncolumns + N2)
        ax_sub = fig.add_subplot(gs[N1, N2])
            
#        for i in range(len(path_lenssel)):
            
        #(error1_l[N])[(error1_l[N])>=data1_y[N]] = ((data1_y[N][(error1_l[N])>=data1_y[N]])*0.9999999999)
        #(error2_l[N])[(error2_l[N])>=data2_y[N]] = ((data2_y[N][(error2_l[N])>=data2_y[N]])*0.9999999999)
        
        # Plot the data and title
        title = r'Bin %i'%(N+1)
        
        #data_x_plot = data_x[N]*(1+0.1*N)
        ax_sub.errorbar(data_x[N], data_y[N], yerr=[error_l[N], error_h[N]], \
        ls='', marker='.', label=datalabels[N])
    
        ax_sub.axvline(x=[5.])
        
        ax_sub.axhline(y=0., ls=':')

        # Plot the models
        
        
        ax_sub.xaxis.set_label_position('top')
        ax_sub.yaxis.set_label_position('right')

        ax.tick_params(labelleft='off', labelbottom='off', top='off', bottom='off', left='off', right='off')
                
        if N2 != 0:
            ax_sub.tick_params(axis='y', labelleft='off')
        
        plt.autoscale(enable=False, axis='both', tight=None)
        plt.axis([1e0,3.5e2,-1.5e-3,1.5e-3])
        plt.ylim(-1.5e-3,1.5e-3)
#        plt.ylim(0,1e2)

        plt.xscale('log')
        #plt.yscale('log')
        
        #plt.title(title, x = 0.6, y = 0.8)
        
        #lgd = ax_sub.legend(bbox_to_anchor=(2.1, 1.4)) # side
        #lgd = ax_sub.legend(bbox_to_anchor=(0.6, 2.7)) # top
        plt.legend()

# Define the labels for the plot
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

ax.xaxis.set_label_coords(0.5, -0.07)
ax.yaxis.set_label_coords(-0.07, 0.5)

#ax.xaxis.label.set_size(17)
#ax.yaxis.label.set_size(17)


plt.tight_layout()

# Save plot
for ext in ['png', 'pdf']:
    plotname = '%s.%s'%(plotfilename, ext)
    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written: ESD profile plot:', plotname)

plt.show()
plt.clf

