
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



centering = 'BCG'
purpose = 'shearcatalog'
Ncat = 1
Runit = 'kpc'
lens_weights = {'None': ''}

O_matter = 0.315
O_lambda = 1.-O_matter
Ok = 0.
h = 0.7

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

path_lenssel = ['No_bins_gama_all/Ptheta5_0_0p2', 'No_bins_gama_old/Ptheta5_0_0p2']
#path_lenssel = ['No_bins/Ptheta5_0_0p2-Stheta5_1', 'No_bins/Ptheta10_0_0p2-Stheta10_1']

path_cosmo = ['ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins25_1_300_arcmin/shearcovariance',\
              'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins25_1_300_arcmin/shearcovariance']
path_filename = 'No_bins_%s.txt'%(blind)

datalabels = [r'$\theta=5$, All', r'$\theta=5$, Selected']
#datalabels = [r'$\theta=5$, Ellipticals', r'$\theta=10$, Ellipticals']

esdfiles = np.array([('/%s/%s/%s/%s'%(path_sheardata, path_lenssel[i], path_cosmo[i], path_filename)) for i in range(len(path_lenssel))])

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
fig = plt.figure(figsize=(5,3.5))
canvas = FigureCanvas(fig)

gs_full = gridspec.GridSpec(1,1)
gs = gridspec.GridSpecFromSubplotSpec(Nrows, Ncolumns, wspace=0, hspace=0, subplot_spec=gs_full[0,0])

ax = fig.add_subplot(gs_full[0,0])

for N1 in range(Nrows):
    for N2 in range(Ncolumns):

        N = np.int(N1*Ncolumns + N2)
        ax_sub = fig.add_subplot(gs[N1, N2])
            
        for i in range(len(path_lenssel)):
            
            #(error1_l[N])[(error1_l[N])>=data1_y[N]] = ((data1_y[N][(error1_l[N])>=data1_y[N]])*0.9999999999)
            #(error2_l[N])[(error2_l[N])>=data2_y[N]] = ((data2_y[N][(error2_l[N])>=data2_y[N]])*0.9999999999)
            
            # Plot the data and title
            title = r'Bin %i'%(N+1)
            ax_sub.errorbar(data_x[i]*(1+0.1*i), data_y[i], yerr=[error_l[i], error_h[i]], ls='', marker='.', color=colors[i,N1], label=datalabels[i])
        
            ax_sub.axvline(x=[5.,10][i], color=colors[i,N1])

        # Plot the models
        
        
        ax_sub.xaxis.set_label_position('top')
        ax_sub.yaxis.set_label_position('right')

        ax.tick_params(labelleft='off', labelbottom='off',                        top='off', bottom='off', left='off', right='off')
                
        if N2 != 0:
            ax_sub.tick_params(axis='y', labelleft='off')
        
        plt.autoscale(enable=False, axis='both', tight=None)
        plt.axis([1e0,3e2,-3e-3,3e-3])
        plt.ylim(-3e-3,3e-3)
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
plotfilename = 'trough_test_all'
plotname = '%s.png'%plotfilename

plt.savefig(plotname, format='png', bbox_inches='tight')
print('Written: ESD profile plot:', plotname)

plt.show()
plt.clf

