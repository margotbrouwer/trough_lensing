
# coding: utf-8

# In[17]:

#!/usr/bin/python

# Import the necessary libraries
import sys

import numpy as np
import pyfits
import os

import scipy.optimize as optimization
import trough_modules_all as utils

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
from matplotlib import rc, rcParams

from scipy.stats import chi2
from matplotlib import gridspec
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# Make use of TeX
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})

# Model to fit the troughs
def trough_model(x, A):
    model_y = A/x
    return model_y



# Colours
# Blue, green, turquoise, cyan
blues = ['#332288', '#44AA99', '#117733', '#88CCEE']

# Light red, Red, light pink, pink
reds = ['#CC6677', '#882255', '#CC99BB', '#AA4499']
colors = np.array([reds,blues])


# Defining the paths to the data
blind = 'A'
thetalist = np.array([5., 10., 15., 20.]) # in arcmin
#thetalist = np.array([3.163, 6.326, 9.490, 12.653]) # in arcmin


# Import shear and random profiles


# Weighted troughs: sizes

Runit = 'arcmin'
path_sheardata = 'data2/brouwer/shearprofile/trough_results_Apr'

path_lenssel = np.array([ ['No_bins_gama_absmag_complex/Pmasktheta%g_0p8_1-delta%g_%s_lw-Wtheta%g'%(theta,theta,delta,theta) \
                for theta in thetalist ] for delta in ['minf_0', '0_inf'] ])

path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datatitles = [r'$\theta_{\rm A} = %g$ arcmin'%theta for theta in thetalist]
datalabels = ['Troughs', 'Ridges']
plotfilename = '/data2/brouwer/shearprofile/trough_results_Apr/Plots/troughs_gama_weighted'
Nrows = 2

path_randoms = np.array([ ['No_bins_gama_randoms/Pmasktheta%g_0p8_1'%theta
                for theta in thetalist ] for delta in ['minf_0', '0_inf'] ])


"""

# Weighted troughs: Redshifts

Runit = 'Mpc'
h=0.7

thetalist = np.array(['10', '6p326'])
samplelist = np.array(['lowZ', 'highZ'])

path_sheardata = 'data2/brouwer/shearprofile/trough_results_Apr'

path_lenssel = np.array([ ['No_bins_gama_%s_complex/Pmasktheta%s_0p8_1-delta%s_%s_lw-Wtheta%s'%(samplelist[i], thetalist[i], thetalist[i], delta, thetalist[i]) \
                for i in range(len(thetalist)) ] for delta in ['minf_0', '0_inf'] ])

path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins10_0p5_20_Mpc/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datalabels = [r'$0.1<z<0.197$', r'$0.197<z<0.3$']
plotfilename = '/data2/brouwer/shearprofile/trough_results_Apr/Plots/troughs_gama_redshifts_weighted'

path_randoms = np.array([ ['No_bins_gama_randoms/Pmasktheta%s_0p8_1'%theta
                for theta in thetalist ] for delta in ['minf_0', '0_inf'] ])



# Randoms

Runit = 'arcmin'

path_sheardata = 'data2/brouwer/shearprofile/trough_results_Apr'
path_lenssel = [['No_bins_gama_randoms/Pmasktheta5_0p%s_1'%masknum for masknum in range(6,10)]]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datalabels = [r'GAMA ($A_{\rm eff} > %g$)'%n for n in np.arange(0.6, 1.0, 0.1)]

plotfilename = '/data2/brouwer/shearprofile/trough_results_Apr/Plots/troughs_gama_randoms'



# KiDS vs GAMA

Runit = 'arcmin'
path_sheardata = 'data2/brouwer/shearprofile/trough_results_Apr'

# Troughs
#path_lenssel = [['No_bins_%s/Pmasktheta5_0p%i_inf-Ptheta5_0_0p2'%(sample, n) for n in np.arange(6, 10)] \
#                for sample in ['gama_absmag_complex', 'gama_kids_absmag_complex', 'kids_absmag_nomask'] ]
#plotfilename = '/data2/brouwer/shearprofile/trough_results_Apr/Plots/troughs_kids_gama'
#Nrows = 2

# Randoms
path_lenssel = [['No_bins_%s/Pmasktheta5_0p%i_inf'%(sample, n) for n in np.arange(8, 9)] \
                for sample in ['gama_absmag_randoms', 'gama_kids_randoms', 'kids_nomask_randoms'] ]
plotfilename = '/data2/brouwer/shearprofile/trough_results_Apr/Plots/troughs_kids_gama_randoms'
Nrows = 1

path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance/'
path_filename = 'No_bins_%s.txt'%(blind)

datatitles = [r'($A_{\rm eff} > %g$)'%n for n in np.arange(0.8, 0.9, 0.1)]
datalabels = ['GAMA', 'KiDS (GAMA area)', 'KiDS (full, no mask)']


#path_randoms = ['No_bins_kids_randoms/Pmasktheta5_0p6_1'%n for n in np.arange(5, 9)]
#path_randoms = ['No_bins_gama_randoms/Pmasktheta5_0p8_1']


"""

print(path_lenssel)
print(datalabels)

Nbins = np.shape(path_lenssel)
Nsize = np.size(path_lenssel)
path_lenssel = np.reshape(path_lenssel, [Nsize])

try:
    path_randoms = np.reshape(path_randoms, [Nsize])
except:
    pass
    
Nlabels = np.size(datalabels)
datalabels = np.reshape(datalabels, [Nlabels])

print('Profiles, Bins:', Nbins)


esdfiles = np.array([('/%s/%s/%s/%s'%(path_sheardata, path_lenssel[i], path_cosmo, path_filename)) \
           for i in range(len(path_lenssel))])

lensIDfiles = np.array([e.replace('_%s.txt'%blind, '_lensIDs.txt') for e in esdfiles])
covfiles = np.array([e.replace('bins_%s.txt'%blind, 'matrix_%s.txt'%blind) for e in esdfiles])

# Importing the shearprofiles and lens IDs
data_x, data_y, error_h, error_l = utils.read_esdfiles(esdfiles)
lensIDs = np.array([np.loadtxt(x) for x in lensIDfiles])


print('Import random signal:')

random_esdfile = np.array(['/%s/%s/%s/%s'%(path_sheardata, path_random, path_cosmo, path_filename) for path_random in path_randoms])
random_data_x, random_data_y, random_error_h, random_error_l = utils.read_esdfiles(random_esdfile)


# Subtract random signal
data_y = data_y-random_data_y
error_h = np.sqrt(error_h**2. + random_error_h**2)
error_l = np.sqrt(error_l**2. + random_error_l**2)

"""
# Plot random signal
for N in range(len(path_randoms)):
    plt.errorbar((1.+N/10.)*data_x[N], random_data_y[N], yerr=[random_error_l[N], random_error_h[N]], \
    ls='', marker='.', label=datalabels[N])
    plt.axhline(y=0., ls=':', color='black')


plt.autoscale(enable=False, axis='both', tight=None)
plt.axis([2,100,-1.5e-3,1.5e-3])
plt.ylim(-1.5e-3,1.5e-3)
plt.xscale('log')

plt.legend(loc='best')
plt.xlabel(xlabel)
plt.ylabel(ylabel)

for ext in ['png', 'pdf']:
    plotname = '%s_randoms.%s'%(plotfilename, ext)
    plt.savefig(plotname, format=ext, bbox_inches='tight')

plt.show()
plt.close()
"""

# Create the plot

Ncolumns = int(Nbins[1]/Nrows)

# Plotting the ueber matrix
fig = plt.figure(figsize=(Ncolumns*5,Nrows*5))
canvas = FigureCanvas(fig)

gs_full = gridspec.GridSpec(1,1)
gs = gridspec.GridSpecFromSubplotSpec(Nrows, Ncolumns, wspace=0, hspace=0, subplot_spec=gs_full[0,0])

ax = fig.add_subplot(gs_full[0,0])

for N1 in range(Nrows):
    for N2 in range(Ncolumns):

        ax_sub = fig.add_subplot(gs[N1, N2])
        
        N = np.int(N1*Ncolumns + N2)
        
        for Nplot in range(Nbins[0]):
            
            Ndata = N + Nplot*(Nbins[1])
            data_x_plot = data_x[Ndata]# * (1.+0.05*Nplot)

            if Nsize==Nbins:
                ax_sub.errorbar(data_x_plot, data_y[Ndata], yerr=[error_l[Ndata], error_h[Ndata]], \
                ls='', marker='.')
            else:
                ax_sub.errorbar(data_x_plot, data_y[Ndata], yerr=[error_l[Ndata], error_h[Ndata]], \
                ls='', marker='.', label=datalabels[Nplot])

            
        #ax_sub.errorbar(data_x[N], -data_y[N], yerr=[error_l[N], error_h[N]], \
        #ls='', marker='.', alpha=0.25, color='blue')
            
#        for i in range(len(path_lenssel)):
            
        #(error1_l[N])[(error1_l[N])>=data1_y[N]] = ((data1_y[N][(error1_l[N])>=data1_y[N]])*0.9999999999)
        #(error2_l[N])[(error2_l[N])>=data2_y[N]] = ((data2_y[N][(error2_l[N])>=data2_y[N]])*0.9999999999)
        
        # Plot the data and title
    
        ax_sub.axvline(x=thetalist[N], color='blue', ls='--')
        ax_sub.axvline(x=1.2*thetalist[N], color='green', ls='--')
        
        ax_sub.axhline(y=0., ls=':', color='black')

        # Plot the models
        
        ax_sub.xaxis.set_label_position('top')
        ax_sub.yaxis.set_label_position('right')

        ax.tick_params(labelleft='off', labelbottom='off', top='off', bottom='off', left='off', right='off')

        if (N1+1) != Nrows:
            ax_sub.tick_params(axis='x', labelbottom='off')
        if N2 != 0:
            ax_sub.tick_params(axis='y', labelleft='off')
        
        plt.autoscale(enable=False, axis='both', tight=None)

        # Define the labels for the plot
        if Runit == 'Mpc':
            plt.axis([0.5,20,-3,9])
            plt.ylim(-3,9)
            
            xlabel = r'Radial distance $R$ (%s/h$_{%g}$)'%(Runit, h*100)
            ylabel = r'ESD $\langle\Delta\Sigma\rangle$ [h$_{%g}$ M$_{\odot}$/pc$^2$]'%(h*100)
            
            ax.xaxis.set_label_coords(0.5, -0.15)
            ax.yaxis.set_label_coords(-0.05, 0.5)
        else:
            plt.axis([2,100,-1.5e-3,3e-3])
            plt.ylim(-1.5e-3,1.99e-3)

            xlabel = r'Opening angle $\theta$ (arcmin)'
            ylabel = r'Shear $\gamma$'
            
            ax.xaxis.set_label_coords(0.5, -0.07)
            ax.yaxis.set_label_coords(-0.2, 0.5)

        plt.xscale('log')
        #plt.yscale('log')
        
        if Nbins[1]>1:
            plt.title(datatitles[N], x = 0.17, y = 0.9, fontsize=16)
        
        #lgd = ax_sub.legend(bbox_to_anchor=(2.1, 1.4)) # side
        #lgd = ax_sub.legend(bbox_to_anchor=(0.6, 2.7)) # top
        plt.legend(loc='upper right') # handlelength=0, handletextpad=0, numpoints=1)

# Define the labels for the plot
ax.set_xlabel(xlabel, fontsize=14)
ax.set_ylabel(ylabel, fontsize=14)



#ax.xaxis.label.set_size(17)
#ax.yaxis.label.set_size(17)

"""
else:
    
    for N in range(Nbins):
        
        plt.errorbar((1.+N/10.)*data_x[N], data_y[N], yerr=[error_l[N], error_h[N]], \
            ls='', marker='.', label=datalabels[N])
        
    plt.axhline(y=0., ls=':', color='black')

    plt.autoscale(enable=False, axis='both', tight=None)
    plt.axis([2,100,-1.5e-3,1.5e-3])
    plt.ylim(-1.5e-3,1.5e-3)

    plt.xscale('log')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='best')
"""

plt.tight_layout()

# Save plot
for ext in ['png', 'pdf']:
    plotname = '%s.%s'%(plotfilename, ext)
    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written: ESD profile plot:', plotname)

plt.show()
plt.clf


# Calculate detection significance

model = np.zeros(len(data_x[0]))

chi2list = np.zeros(Nsize)
chi2covlist = np.zeros(Nsize)
Alist = []
Alist_error = []

for N in range(Nsize):
    
    datax = data_x[N]
    datay = data_y[N]
    
    error = error_h[N]
    covariance = np.array(np.loadtxt(covfiles[N])).T

    chi2list[N] = np.sum(datay**2/error**2)
    chi2covlist[N] = utils.calc_chi2(datay, model, covariance, 1)

    # Signal to noise of the first bin beyond the trough size
    thetaNlist = np.append(thetalist, thetalist)
    
    if Runit == 'arcmin':
        xmin = thetaNlist[N]*1.2
        xmax = 70.
    if Runit == 'Mpc':
        xmin = Rlist[N]*1.2
        xmax = 10.

    xmask = (xmin < datax) & (datax < xmax)
    xwhere = np.array(np.where(xmask))[0]

    # Signal to noise of the Amplitude fit
    
    # With covariance
    ind = np.lexsort((covariance[3,:], covariance[1,:], covariance[2,:], covariance[0,:]))
    covmatrix = np.reshape(covariance[4][ind], [len(datax), len(datax)])
    covmatrix = covmatrix[int(xwhere[0]):int(xwhere[-1]+1), int(xwhere[0]):int(xwhere[-1]+1)]

    A, Acov = optimization.curve_fit(f=trough_model, xdata=datax[xmask], ydata=datay[xmask], p0=[0.], \
    sigma=covmatrix, absolute_sigma=True)
    
    Alist = np.append(Alist, A)
    Alist_error = np.append(Alist_error, np.sqrt(Acov[0,0]))
    
    print('First bin S/N):', datay[xwhere[0]]/error[xwhere[0]])
    print('Amplitude fit S/N', A/np.sqrt(Acov[0,0]))
    print

# Calculate the detection significance
chilist, chicovlist = np.sqrt(chi2list), np.sqrt(chi2covlist)

problist = [chi2.sf(chi2covlist[i], len(data_x[0])) for i in range(Nsize)] # Survival Function (SF = 1 - CDF)
sigmalist = [chi2.ppf(1.-problist[i], len(data_x[0])) for i in range(Nsize)]

print('Chi2 (without covariance):', chi2list)
print('Chi2 (with covariance):', chi2covlist)
print('Probability:', problist)
print('Sigma:', sigmalist)
