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

from matplotlib import gridspec
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

import matplotlib.style
matplotlib.style.use('classic')

# Make use of TeX
rc('text',usetex=True)

# Change all fonts to 'Computer Modern'
rc('font',**{'family':'serif','serif':['Computer Modern']})

# Model to fit the troughs
def trough_model(x, A):
    model_y = A*x**-0.5
    return model_y

# Model to fit the mocks
def mock_model(x, A, B):
    model_y = A*x**B
    return model_y


# Colours
# Blue, green, turquoise, cyan
blues = ['#332288', '#44AA99', '#117733', '#88CCEE']

# Light red, Red, light pink, pink
reds = ['#CC6677', '#882255', '#CC99BB', '#AA4499']
#colors = np.array([reds,blues])

#colors = ['#0571b0', '#92c5de', '#d7191c', '#fdae61']
colors = ['blue', 'green', 'red', 'orange']


# Defining the paths to the data
blind = 'A'
thetalist = np.array([5., 10., 15., 20.]) # in arcmin
#thetalist = np.array([5.]) # in arcmin
#thetalist = np.array([3.163, 6.326, 9.490, 12.653]) # in arcmin


# Import shear and random profiles

"""


# Fiducial troughs:sizes

Runit = 'arcmin'
path_sheardata = 'data2/brouwer/shearprofile/trough_results'

path_lenssel = np.array([ ['No_bins_kids_mice_complex/Pmasktheta%g_0p8_inf-Ptheta%g_%s'%(theta,theta,perc) \
                for theta in thetalist ] for perc in ['0_0p2', '0p8_1'] ])

path_cosmo = 'ZB_0p1_0p9-Om_0p25-Ol_0p75-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datatitles = [r'$\theta_{\rm A} = %g$ arcmin'%theta for theta in thetalist]
datalabels = [r'Troughs $(P(\theta_{\rm A}) < 0.2)$', r'Ridges $(P(\theta_{\rm A}) > 0.8)$']
plotfilename = '/data2/brouwer/shearprofile/trough_results_May/Plots/troughs_kids_fiducial'
Nrows = 2
plotfit = True

path_randoms = np.array([ ['No_bins_kids_randoms_complex/Pmasktheta%g_0p8_inf'%theta
                for theta in thetalist ] for perc in ['0_0p2', '0p8_1'] ])

thetalist = np.array([5., 10., 15., 20.]*2)

"""

# Weighted troughs: sizes

Runit = 'arcmin'
path_sheardata = 'data2/brouwer/shearprofile/trough_results_final'

path_lenssel = np.array([ ['No_bins_kids_mice_complex/Pmasktheta%g_0p8_inf-delta%g_%s_lw-Wtheta%g'%(thetalist[i], thetalist[i], delta, thetalist[i]) \
                for i in range(len(thetalist)) ] for delta in ['minf_0', '0_inf'] ])

path_cosmo = 'ZB_0p1_0p9-Om_0p25-Ol_0p75-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datatitles = [r'$\theta_{\rm A} = %g$ arcmin'%theta for theta in thetalist]
datalabels = [r'Weighted trough profile', r'Weighted ridge profile']
plotfilename = '/%s/Plots/troughs_kids_weighted'%path_sheardata
Nrows = 2
plotfit = True

path_randoms = np.array([ ['No_bins_kids_randoms_complex/Pmasktheta%g_0p8_inf'%theta
                for theta in thetalist ] for delta in ['minf_0', '0_inf'] ])

path_mockdata = 'data2/brouwer/shearprofile/trough_results_final'
path_mocksel = np.array([ ['No_bins_mice_all_nomask-1/Pmasktheta%g_0p8_inf-delta%g_%s.txt'%(thetalist[i], thetalist[i], delta) \
                for i in range(len(thetalist)) ] for delta in ['minf_0', '0_inf'] ])
mocklabels = [r"MICE"]
thetalist = np.array([5., 10., 15., 20.]*2)

"""

# Weighted troughs: Redshifts

h=0.7
Runit = 'Mpc'
#Rlist = [0.82106757, 1.64213514, 2.46320271, 3.28427028]
Rlist = [1.64]*4

plotfit = False
Nrows = 1


thetalist = np.array(['10', '6p826'])
samplelist = np.array(['lowZ', 'highZ'])

path_sheardata = 'data2/brouwer/shearprofile/trough_results_final'

path_lenssel = np.array([ ['No_bins_kids_%s_complex/Pmasktheta%s_0p8_inf-delta%s_%s_lw-Wtheta%s'%(samplelist[i], thetalist[i], thetalist[i], delta, thetalist[i]) \
                for i in range(len(thetalist)) ] for delta in ['minf_0', '0_inf'] ])

path_cosmo = 'ZB_0p1_0p9-Om_0p25-Ol_0p75-Ok_0-h_0p7/Rbins10_0p5_20_Mpc/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datatitles = [r'$0.1<z<0.198$', r'$0.198<z<0.3$']
datalabels = ['Weighted troughs', 'Weighted ridges']

plotfilename = '/data2/brouwer/shearprofile/trough_results_final/Plots/troughs_kids_redshifts_weighted'

path_randoms = np.array([ ['No_bins_kids_randoms_complex/Pmasktheta%s_0p8_inf'%theta
                for theta in thetalist ] for delta in ['minf_0', '0_inf'] ])

path_mockdata = 'data2/brouwer/shearprofile/trough_results_final'
path_mocksel = np.array([ ['No_bins_mice_%s_nomask-1/Pmasktheta%s_0p8_inf-delta%s_%s.txt'%(samplelist[i], thetalist[i], thetalist[i], delta) \
                for i in range(len(thetalist)) ] for delta in ['minf_0', '0_inf'] ])
mocklabels = [r"MICE"]





# Randoms

Runit = 'arcmin'

path_sheardata = 'data2/brouwer/shearprofile/trough_results_May'
path_lenssel = [['No_bins_gama_randoms/Pmasktheta5_0p%s_1'%masknum for masknum in range(6,10)]]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datalabels = [r'GAMA ($A_{\rm eff} > %g$)'%n for n in np.arange(0.6, 1.0, 0.1)]

plotfilename = '/data2/brouwer/shearprofile/trough_results_May/Plots/troughs_gama_randoms'




## KiDS vs GAMA

# Fiducial troughs

Runit = 'arcmin'
plotfit = True
thetalist = np.array([5., 5., 5., 5.]) # in arcmin

path_sheardata = 'data2/brouwer/shearprofile/trough_results_final'
path_lenssel = [ [ 'No_bins_%s/Pmasktheta5_0p8_inf-Ptheta5_%s'%(cat, perc) ] \
    for perc in ['0_0p2', '0p8_1'] for cat in ['gama_mice_complex', 'kids_mice_complex'] ]

path_cosmo = 'ZB_0p1_0p9-Om_0p25-Ol_0p75-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datatitles = [r'$\theta_{\rm A} = %g$ arcmin'%theta for theta in thetalist]
datalabels = [r"GAMA: $P(5')<0.2$", r"KiDS: $P(5')<0.2$", r"GAMA: $P(5')>0.8$", r"KiDS: $P(5')>0.8$"]

plotfilename = '/%s/Plots/troughs_gama_kids_mocks_complex'%path_sheardata
Nrows = 1

path_randoms = [ ['No_bins_%s/Pmasktheta5_0p8_inf'%(cat) \
    for perc in ['0_0p2', '0p8_1'] for cat in ['gama_randoms_complex', 'kids_randoms_complex'] ]]

# MICE mocks
path_mockdata = 'data2/brouwer/shearprofile/trough_results_final'
path_mocksel = [ ['No_bins_mice_all_nomask-%g/Pmasktheta5_0p8_inf-Ptheta5_%s.txt'%(cat, perc) \
    for cat in np.arange(16)+1 for perc in ['0_0p2', '0p8_1'] ]]
mocklabels = [r"MICE-GC mocks"]

# SLICS mocks
path_slics = 'slics_mocks_nomask'

slicsfiles = 'Ptheta5-0_0p2.npy'
slicsdata = np.load('/%s/%s/%s'%(path_sheardata, path_slics, slicsfiles))

error_factor = 100./360.3
covfiles_slics = 'Ptheta%s_cov.npy'%(('%g'%theta).replace('.','p'))
covariance_slics = error_factor * np.array((np.load('/%s/%s/%s'%(path_sheardata, path_slics, covfiles_slics))[0]).values())

slicsdata_x = slicsdata[0]
slicsdata_y = np.array(slicsdata[1])

print(slicsdata_x)
print(slicsdata_y)

slicserrors = [np.sqrt(np.diag(covariance_slics[x])) for x in range(len(covariance_slics))]
slicserror_h = np.array(slicserrors)
slicserror_l = np.array(slicserrors)




# Troughs (with different completeness)

Runit = 'arcmin'
plotfit = False

path_sheardata = 'data2/brouwer/shearprofile/trough_results_May'
path_lenssel = [ ['No_bins_kids_all_complex/Pmasktheta5_0p%g_inf-Ptheta5_0_0p2'%(compnum)] for compnum in np.arange(6, 10, 1) ]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datatitles = []
datalabels = [r'comp(A) $> 0.%g$'%(compnum) for compnum in np.arange(6, 10, 1)]
plotfilename = '/data2/brouwer/shearprofile/trough_results_May/Plots/troughs_kids_completeness'
Nrows = 1



# Randoms (with different sizes)

Runit = 'arcmin'
plotfit = False

thetalist = np.array([5.,5.]) # in arcmin

path_sheardata = 'data2/brouwer/shearprofile/trough_results_final'
path_lenssel = [ ['No_bins_%s/Pmasktheta%g_0p8_inf'%(cat,theta) for theta in [thetalist[0]]] \
                                        for cat in ['gama_randoms_complex', 'kids_randoms_complex'] ]
path_cosmo = 'ZB_0p1_0p9-Om_0p25-Ol_0p75-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datatitles = [r'$\theta_{\rm A} = %g$ arcmin'%theta for theta in thetalist]
datalabels = [r'GAMA random signal', r'KiDS random signal']
plotfilename = '/%s/Plots/kids_gama_randoms_complex'%path_sheardata
Nrows = 1

path_mocksel = []



# Randoms (with different completeness)

Runit = 'arcmin'
plotfit = False

path_sheardata = 'data2/brouwer/shearprofile/trough_results_May'
path_lenssel = [ ['No_bins_kids_randoms/Pmasktheta5_0p%g_inf'%(compnum)] for compnum in np.arange(6, 10, 1) ]
path_cosmo = 'ZB_0p1_0p9-Om_0p315-Ol_0p685-Ok_0-h_0p7/Rbins20_2_100_arcmin/shearcovariance'
path_filename = 'No_bins_%s.txt'%(blind)

datatitles = []
datalabels = [r'comp(A) $> 0.%g$'%(compnum) for compnum in np.arange(6, 10, 1)]
plotfilename = '/data2/brouwer/shearprofile/trough_results_May/Plots/troughs_kids_randoms_completeness'
Nrows = 1


"""


Nbins = np.shape(path_lenssel)
Nsize = np.size(path_lenssel)
path_lenssel = np.reshape(path_lenssel, [Nsize])

print('Profiles, Bins:', Nbins)

esdfiles = np.array([('/%s/%s/%s/%s'%(path_sheardata, path_lenssel[i], path_cosmo, path_filename)) \
           for i in range(len(path_lenssel))])

lensIDfiles = np.array([e.replace('_%s.txt'%blind, '_lensIDs.txt') for e in esdfiles])
covfiles = np.array([e.replace('bins_%s.txt'%blind, 'matrix_%s.txt'%blind) for e in esdfiles])

# Importing the shearprofiles and lens IDs
data_x, data_y, error_h, error_l = utils.read_esdfiles(esdfiles)
lensIDs = np.array([np.loadtxt(x) for x in lensIDfiles])


try:
    # Importing the mock shearprofiles
    Nmocks = np.shape(path_mocksel)
    path_mocksel = np.reshape(path_mocksel, np.size(path_mocksel))

    if Nmocks[1] > 5:
        valpha = 0.3
    else:
        valpha = 1.

    esdfiles_mock = np.array([('/%s/%s'%(path_mockdata, path_mocksel[i])) \
           for i in range(len(path_mocksel))])

    data_x_mock, data_y_mock, error_h_mock, error_l_mock = utils.read_esdfiles(esdfiles_mock)

except:
    pass

try:
    print('Import random signal:')

    path_randoms = np.reshape(path_randoms, [Nsize])
    random_esdfile = np.array(['/%s/%s/%s/%s'%(path_sheardata, path_random, path_cosmo, path_filename) for path_random in path_randoms])
    random_data_x, random_data_y, random_error_h, random_error_l = utils.read_esdfiles(random_esdfile)

    # Subtract random signal
    data_y = data_y-random_data_y
    error_h = np.sqrt(error_h**2. + random_error_h**2)
    error_l = np.sqrt(error_l**2. + random_error_l**2)

except:
    print()
    print('No randoms subtracted!')
    print()
    pass
    



# Calculate detection significance

model = np.zeros(len(data_x[0]))

chi2list = np.zeros(Nsize)
chi2covlist = np.zeros(Nsize)
Alist = []
Alist_error = []

for n in range(Nsize):
    
    datax = data_x[n]
    datay = data_y[n]
    
    error = error_h[n]
    covariance = np.array(np.loadtxt(covfiles[n])).T

    chi2list[n] = np.sum(datay**2/error**2)
    chi2covlist[n] = utils.calc_chi2(datay, model, covariance, 1)

    # Signal to noise of the first bin beyond the trough size
    thetaNlist = np.append(thetalist, thetalist)
    
    if 'arcmin' in Runit:
        xmin = thetalist[n]*1.2
        xmax = 70.
    if 'pc' in Runit:
        xmin = Rlist[n]*1.2
        xmax = 10.
        thetalist = Rlist

    xmask = (xmin < datax) & (datax < xmax)
    xwhere = np.array(np.where(xmask))[0]

    # Signal to noise of the Amplitude fit

    # With covariance
    ind = np.lexsort((covariance[3,:], covariance[1,:], covariance[2,:], covariance[0,:]))
    covmatrix = np.reshape(covariance[4][ind], [len(datax), len(datax)])
    covmatrix = covmatrix[int(xwhere[0]):int(xwhere[-1]+1), int(xwhere[0]):int(xwhere[-1]+1)]
        
    A, Acov = optimization.curve_fit(f=trough_model, xdata=datax[xmask], ydata=datay[xmask], p0=[0.], \
    sigma=covmatrix, absolute_sigma=True)
    A = A[0]
    
    Alist = np.append(Alist, A)
    Alist_error = np.append(Alist_error, np.sqrt(Acov[0,0]))
    
    print
    print((datatitles*Nbins[0])[n], datalabels[n/Nbins[1]])
    print('First bin S/N):', datay[xwhere[0]]/error[xwhere[0]])
    print('Amplitude fit S/N', A/np.sqrt(Acov[0,0]))
    print('Amplitude:', A)
    print('Error(Amplitude):', np.sqrt(Acov[0,0]) )
    print

try:
    mock_indices = []
    mock_amplitudes = []
    # Fit to mock profile
    for m in range(len(path_mocksel)):
        
        datax_mock = data_x_mock[m]
        datay_mock = data_y_mock[m]
        
        # Without covariance
        #Avals_mock, Acov_mock = optimization.curve_fit(f=mock_model, xdata=datax_mock[xmask], ydata=datay_mock[xmask], p0=[0., -1.])
        Avals_mock, Acov_mock = optimization.curve_fit(f=trough_model, xdata=datax_mock[xmask], ydata=datay_mock[xmask], p0=[0.])
                
        A_mock = Avals_mock[0]
        #B_mock = Avals_mock[1]
        #mock_indices = np.append(mock_indices, B_mock)
        mock_amplitudes = np.append(mock_amplitudes, A_mock)
        
    mock_troughs = mock_amplitudes[mock_amplitudes < 0.]
    mock_ridges =  mock_amplitudes[mock_amplitudes > 0.]


    print('Mock troughs:')
    print(np.sort(mock_troughs))
    print(np.argsort(mock_troughs))
    print
    print('Mock ridges:')
    print(np.sort(mock_ridges))
    print(np.argsort(mock_ridges))


    #print('Mean mock index:', np.mean(mock_indices))

    # Calculate the detection significance
    chilist, chicovlist = np.sqrt(chi2list), np.sqrt(chi2covlist)

    print('Chi2 (without covariance):', chi2list)
    print('Chi2 (with covariance):', chi2covlist)
except:
    pass

# Create the plot

Ncolumns = int(Nbins[1]/Nrows)

# Plotting the ueber matrix
if Nbins[1] > 1:
    fig = plt.figure(figsize=(Ncolumns*4.,Nrows*3.5))
else:
    fig = plt.figure(figsize=(5,4))
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
            
            if Nbins[0] > 2:
                data_x_plot = data_x[Ndata] * (1.+0.03*Nplot)
            else:
                data_x_plot = data_x[Ndata]
            
            # Plot the fit
            if plotfit:
                xmin = 1.2*thetalist[N]

                model_x = np.linspace(xmin, xmax, 50)
                model_y = trough_model(model_x, Alist[Ndata])
                ax_sub.plot(model_x, model_y, ls='-', color=colors[Nplot])

            if Nsize==Nbins:
                ax_sub.errorbar(data_x_plot, data_y[Ndata], yerr=[error_l[Ndata], error_h[Ndata]], \
                ls='', marker='.', color = colors[Nplot], zorder=3)
            else:
                ax_sub.errorbar(data_x_plot, data_y[Ndata], yerr=[error_l[Ndata], error_h[Ndata]], \
                ls='', marker='.', label=datalabels[Nplot], color = colors[Nplot], zorder=3)
        
        try:
            # Plot mock shearprofiles
            for Nmock in range(Nmocks[1]):

                Ndata = N + Nmock*(Nbins[1])
                
                if Ndata==0:
                    ax_sub.plot(data_x_plot, data_y_mock[Ndata], marker='', ls='-', \
                    color='grey', label=mocklabels[0], alpha=valpha, zorder=1)
                else:
                    ax_sub.plot(data_x_plot, data_y_mock[Ndata], marker='', ls='-', \
                    color='grey', alpha=valpha, zorder=1)
        except:
            pass
        
        #ax_sub.plot(slicsdata_x, slicsdata_y, ls='-', color='blue')
        
        if 'gama' not in path_lenssel[0]:
            # Negative troughs for comparison
            ax_sub.plot(data_x[N], -data_y[N], ls='', marker='.', label='Flipped trough profile',\
            alpha=0.5, color='green')
        
        # Plot the axes and title
    
        # Vertical lines
        ax_sub.axvline(x=1.2*thetalist[N], color='black', ls=':')
        ax_sub.axvline(x=xmax, color='black', ls=':')
        
        ax_sub.axhline(y=0., ls=':', color='black')
        
        ax_sub.xaxis.set_label_position('top')
        ax_sub.yaxis.set_label_position('right')

        ax.tick_params(labelleft='off', labelbottom='off', top='off', bottom='off', left='off', right='off')

        if (N1+1) != Nrows:
            ax_sub.tick_params(axis='x', labelbottom='off')
        if N2 != 0:
            ax_sub.tick_params(axis='y', labelleft='off')
        
        plt.autoscale(enable=False, axis='both', tight=None)
        
        # Define the labels for the plot
        if 'pc' in Runit:
            plt.axis([0.5,20,-2.5,7])
            
            xlabel = r'Radius $R$ (%s/h$_{%g}$)'%(Runit, h*100)
            ylabel = r'ESD $\langle\Delta\Sigma\rangle$ [h$_{%g}$ M$_{\odot}$/pc$^2$]'%(h*100)
            
            ax.xaxis.set_label_coords(0.5, -0.15)
            ax.yaxis.set_label_coords(-0.05, 0.5)
        else:
            #plt.axis([2,110,-1.5e-3,2e-3])
            plt.axis([2,110,-0.9e-3,1.5e-3])

            xlabel = r'Angular separation $\theta$ (arcmin)'
            ylabel = r'Shear $\gamma$'
            if 'random' in datalabels[0]:
                ylabel = r'Random shear $\gamma_0$'
                plt.axis([2,110,-1e-3,1e-3])
                
            if Nbins[1] > 1:
                ax.xaxis.set_label_coords(0.5, -0.025*Ncolumns)
                ax.yaxis.set_label_coords(-0.06*Nrows, 0.5)
            else:
                ax.xaxis.set_label_coords(0.5, -0.07)
                ax.yaxis.set_label_coords(-0.15, 0.5)
        
        plt.xscale('log')
        #plt.yscale('log')
        
        if Nbins[1]>1:
            plt.title(datatitles[N], x = 0.73, y = 0.86, fontsize=16)


# Define the labels for the plot
ax.set_xlabel(xlabel, fontsize=14)
ax.set_ylabel(ylabel, fontsize=14)

#ax.xaxis.label.set_size(17)
#ax.yaxis.label.set_size(17)

handles, labels = ax_sub.get_legend_handles_labels()

# Plot the legend
if Nbins[1] > 1:
    lgd = ax_sub.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.95*Ncolumns, 0.6*Nrows)) # side
else:
    lgd = ax_sub.legend(handles[::-1], labels[::-1], bbox_to_anchor=(0.85, 1.55)) # top
    if 'random' in datalabels[0]:
        plt.legend(handles[::-1], labels[::-1], loc='upper center')


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
for ext in ['pdf']:
    plotname = '%s.%s'%(plotfilename, ext)
    plt.savefig(plotname, format=ext, bbox_inches='tight')
    
print('Written: ESD profile plot:', plotname)

plt.show()
plt.clf


