#!/usr/bin/python

import numpy as np
import sys
import subprocess as sub
import shlex

Rmin = 1.
Rmax = 300.
Nbins = 25

Rbins = np.logspace(np.log10(Rmin), np.log10(Rmax), Nbins+1)

Rlim_low = Rbins[0:-1]
Rlim_high = Rbins[1:len(Rbins)]
Rcenters = Rlim_low + np.diff(Rbins)/2.

print(Rcenters)

xmin = 6.
xmax = 60.
xmask = (xmin < Rcenters) & (Rcenters < xmax)

filename = 'Rbins_troughs_cut.txt'
file_header = 'Radial bins of the troughs (in arcmin), cut within 6 and 60 arcmin'

Rbins_tot = np.array([Rlim_low[xmask], Rcenters[xmask], Rlim_high[xmask]])

np.savetxt(filename, Rbins_tot.T, delimiter='    ',  fmt='%.18g', header = file_header)
print('Written:', filename)
