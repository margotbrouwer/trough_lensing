#!/usr/bin/python

"Module to determine the local overdensity delta_r within a sphere of radius r."

# Import the necessary libraries
import astropy.io.fits as pyfits
import gc
import numpy as np
import sys
import os
import time
from glob import glob

from astropy import constants as const, units as u
from astropy.coordinates import SkyCoord
from collections import Counter

import trough_modules_all as utils
from matplotlib import pyplot as plt


Nsrc = 50
dist = 10.

RAs1 = np.ones(Nsrc)*(10.)
RAs2 = np.ones(Nsrc)*(10.+dist)
DECs = np.linspace(-90, 90, Nsrc)

coords1 = SkyCoord(ra=RAs1*u.deg, dec=DECs*u.deg)
coords2 = SkyCoord(ra=RAs2*u.deg, dec=DECs*u.deg)

idx, d2d, d3d = coords1.match_to_catalog_sky(coords2)

print(idx)
print(d2d.degree/dist)

galRA = RAs1
srcRA = RAs2
galDEC = DECs
srcDEC = DECs

srcR = np.degrees(np.arccos(np.cos(np.radians(galDEC))*\
                np.cos(np.radians(srcDEC))*\
                np.cos(np.radians(galRA-srcRA))+\
                np.sin(np.radians(galDEC))*\
                np.sin(np.radians(srcDEC))))

print(srcR/dist)
print(np.cos(np.radians(DECs)))

#plt.plot(DECs, d2d.degree/dist)
#plt.plot(DECs, srcR/dist)
plt.plot(DECs, np.cos(np.radians(DECs)))

plt.show()
