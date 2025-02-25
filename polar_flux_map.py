# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:27:12 2025

@author: PaulMoran
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from acstools import polarization_tools
from photutils.aperture import RectangularAperture, ApertureStats
from astropy.wcs import WCS

# Load and display the images
file_0 = glob.glob('CRAB_HST_ACS_HESTER/hst_10526_01_acs_wfc_f606wpol0v/*drz.fits')
img_0 = fits.open(file_0[0])['SCI',1].data

file_60 = glob.glob('CRAB_HST_ACS_HESTER/hst_10526_01_acs_wfc_f606wpol60v/*drz.fits')
ref_img_60 = file_60[0]
img_60 = fits.open(ref_img_60)['SCI',1].data

file_120 = glob.glob('CRAB_HST_ACS_HESTER/hst_10526_01_acs_wfc_f606wpol120v/*drz.fits')
ref_img_120 = file_120[0]
img_120 = fits.open(ref_img_120)['SCI',1].data

# Get WCS from header
hdul = fits.open(file_0[0])
wcs = WCS(hdul[1].header)
#print(hdul[0].header["PA_V3"])
pa_v3 = 87.236
chi = -38.200

# Loop through the images scan the image to get photometry in rectangle
pol_map = []

for i in range(4,2046,8):
    for j in range(4,2048,8):
        aperture = RectangularAperture((i,j), 8, 8)
        pol_data = polarization_tools.Polarization(ApertureStats(img_0, aperture).median, 
                    ApertureStats(img_60, aperture).median, 
                    ApertureStats(img_120, aperture).median, 'F606W', 'WFC', pa_v3)
        
        pol_data.calc_stokes()    
        if pol_data.stokes_i > 1.0: 
            pol_frac = np.sqrt(pol_data.stokes_q**2 + pol_data.stokes_u**2)/pol_data.stokes_i
            pol_map.append([i, j, pol_frac * 5000.0])
        else:
            pol_map.append([i, j, 0.0])
        
pol_array = np.array(pol_map)

x = pol_array[:,0]
y = pol_array[:,1]
pf = pol_array[:,2]

X = np.asarray(x)
Y = np.asarray(y)
PF = np.asarray(pf)

Xu = np.unique(x)
Yu = np.unique(y)

PFimg = PF[np.lexsort((X, Y))].reshape(Xu.size, Yu.size)
plt.imshow(PFimg, vmin=-0.5, vmax=1.2, origin='lower')
plt.show()
plt.imsave('polar_flux_map.png', PFimg, origin='lower')