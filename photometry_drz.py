# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 12:35:48 2025

@author: PaulMoran
"""

import glob
import matplotlib.pyplot as plt
from acstools import acszpt
import numpy as np
#from astropy.table import Table
from astropy.io import fits
from astropy.visualization import simple_norm
from acstools import polarization_tools
from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAnnulus
from photutils.aperture import CircularAperture, ApertureStats
from astropy.stats import sigma_clipped_stats

# Load and display the images
files = glob.glob('CRAB_HST_ACS_HESTER/f550m/*drz.fits')
ref_image = files[0]
image_data = fits.open(ref_image)['SCI',1].data

# Load the polarization images: 0, 60, 120
file_0 = glob.glob('CRAB_HST_ACS_HESTER/hst_10526_01_acs_wfc_f606wpol0v/*drz.fits')
ref_img_0 = file_0[0]
img_0 = fits.open(ref_img_0)['SCI',1].data

file_60 = glob.glob('CRAB_HST_ACS_HESTER/hst_10526_01_acs_wfc_f606wpol60v/*drz.fits')
ref_img_60 = file_60[0]
img_60 = fits.open(ref_img_60)['SCI',1].data

file_120 = glob.glob('CRAB_HST_ACS_HESTER/hst_10526_01_acs_wfc_f606wpol120v/*drz.fits')
ref_img_120 = file_120[0]
img_120 = fits.open(ref_img_60)['SCI',1].data

# Display the F550M CLEAR image
norm1 = simple_norm(image_data,stretch='log',min_cut=0,max_cut=10)
plt.imshow(image_data, vmin=-0.5, vmax=1.2, origin='lower')
plt.show()

# Display the F606W POL_0 image
plt.imshow(img_0, vmin=-0.5, vmax=1.2, origin='lower')
plt.show()

'Perform the Photometry on the target'

# Get the phtometric zeropoints for the filter for the observation date
q = acszpt.Query(date='2005-09-06', detector='WFC', filt='F550M')
filter_zpt = q.fetch()

print(filter_zpt)

positions = [1120.000,825.00]
sky_position = [250.0, 1500.0]

# Photometry with circular aperture and annulii
apertures = CircularAperture(positions, r = 4)
rawflux = aperture_photometry(image_data, apertures)

annulus_apertures = CircularAnnulus(sky_position, r_in = 4, r_out = 15)
annulus_masks = annulus_apertures.to_mask(method = 'center')

# Get the sky background from the outer annulus using median
bkg_median = []
annulus_data = annulus_masks.multiply(image_data)
annulus_data_1d = annulus_data[annulus_masks.data > 0]
_, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)

bkg_median.append(median_sigclip)
bkg_median = np.array(bkg_median)

#Perform background subtraction
rawflux['annulus_median'] = bkg_median / annulus_apertures.area
rawflux['aper_bkg'] = bkg_median * apertures.area
rawflux['final_phot'] = rawflux['aperture_sum'] - rawflux['aper_bkg']

# Apply the correction from 0.2" to infinity.
# For F555W, the correction is 0.841
correction_inf = 0.841
flux_inf = rawflux['final_phot'] / correction_inf

# Now convert instrumental fluxes to physical 
# fluxes and magnitudes along with correction for WFC for 0.5" of 0.914.
m_st = -2.5 * np.log10(flux_inf) + filter_zpt['STmag'][0].value - 0.914
m_ab = -2.5 * np.log10(flux_inf) + filter_zpt['ABmag'][0].value - 0.914
m_vega = -2.5 * np.log10(flux_inf) + filter_zpt['VEGAmag'][0].value - 0.914

print(m_st, m_ab, m_vega)

# PA_V3 angle
pa_v3 = 87.236
chi = -38.200

pol_aperture = CircularAperture(positions, r = 5)
sky_aperture = CircularAperture(sky_position, r = 5)

pol_data = polarization_tools.Polarization(ApertureStats(img_0, pol_aperture).sum  
                        - ApertureStats(img_0, sky_aperture).median, 
            ApertureStats(img_60, pol_aperture).sum 
                    - ApertureStats(img_60, sky_aperture).median, 
            ApertureStats(img_120, pol_aperture).sum 
                - ApertureStats(img_120, sky_aperture).median, 'F606W', 'WFC', pa_v3)

pol_data.calc_stokes()
pol_data.calc_polarization()
angle = 0.5*np.arctan2(pol_data.stokes_u, pol_data.stokes_q)*180.0/np.pi + pa_v3 + chi 
pol_frac = np.sqrt(pol_data.stokes_q**2 + pol_data.stokes_u**2)/pol_data.stokes_i

print(f'Stokes I = {pol_data.stokes_i}, Stokes Q = {pol_data.stokes_q}, Stokes U = {pol_data.stokes_u}')

print(f'polarization = {pol_frac:0.2%}, Position Angle = {pol_data.angle:0.5}')

# Get header data
#hdul = fits.open(ref_img_0)
#print(hdul[0].header["PA_V3"])







