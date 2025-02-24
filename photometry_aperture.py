# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 12:35:48 2025

@author: PaulMoran
"""

import glob
import matplotlib.pyplot as plt
from acstools import acszpt
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.visualization import simple_norm
from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAnnulus
from photutils.aperture import CircularAperture
from astropy.stats import sigma_clipped_stats

# Load and display the images
files = glob.glob('CRAB_HST_ACS/*flt.fits')
ref_image = files[0]
ref_fits = fits.open(ref_image)
image_data = fits.open(ref_image)['SCI',1].data

norm1 = simple_norm(image_data,stretch='log',min_cut=0,max_cut=50000)

plt.imshow(image_data, origin='lower',
                      norm=norm1,cmap='gray')
#plt.gca().tick_params(labelcolor='none',axis='both',color='none')
plt.show()


'Perform the Photometry on the target in F550M'

# Get the phtometric zeropoints for the filter for the observation date
q = acszpt.Query(date='2005-09-06', detector='WFC', filt='F550M')
filter_zpt = q.fetch()

print(filter_zpt)

positions = [1317.504,1044.280]

# Photometry with circular aperture and annulii

apertures = CircularAperture(positions, r = 6)
rawflux = aperture_photometry(image_data, apertures)

annulus_apertures = CircularAnnulus(positions, r_in = 6, r_out = 100)
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

#rawflux['final_phot'] = rawflux['aperture_sum'] - 150/1150

# Apply the correction from 0.2" to infinity.
# For F555W, the correction is 0.841
correction_inf = 0.841
flux_inf = rawflux['final_phot'] / correction_inf


# Now convert instrumental fluxes to physical fluxes and magnitudes.
# F_lambda is the flux density in units of erg/sec/cm^2/Angstrom.
#F_lambda = flux_inf * filter_zpt['PHOTFLAM']
m_st = -2.5 * np.log10(flux_inf/1150) + filter_zpt['STmag'][0].value
m_ab = -2.5 * np.log10(flux_inf/1150) + filter_zpt['ABmag'][0].value
m_vega = -2.5 * np.log10(flux_inf/1150) + filter_zpt['VEGAmag'][0].value

# Assemble the values into an Astropy Table. Note that we could
# attach units to these columns, however advanced Astropy
# Tables use is outside the scope of this example.
phot_table = Table({'Measured Flux': rawflux['final_phot'],
                    'ST Mag': m_st, 'AB Mag': m_ab, 'Vega Mag': m_vega}, 
                   names=['Measured Flux','ST Mag', 'AB Mag', 'Vega Mag'])

print(phot_table)











