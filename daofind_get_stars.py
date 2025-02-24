# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:35:52 2025

@author: PaulMoran
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from acstools import polarization_tools
from photutils.aperture import RectangularAperture, ApertureStats
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

# Load and display the images
files = glob.glob('CRAB_HST_ACS_HESTER/f550m/*drz.fits')
ref_image = files[0]
image_data = fits.open(ref_image)['SCI',1].data

# Get WCS from header
hdul = fits.open(ref_image)
wcs = WCS(hdul[1].header)
  
data = hdul[1].data
mean, median, std = sigma_clipped_stats(data, sigma=3.0)  
print(np.array((mean, median, std)))  

daofind = DAOStarFinder(fwhm=2.0, threshold=400.*std)  
sources = daofind(data - median)  
for col in sources.colnames:  
    if col not in ('id', 'npix'):
        sources[col].info.format = '%.2f'  # for consistent table output
sources.pprint(max_width=76) 