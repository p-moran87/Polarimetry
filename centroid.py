# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:17:38 2025

@author: PaulMoran
"""

import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.modeling.models import Gaussian2D
from astropy.visualization import simple_norm
from photutils.centroids import centroid_quadratic
from photutils.datasets import make_noise_image
from photutils.profiles import RadialProfile
from astropy.table import Table
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

xycen = centroid_quadratic(image_data, xpeak=1117.5, ypeak=823.5)
edge_radii = np.arange(20)
rp = RadialProfile(image_data, xycen, edge_radii, mask=None)

print(rp.radius)
print(rp.profile)

plt.plot(rp.radius, rp.profile)