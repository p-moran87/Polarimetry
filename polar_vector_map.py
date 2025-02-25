# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:38:25 2025

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
img_60 = fits.open(file_60[0])['SCI',1].data

file_120 = glob.glob('CRAB_HST_ACS_HESTER/hst_10526_01_acs_wfc_f606wpol120v/*drz.fits')
img_120 = fits.open(file_120[0])['SCI',1].data

# Get WCS from header
hdul = fits.open(file_0[0])
wcs = WCS(hdul[1].header)
#print(hdul[0].header["PA_V3"])
pa_v3 = 87.236
chi = -38.200

# Loop through images - scan to get photometry in rectangle
pol_map = []
pixX = []
pixY = []

for i in range(16,2046,32):
    for j in range(16,2048,32):
        aperture = RectangularAperture((i,j), 32, 32)
        pol_data = polarization_tools.Polarization(ApertureStats(img_0, aperture).median, 
                    ApertureStats(img_60, aperture).median, 
                    ApertureStats(img_120, aperture).median, 'F606W', 'WFC', pa_v3)
        
        pol_data.calc_stokes()
        angle = np.degrees(0.5*np.arctan2(pol_data.stokes_u, pol_data.stokes_q)) + pa_v3 + chi
        angle = np.mod(angle, 180)
        
        if pol_data.stokes_i > 1.0: 
            pol_frac = np.sqrt(pol_data.stokes_q**2 + pol_data.stokes_u**2)/pol_data.stokes_i
            if pol_frac <= 1.0:
                pol_map.append([i, j, pol_frac*100.0, angle])
                # Angles are measured E of Astronomical N
                pixX.append(pol_frac*np.cos(np.deg2rad(angle + 90.0))) # X-vector 
                pixY.append(pol_frac*np.sin(np.deg2rad(angle + 90.0))) # Y-vector    
            else:
                pol_map.append([i, j, 0.0, angle])
                pixX.append(0.0)
                pixY.append(0.0)
        else:
            pol_map.append([i, j, 0.0, angle])
            pixX.append(0.0)
            pixY.append(0.0)
        
# First, we plot the background image
fig = plt.figure(figsize=(10,10))
i_plot = fig.add_subplot(111, projection=wcs)
plt.imshow(img_0, vmin=-0.5, vmax=1.2, origin='lower')
overlay = i_plot.get_coords_overlay('fk5')
overlay.grid(color='white', ls='dotted')

# ranges of the axis
xx0, xx1 = i_plot.get_xlim()
yy0, yy1 = i_plot.get_ylim()

# re-binned number of points in each axis
factor = [32,32]
nx_new = 2048 // factor[0]
ny_new = 2048 // factor[1]

# These are the positions of the quivers
X,Y = np.meshgrid(np.linspace(xx0,xx1,nx_new,endpoint=True),
                  np.linspace(yy0,yy1,ny_new,endpoint=True))

# keyword arguments for quiverplots
quiveropts = dict(headlength=0, headwidth=1, pivot='middle')
i_plot.quiver(X, Y, pixX, pixY, scale=10, **quiveropts)

# Draw a horiz. line as ref. distance & add label above line
x_start, x_end = 2046, 2096
i_plot.plot([x_start, x_end], [50, 50], 'w-', lw=2)  # white line at y=0
i_plot.text((x_start + x_end) / 2, 100, r'$50\%$', 
            fontsize=18, color='white', ha='center')

plt.savefig('Polarisation Map')
plt.show(block=False)

pol_array = np.array(pol_map)
print("Mean Pol %:", np.mean(pol_array[:,2]))
print("Mean Angle (Deg):", np.mean(pol_array[:,3]))

#histogram of polarization %
plt.hist(pol_array[:,2],24,range=(0,100))
plt.xlabel('Polarisation (%)')
plt.ylabel('Frequency')
plt.title('Histogram of Polarisation %')
plt.savefig('Histogram Pol %')
plt.show(block=False)

# histogram of polarization angle
plt.hist(pol_array[:,3],18, range=(0,180))
plt.xlabel('Angle ($^\circ$)')
plt.ylabel('Frequency')
plt.title('Histogram of Position Angle ($^\circ$)')
plt.savefig('Histogram Angle')
plt.show(block=False)