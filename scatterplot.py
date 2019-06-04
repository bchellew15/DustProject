#load in alphas and make a scatterplot

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

#load wavelengths
hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
wavelength = np.array(hdulist[1].data)  # Angstroms

#load in npy files
alphas1 = np.load('alphas_1d.npy')
alpha_std1 = np.load('alphas_1d_stds.npy')
alphas2 = np.load('alphas_2d.npy')
alpha_std2 = np.load('alphas_2d_stds.npy')
alphas3 = np.load('alphas_iris.npy')
alpha_std3 = np.load('alphas_iris_stds.npy')

plt.subplot(2, 3, 1)
plt.plot(wavelength[50:-100], alphas1[50:-100], c='k') #s=5
plt.plot(wavelength[50:-100], alpha_std1[50:-100], c='k') #s=5
plt.xlim(4830, 5040)
plt.ylim(0, 0.6)
#plt.xlabel("wavelength")
#plt.ylabel("alpha")
plt.title("Original")

plt.subplot(2, 3, 4)
plt.plot(wavelength[50:-100], alphas1[50:-100], c='k') #s=5
plt.plot(wavelength[50:-100], alpha_std1[50:-100], c='k') #s=5
plt.xlim(6530, 6770)
plt.ylim(0, 1.1)
#plt.xlabel("wavelength")
#plt.ylabel("alpha")
plt.title("Original")

plt.subplot(2, 3, 2)
plt.plot(wavelength[50:-100], alphas2[50:-100], c='k') #s=5
plt.plot(wavelength[50:-100], alpha_std2[50:-100], c='k') #s=5
plt.xlim(4830, 5040)
plt.ylim(0, 0.6)
#plt.xlabel("wavelength")
#plt.ylabel("alpha")
plt.title("With Tao")

plt.subplot(2, 3, 5)
plt.plot(wavelength[50:-100], alphas2[50:-100], c='k') #s=5
plt.plot(wavelength[50:-100], alpha_std2[50:-100], c='k') #s=5
plt.xlim(6530, 6770)
plt.ylim(0, 1.1)
#plt.xlabel("wavelength")
#plt.ylabel("alpha")
plt.title("With Tao")

plt.subplot(2, 3, 3)
plt.plot(wavelength[50:-100], alphas3[50:-100], c='k') #s=5
plt.plot(wavelength[50:-100], alpha_std3[50:-100], c='k') #s=5
plt.xlim(4830, 5040)
plt.ylim(0, 0.6)
#plt.xlabel("wavelength")
#plt.ylabel("alpha")
plt.title("With IRIS")

plt.subplot(2, 3, 6)
plt.plot(wavelength[50:-100], alphas3[50:-100], c='k') #s=5
plt.plot(wavelength[50:-100], alpha_std3[50:-100], c='k') #s=5
plt.xlim(6530, 6770)
plt.ylim(0, 1.1)
#plt.xlabel("wavelength")
#plt.ylabel("alpha")
plt.title("With IRIS")


plt.tight_layout()
plt.show()




