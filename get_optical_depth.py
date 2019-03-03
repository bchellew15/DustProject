from ebv_to_tao import f99
import numpy as np
from astropy.io import fits

#using E(B-V) values calculated from the SFD map
#see get_dustmaps.py

hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
wavelength = np.array(hdulist[1].data)  # Angstroms

ebv = np.load("/Users/blakechellew/Documents/DustProject/ebv.npy")

taos = np.array([f99(wavelength, item) for item in ebv])
print("shape: ", taos.shape)

np.save("/Users/blakechellew/Documents/DustProject/taos.npy", taos)





