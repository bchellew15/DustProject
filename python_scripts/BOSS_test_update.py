#look at update BOSS data files
#also generate padded flambda and ivar

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys #for numpy threshold

filenames = ['skyfibers_lam0.fits', 'skyfibers_lam1.fits', 'skyfibers_lam2.fits', 'skyfibers_lam3.fits', 'skyfibers_lam4.fits', \
             'skyfibers_lam5.fits', 'skyfibers_lam6.fits', 'skyfibers_lam7.fits', 'skyfibers_lam8.fits', 'skyfibers_lam9.fits']

hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[0])
print(hdulist[0].data)
print(np.unique(hdulist[2].data))

hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[1])
print(hdulist[0].data)
print(np.unique(hdulist[2].data))

#column names;
#0: flambda
#1: ivar
#2: lambda
#3: RA
#4: dec
#5: mjd
#6: plate
#7: fiber
