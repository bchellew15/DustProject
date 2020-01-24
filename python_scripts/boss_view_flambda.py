#look at flambda for BOSS
#goal: see if any of the spectra look strange

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

filenames = ['skyfibers_lam0.fits', 'skyfibers_lam1.fits', 'skyfibers_lam2.fits', 'skyfibers_lam3.fits', 'skyfibers_lam4.fits', \
                 'skyfibers_lam5.fits', 'skyfibers_lam6.fits', 'skyfibers_lam7.fits', 'skyfibers_lam8.fits', 'skyfibers_lam9.fits']
num_files = len(filenames)

#get list of plates and fiber ids
hdulist = fits.open("/Volumes/TOSHIBA/Dust_Overflow/" + filenames[0])
plate = hdulist[6].data
fiber_id = hdulist[7].data
flambda = hdulist[0].data

#flambdas_on_plate = []
#ivars_on_plate = []

#make plots:
for q in range(plate.shape[0]):
    if plate[q] == 4900: #started with 5896
        flambda = np.zeros(0)
        ivar = np.zeros(0)
        for f in filenames:
            flambda = np.append(flambda, fits.open("/Volumes/TOSHIBA/Dust_Overflow/" + f)[0].data[q])
            ivar = np.append(ivar, fits.open("/Volumes/TOSHIBA/Dust_Overflow/" + f)[1].data[q])
        wavelength = 10**( min(hdulist[2].data) + np.arange(flambda.shape[0])*1e-4 ).astype('float32')
        plt.plot(wavelength, flambda, 'k', drawstyle='steps')
        plt.plot(wavelength, 1/np.sqrt(ivar), 'r', drawstyle='steps')
        plt.title(fiber_id[q])
        plt.ylim(-5, 5)
        plt.xlim(7800, 8200)
        plt.show()

        #print(np.mean(flambda), fiber_id[q])
        #flambdas_on_plate.append(np.mean(flambda))

        #ivars_on_plate.append(np.mean(ivar))
        #print(np.mean(ivar), fiber_id[q])
        

#plt.hist(flambdas_on_plate, 20)
#plt.hist(ivars_on_plate, 20)
#plt.show()


