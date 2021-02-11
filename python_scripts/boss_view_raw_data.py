# look at flambda and make plots.
# look at indices of plates and fibers.
# view a very basic plot of alpha vs wavelength for boss (given alpha)
# goal: look at 100 micron data by fiber / plate
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

i100_old = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/SFD_i100_at_BOSS_locations.txt")[:,2]
i100 = np.load("/Volumes/TOSHIBA/Dust_Overflow/i100_tao_boss_iris.npy", mmap_mode='r')

wavelength_boss = np.load('../alphas_and_stds/wavelength_boss.npy')

# coords = np.loadtxt("/Users/blakechellew/Documents/DustProject/BrandtFiles/BOSS_locations_galactic.txt")
# see "sky_locations.py" for sdss locations

i100 = i100[:, -1]
p1 = plate[2384]
p2 = plate[2385]
print(p1)
print(p2)
relevant_fibers = np.where(plate == p1)[0]
print(relevant_fibers)
rel_i100 = i100[relevant_fibers]
print(rel_i100)
print(np.mean(rel_i100))
print(i100[2389])
print(i100[2390])
print(i100[2384])
print(i100[2385])





exit(0)

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


