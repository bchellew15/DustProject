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
plate0 = hdulist[6].data
mjd = hdulist[5].data
fiber_id = hdulist[7].data
flambda = hdulist[0].data

plate = 10000 * plate0 + mjd % 10000 #plate is now a unique identifier

i100_old = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/SFD_i100_at_BOSS_locations.txt")[:,2]
i100 = np.load("/Volumes/TOSHIBA/Dust_Overflow/i100_tao_boss_iris.npy", mmap_mode='r')

wavelength_boss = np.load('../alphas_and_stds/wavelength_boss.npy')

# coords = np.loadtxt("/Users/blakechellew/Documents/DustProject/BrandtFiles/BOSS_locations_galactic.txt")
# see "sky_locations.py" for sdss locations

# check i100 variance on certain plates:

# get plate numbers for masked plates

print("plate numbers")
for i in [36, 938, 1509, 1265, 1786, 2141, 2380, 2383, 2388]:
    plate_num = plate0[np.where(plate == np.unique(plate)[i])]
    print(plate_num)
exit(0)

print("check 1")

i100_1d = np.copy(i100[:, 2000])

print("check 2")

# masking:
i100_1d[i100_old > 10] = np.nan
for p in np.unique(plate):
    if np.mean(i100_old[plate == p]) > 10:
        i100_1d[plate == p] = np.nan
i100_1d[plate == np.unique(plate)[36]] = np.nan
i100_1d[plate == np.unique(plate)[938]] = np.nan
i100_1d[plate == np.unique(plate)[1509]] = np.nan
i100_1d[plate == np.unique(plate)[1265]] = np.nan
i100_1d[plate == np.unique(plate)[1786]] = np.nan
i100_1d[plate == np.unique(plate)[2141]] = np.nan
i100_1d[plate == np.unique(plate)[2380]] = np.nan
i100_1d[plate == np.unique(plate)[2383]] = np.nan
i100_1d[plate == np.unique(plate)[2388]] = np.nan

print("check 3")

variances = np.zeros((len(np.unique(plate))))
print(variances.shape)
for i in range(len(np.unique(plate))):
    if i % 50 == 0:
        print("PROGRESS:", i)
    i100_plate = i100_1d[plate == np.unique(plate)[i]]
    variance = np.nanvar(i100_plate)
    variances[i] = variance
mean_variance = np.nanmean(variances)
print("mean:", mean_variance)
print("mean std:", np.sqrt(mean_variance))
med_variance = np.nanmedian(variances)
print("median:", med_variance)
print("median std:", np.sqrt(med_variance))
print("std:", np.nanstd(variances))
plt.hist(variances)
plt.show()

exit(0)


# view spectra for certain fibers:
"""
print("plate id for idx 3178:")
plate_3178 = plate[3178]
idx_3178 = np.argwhere(np.unique(plate) == plate_3178)
print(idx_3178)
exit(0)
"""

def view_fiber(idx):
    flambda_fiber = np.array([])
    ivar_fiber = np.array([])
    for i in range(len(filenames)):
        hdulist = fits.open("/Volumes/TOSHIBA/Dust_Overflow/" + filenames[i])
        flambda_fiber = np.append(flambda_fiber, hdulist[0].data[idx])
        ivar_fiber = np.append(ivar_fiber, hdulist[1].data[idx])
    print(np.mean(ivar_fiber[(wavelength_boss > 5000) & (wavelength_boss < 9000)]))
    plt.plot(wavelength_boss, flambda_fiber, 'k')
    plt.plot(wavelength_boss, 1 / np.sqrt(ivar_fiber), 'r')
    plt.ylim(-7, 7)
    # plt.xlim(6200, 6600)  # TEMP just to check one plate
    plt.show()

# look at masked plates and fibers
indices = np.argwhere(plate == np.unique(plate)[1265])
for i in indices:
    view_fiber(i[0])
"""
indices = np.argwhere(plate == np.unique(plate)[1265])
for i in indices:
    view_fiber(i[0])
indices = np.argwhere(plate == np.unique(plate)[1786])
for i in indices:
    view_fiber(i[0])
indices = np.argwhere(plate == np.unique(plate)[938])
print("fiber ids:")
for i in indices:
    view_fiber(i[0])
"""

"""
# a random plate that should be fine:
indices = np.argwhere(plate == np.unique(plate)[1000])
for i in indices:
    view_fiber(i[0])
"""
exit(0)

view_fiber(3000)  # a typical spectrum
view_fiber(3178)  # bad data
idx = np.argwhere((plate == np.unique(plate)[938]) * (fiber_id == 136))[0][0]  # was fiber 252
view_fiber(idx)  # this one was also masked

exit(0)

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


