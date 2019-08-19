#look at updated BOSS data files
#these arrays have already been padded
#BOSS_test.py has code to pad arrays and load padded arrays

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys #for numpy threshold

filenames = ['skyfibers_lam0.fits', 'skyfibers_lam1.fits', 'skyfibers_lam2.fits', 'skyfibers_lam3.fits', 'skyfibers_lam4.fits', \
             'skyfibers_lam5.fits', 'skyfibers_lam6.fits', 'skyfibers_lam7.fits', 'skyfibers_lam8.fits', 'skyfibers_lam9.fits']

'''
hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[0])
print(hdulist[0].data)
print(np.unique(hdulist[2].data))

hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[1])
print(hdulist[0].data)
print(np.unique(hdulist[2].data))
'''

for f in filenames:
    hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + f)
    #print(hdulist[1].data)
    #print("mean ivar:")
    #print(np.mean(hdulist[1].data))

#column names;
#0: flambda
#1: ivar
#2: lambda
#3: RA
#4: dec
#5: mjd
#6: plate
#7: fiber

#compare locations (RA, DEC) to the ones I used for SFD:
#correct RA:
'''
hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[0])
dec = hdulist[4].data
#RA I used:
my_dec = np.loadtxt("/Users/blakechellew/Documents/DustProject/BrandtFiles/BOSS_locations.txt")[:,1]

for i in range(len(dec)):
    if dec[i] != my_dec[i]:
        print("error: values don't match")

plt.plot(dec, my_dec, 'k.')
plt.plot(np.arange(360), np.arange(360), 'r')
plt.show()
'''



#check whether fibers on same plate start with same lambda
#which means: check whether fibers on the same plate have same padding
'''
hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[0])
plate = hdulist[6].data
mjd = hdulist[5].data
plate = 10000*plate + mjd%10000 #plate is now a unique identifier

#for each unique plates, get all the lambdas, see if any don't match
unique_plates = np.unique(plate, return_index=True)[1] #2512 unique plates

for i in unique_plates:
    p = hdulist[6].data[i]
    d = mjd[i]
    idces =  np.nonzero((hdulist[6].data == p) * (mjd == d))[0] #indices of rows with this p and d

    nonzero_idces = [np.nonzero(hdulist[0].data[j])[0] for j in idces]
    paddings = [k[0] for k in nonzero_idces if k.size != 0]
    
    if len(np.unique(paddings)) != 1:
        print("fail")
'''

#histograms of ivar:
'''
for f in filenames:
    hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + f)

    ivar = hdulist[1].data.flatten()
    plt.hist(ivar, bins=50, range=(0, 30))
    plt.show()
'''

#histograms of flam::

for f in filenames:
    hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + f)

    flam = hdulist[0].data.flatten()
    plt.hist(flam, bins=50, range=(-3, 3))
    plt.show()


