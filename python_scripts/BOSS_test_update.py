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

#open all the files:
'''
for f in filenames:
    hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + f)
'''

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
'''
for f in filenames:
    hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + f)

    flam = hdulist[0].data.flatten()
    plt.hist(flam, bins=50, range=(-3, 3))
    plt.show()
'''

#look at means over plates:
'''
#i100 = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/SFD_i100_at_BOSS_locations.txt")[:,2].astype('float32')
i100 = np.load("/Users/blakechellew/Documents/DustProject/IRIS/iris_i100_at_boss.npy", mmap_mode='r').astype('float32')

hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[0])
plate = hdulist[6].data
mjd = hdulist[5].data
plate = 10000*plate + mjd%10000 #plate is now a unique identifier

#for each unique plates, get all the lambdas, see if any don't match
unique_plate_idces = np.unique(plate, return_index=True)[1] #2512 unique plates

for i in unique_plate_idces:
    p = hdulist[6].data[i]
    d = mjd[i]
    idces =  np.nonzero((hdulist[6].data == p) * (mjd == d))[0] #indices of rows with this p and d

    #i100 on this plate:
    i100_sub = i100[idces]
    i100_mean = np.mean(i100_sub)

    if i100_mean > 10:
        print("plate", plate[i])
'''

#compare xxsig and yxsig of SFD and IRIS
'''
xxsig = np.load("xxsig_82119.npy")
xxsig_iris = np.load("xxsig_iris_1d_82119.npy")
yxsig = np.load("yxsig_82119.npy")
yxsig_iris = np.load("yxsig_iris_1d_82119.npy")

xx_difs = xxsig_iris - xxsig

xx_difs.sort()

print(xx_difs[0])
print(xx_difs[-1])

yx_difs = yxsig_iris - yxsig
sorted_idx = np.argsort(yx_difs)[0]
print("index:", sorted_idx)
print("iris at index:", yxsig_iris[sorted_idx])
print("sfd at index:", yxsig[sorted_idx])

print("iris at another index:", yxsig_iris[sorted_idx-1])
print("sfd at another index:", yxsig[sorted_idx-1])
print("iris at another index:", yxsig_iris[sorted_idx+1])
print("sfd at another index:", yxsig[sorted_idx+1])

yx_difs.sort()
print(yx_difs[0:10])
print(yx_difs[-1])
'''


#compare SFD and IRIS at same index:
i100_sfd = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/SFD_i100_at_BOSS_locations.txt")[:,2].astype('float32')
i100_iris = np.load("/Users/blakechellew/Documents/DustProject/IRIS/iris_i100_at_boss.npy", mmap_mode='r').astype('float32')
diffs = i100_iris - i100_sfd
print(diffs)

print("mean SFD:", np.mean(i100_sfd))
print("mean IRIS:", np.mean(i100_iris))

neg_diffs = diffs[diffs<0]
big_diffs = diffs[diffs>2]

print("sfd value:", i100_sfd[3178])
print("iris value:", i100_iris[3178])



print(diffs[3178])

print("neg_diffs")
print(np.sort(neg_diffs))
print("big_diffs")
print(np.sort(big_diffs))


plt.hist(diffs, bins=50, range=(0, 2))
plt.title("IRIS - SFD i100")
plt.show()

