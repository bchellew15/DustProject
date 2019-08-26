import numpy as np
from time import time
from astropy.io import fits
import matplotlib.pyplot as plt

'''
#load data
hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
plate = np.array(hdulist[0].data)
wavelength = np.array(hdulist[1].data)  # Angstroms
ivar = np.array(hdulist[3].data)        # inverse variance, units of 1/flambda^2
ivar *= (ivar > 0) #correct the negative values

#x1 = np.random.random_sample(91847)
#x1 = np.random.random_sample((91847, 4000))
#x2 = np.zeros(x1.shape)
'''

'''
#check time of averaging:

t1 = time()
#new method:
boundaries = np.sort(np.unique(plate, return_index=True)[1])
for idx, b in enumerate(boundaries[:-1]):
    avgs = np.mean(x1[boundaries[idx]:boundaries[idx+1]], axis=0) #mean across plates, not wavelength
    x2[boundaries[idx]:boundaries[idx+1]] = avgs
t2 = time()
print("check 1: ", t2-t1)


#old method:
plate_numbers = np.unique(plate)
for p in plate_numbers:
    mask = np.isin(plate, p)
    mask = mask.reshape(len(mask), 1)
    mask = np.tile(mask, 4000)
    vals_on_plate = x1[mask]
    avg = np.mean(vals_on_plate)
    np.putmask(x2, mask, avg)
t3 = time()
print("check 2: ", t3-t2)
'''

'''
#check time for ivar conversion:

t1 = time()
#unit conversion for ivar:
for i in range(0, len(wavelength)):
    ivar[:,i] /= pow(wavelength[i], 2)

t2 = time()
print("check 1: ", t2-t1)
    
ivar /= np.power(wavelength, 2)

t3 = time()
print("check 2: ", t3-t2)
'''

'''
#check time of alpha

x = np.random.random_sample((91847, 4000))
y = np.random.random_sample((91847, 4000))
#x = np.random.random_sample((91847, 4000))
#x2 = np.zeros(x1.shape)


t1 = time()

#calculate alpha
xx = np.multiply(x, x)

t2 = time()
print("check1: ", t2-t1)

yx = np.multiply(y, x)

t3 = time()
print("check2: ", t3-t2)

yxsig = np.multiply(yx, ivar)

t4 = time()
print("check3: ", t4-t3)

xxsig = np.multiply(xx, ivar)

t5 = time()
print("check4: ", t5-t4)

sums1 = np.sum(yxsig, axis=0)

t6 = time()
print("check5: ", t6-t5)

sums2 = np.sum(xxsig, axis=0)

t7 = time()
print("check6: ", t7-t6)

alphas = np.divide(sums1, sums2)

t8 = time()
print("check7: ", t8-t7)

alpha_std = np.sqrt(1/sums2)

t9 = time()
print("check8: ", t9-t8)
'''
'''
t1 = time()

alphas = np.zeros(wavelength.shape)
alpha_std = np.zeros(wavelength.shape)
for i, lmda in enumerate(wavelength):
    ylambda = y[:,i]
    varlambda = ivar[:,i]
    yx = np.multiply(ylambda, x)
    yxsig = np.multiply(yx, varlambda)
    sum1 = np.sum(yxsig)
    xx = np.multiply(x, x)
    xxsig = np.multiply(xx, varlambda)
    sum2 = np.sum(xxsig)
    alpha = sum1 / sum2
    alphas[i] = alpha
    alpha_std[i] = np.sqrt(1/sum2)
t2 = time()
print("check1: ", t2-t1)
'''

'''
#create array of small taos:
taos_small = np.full((91847, 4000), 0.01)
np.save("/Users/blakechellew/Documents/DustProject/taos_small.npy", taos_small)
'''

'''
#compare two sets of alphas
blake_alphas = np.load('blake_alphas.npy')
jiaoyue_alphas = np.load('jiaoyue_alphas.npy')
print("blake")
print(blake_alphas)
print('jiaoyue')
print(jiaoyue_alphas)

max_frac_diff = 0
for i in range(len(blake_alphas)):
    b = blake_alphas[i]
    j = jiaoyue_alphas[i]
    diff = np.abs(b-j)
    frac_diff = diff/b
    if frac_diff > max_frac_diff:
        max_frac_diff = frac_diff
print("max frac diff: ", max_frac_diff)
'''

plt.scatter(np.arange(100), np.arange(100), s=4, c='blue')
plt.show()


    


