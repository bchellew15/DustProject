from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from time import time

t1 = time()

#load data
fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
l = np.array(fiberinfo[:, 2])     # Galactic longitude, degrees
for i in range(len(l)):
    if l[i] > 180:
        l[i] = l[i] - 360
b = np.array(fiberinfo[:, 3])     # Galactic latitude, degrees
i100 = np.array(fiberinfo[:, 4])  # 100 micron intensity (MJy/sr)
hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
plate = np.array(hdulist[0].data)
wavelength = np.array(hdulist[1].data)  # Angstroms
flambda = np.array(hdulist[2].data)     # df/dlam, units of 1e-17 erg/s/cm^2/A
ivar = np.array(hdulist[3].data)        # inverse variance, units of 1/flambda^2
ivar *= (ivar > 0) #correct the negative values

t2 = time()
print("check 1: ", t2-t1)

print("flambda shape: ", flambda.shape)
print("wavelength shape: ", wavelength.shape)


#mask the high-intensity values:
mask_indices = np.arange(len(i100))[i100>10]
l = np.delete(l, mask_indices)
b = np.delete(b, mask_indices)
i100 = np.delete(i100, mask_indices)
plate = np.delete(plate, mask_indices)
flambda = np.delete(flambda, mask_indices, 0)
ivar = np.delete(ivar, mask_indices, 0)


t3 = time()
print("check 2: ", t3-t2)

def plot_alphas(i100, plate, flambda, ivar, color, bin=False):

    plate_ = np.copy(plate)
    ivar_ = np.copy(ivar)

    t3a = time()
    print("check 3a: ", t3a-t3)
    
    #calculate y
    y = np.copy(flambda)
    for i in range(0, len(wavelength)):
        y[:,i] *= wavelength[i]

    #y = np.multiply(flambda, wavelength)

    t3b = time()
    print("check 3b: ", t3b-t3a)
        
    #calculate x
    lam100 = 100 * pow(10,-6)
    c = 2.998 * pow(10, 8)
    freq100 = c / lam100
    x1 = np.copy(i100)
    x1 *= freq100
    #avg x1 over plates (assuming in ascending order)
    x2 = np.zeros(x1.shape) #will be array of averages

    t4 = time()
    print("check 3c: ", t4-t3b)
    
    boundaries = np.sort(np.unique(plate_, return_index=True)[1])
    for idx, b in enumerate(boundaries[:-1]):
        avgs = np.mean(x1[boundaries[idx]:boundaries[idx+1]], axis=0) #mean across plates, not wavelength
        x2[boundaries[idx]:boundaries[idx+1]] = avgs
        
    x = np.subtract(x1, x2)

    t5 = time()
    print("check 4: ", t5-t4)

    #unit conversion (see notes)
    x *= (1.617 * 10**-10)
    x = x.reshape(len(x), 1)

    #unit conversion for ivar
    ivar_ /= np.power(wavelength, 2)

    t5b = time()
    print("check 4b: ", t5b-t5)
        
    #calculate alpha
    xx = np.multiply(x, x)
    yx = np.multiply(y, x)
    yxsig = np.multiply(yx, ivar_)
    xxsig = np.multiply(xx, ivar_)
    sums1 = np.sum(yxsig, axis=0)
    sums2 = np.sum(xxsig, axis=0)
    alphas = np.divide(sums1, sums2)
    alpha_std = np.sqrt(1/sums2)
        
    t6 = time()
    print("check 5: ", t6-t5b)
        
    if bin:
        #binning:
        binned_lambdas = np.arange(3900, 9000, 50)
        binned_alphas = np.zeros(binned_lambdas.shape)
        binned_std = np.zeros(binned_lambdas.shape)
        for i, lmda in enumerate(binned_lambdas):
            indices = np.where((wavelength > lmda-25) & (wavelength < lmda+25))[0]
            relevant_alphas = alphas[indices]
            relevant_stds = alpha_std[indices]
            #weighted average:
            variance = np.power(relevant_stds, 2)
            numerator = np.sum(np.divide(relevant_alphas, variance))
            denominator = np.sum(np.divide(1, variance))
            avg1 = numerator / denominator
            avg2 = 1 / denominator
            binned_alphas[i] = avg1
            binned_std[i] = np.sqrt(avg2)
        plt.plot(binned_lambdas[0:-2], binned_alphas[0:-2], color)
        plt.plot(binned_lambdas[0:-2], binned_std[0:-2], color)
        
    else: 
        #plot alpha vs. wavelength
        plt.plot(wavelength[50:-100], alphas[50:-100], color)
        plt.plot(wavelength[50:-100], alpha_std[50:-100], color)

    t7 = time()
    print("check 6: ", t7-t6)
        
    print("alphas :", alpha_std)

    plt.xlabel("Wavelength")
    plt.ylabel("Alpha")


#figure 3: 
plot_alphas(i100, plate, flambda, ivar, 'k', bin=True)
plt.show()


'''
#figures 4a and 6a: only certain ranges of b.
mask_indices_b1 = np.arange(len(i100))[np.abs(b) > 35]
i100_b1 = np.delete(i100, mask_indices_b1)
plate_b1 = np.delete(plate, mask_indices_b1)
flambda_b1 = np.delete(flambda, mask_indices_b1, 0)
ivar_b1 = np.delete(ivar, mask_indices_b1, 0)

mask_indices_b2 = np.arange(len(i100))[(np.abs(b) < 35) | (np.abs(b) > 50)]
i100_b2 = np.delete(i100, mask_indices_b2)
plate_b2 = np.delete(plate, mask_indices_b2)
flambda_b2 = np.delete(flambda, mask_indices_b2, 0)
ivar_b2 = np.delete(ivar, mask_indices_b2, 0)

mask_indices_b3 = np.arange(len(i100))[np.abs(b) < 50]
i100_b3 = np.delete(i100, mask_indices_b3)
plate_b3 = np.delete(plate, mask_indices_b3)
flambda_b3 = np.delete(flambda, mask_indices_b3, 0)
ivar_b3 = np.delete(ivar, mask_indices_b3, 0)


#figure 4a:
plot_alphas(i100_b1, plate_b1, flambda_b1, ivar_b1, 'b', bin=True)
plot_alphas(i100_b2, plate_b2, flambda_b2, ivar_b2, 'r', bin=True)
plot_alphas(i100_b3, plate_b3, flambda_b3, ivar_b3, 'g', bin=True)
plot_alphas(i100, plate, flambda, ivar, 'k', bin=True)
plt.show()



#figure 6a:
plot_alphas(i100_b1, plate_b1, flambda_b1, ivar_b1, 'b')
plot_alphas(i100_b2, plate_b2, flambda_b2, ivar_b2, 'r')
plot_alphas(i100_b3, plate_b3, flambda_b3, ivar_b3, 'g')
plot_alphas(i100, plate, flambda, ivar, 'k')
plt.xlim(6530, 6770)
plt.ylim(0, 1.4)
plt.show()
'''
'''
#figures 4b and 6b: certain ranges of l

mask_indices_l1 = np.arange(len(i100))[np.abs(l) > 60]
i100_l1 = np.delete(i100, mask_indices_l1)
plate_l1 = np.delete(plate, mask_indices_l1)
flambda_l1 = np.delete(flambda, mask_indices_l1, 0)
ivar_l1 = np.delete(ivar, mask_indices_l1, 0)

mask_indices_l2 = np.arange(len(i100))[(np.abs(l) < 60) | (np.abs(l) > 120)]
i100_l2 = np.delete(i100, mask_indices_l2)
plate_l2 = np.delete(plate, mask_indices_l2)
flambda_l2 = np.delete(flambda, mask_indices_l2, 0)
ivar_l2 = np.delete(ivar, mask_indices_l2, 0)

mask_indices_l3 = np.arange(len(i100))[np.abs(l) < 120]
i100_l3 = np.delete(i100, mask_indices_l3)
plate_l3 = np.delete(plate, mask_indices_l3)
flambda_l3 = np.delete(flambda, mask_indices_l3, 0)
ivar_l3 = np.delete(ivar, mask_indices_l3, 0)

#figure 4b:
plot_alphas(i100_l1, plate_l1, flambda_l1, ivar_l1, 'b', bin=True)
plot_alphas(i100_l2, plate_l2, flambda_l2, ivar_l2, 'r', bin=True)
plot_alphas(i100_l3, plate_l3, flambda_l3, ivar_l3, 'g', bin=True)
plot_alphas(i100, plate, flambda, ivar, 'k', bin=True)
plt.show()

#figure 6b:
plot_alphas(i100_l1, plate_l1, flambda_l1, ivar_l1, 'b')
plot_alphas(i100_l2, plate_l2, flambda_l2, ivar_l2, 'r')
plot_alphas(i100_l3, plate_l3, flambda_l3, ivar_l3, 'g')
plot_alphas(i100, plate, flambda, ivar, 'k')
plt.xlim(6530, 6770)
plt.ylim(0, 1.4)
plt.show()
'''
'''
#figure 5:
plot_alphas(i100, plate, flambda, ivar, 'k')
plt.xlim(4830, 5040)
plt.ylim(0, 1)
plt.show()

plot_alphas(i100, plate, flambda, ivar, 'k')
plt.xlim(6530, 6770)
plt.ylim(0, 1)
plt.show()
'''


#no longer needed:
'''
def plot_alphas(i100, plate, flambda, ivar, color):

    #calculate y
    y = np.array(flambda)
    for i in range(0, len(wavelength)):
        y[:,i] *= wavelength[i]
    
    #calculate x
    lam100 = 100 * pow(10,-6)
    c = 2.998 * pow(10, 8)
    freq100 = c / lam100
    x1 = np.array(i100)
    x1 *= freq100
    #avg x1 over plates (assuming in ascending order)
    x2 = np.zeros(x1.shape) #will be array of averages
    plate_numbers = np.unique(plate)
    for p in plate_numbers:
        mask = np.isin(plate, p)
        vals_on_plate = x1[mask]
        avg = np.mean(vals_on_plate)
        np.putmask(x2, mask, avg)
    x = np.subtract(x1, x2)
    x = x.reshape(len(x), 1)
    #unit conversion (see notes)
    x *= (1.617 * 10**-10)

    #calculate alpha
    binned_lambdas = np.arange(3800, 9000, 50)
    alphas = np.zeros(binned_lambdas.shape)
    alpha_varianc = np.zeros(binned_lambdas.shape)
    for i, lmda in enumerate(binned_lambdas):
        indices = np.where((wavelength > lmda) & (wavelength < lmda+50))[0]
        ylambda = y[:,indices]
        varlambda = ivar[:,indices]
        yx = np.multiply(ylambda, x)
        yxsig = np.multiply(yx, varlambda)
        sum1 = np.sum(yxsig)
        xx = np.multiply(x, x)
        xxsig = np.multiply(xx, varlambda)
        sum2 = np.sum(xxsig)
        alpha = np.divide(sum1, sum2)
        alphas[i] = alpha
        alpha_variance[i] = 1/sum2
    #plot alpha vs. wavelength
    plt.plot(binned_lambdas[1:-2], alphas[1:-2], color)
    plt.plot(binned_lambdas[1:-2], alpha_variance[1:-2], color)
    plt.xlabel("Wavelength")
    plt.ylabel("Alpha")
'''
