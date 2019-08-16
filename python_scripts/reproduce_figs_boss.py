from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from time import time
import sys

#reproduce figures 3-6 from the paper
#see reproduce_figs.py
#this is a modification that aims to use less RAM

#command line args
if len(sys.argv) == 1:
    print("Usage: reproduce_figs.py mode[1d, 2d, iris, iris_1d] boss[0, 1] save[save]") 
else:
    mode = sys.argv[1]

if len(sys.argv) >= 3:
    boss = int(sys.argv[2])
else:
    boss = False

if len(sys.argv) == 4:
    save = True
else:
    save = False
    

#calculate x, y, alpha. then plot alpha vs. wavelength
def plot_alphas(i100, plate, flambda, ivar, color, bin=False):

    #calculate y

    y = np.memmap('y_mem.dat', dtype=np.float32, mode='w+', shape=flambda.shape)
    for i in range(flambda.shape[0]):
        y[i]  = np.multiply(flambda[i], wavelength, dtype='float32')

    print("check 5")
        
    #calculate x
    lam100 = 100 * pow(10,-6) #100 microns
    c = 2.998 * pow(10, 8) #speed of light
    freq100 = c / lam100

    #copy i100:
    x1 = np.memmap('x1_mem.dat', dtype=np.float32, mode='w+', shape=i100.shape)
    x1[:] = i100[:]
    
    print("check 6")

    for i in range(x1.shape[0]):
        x1[i] *= freq100
    #avg x1 over plates (assuming in ascending order)
    x2 = np.memmap('x2_mem.dat', dtype=np.float32, mode='w+', shape=x1.shape) #x2 will be array of averages
    boundaries = np.sort(np.unique(plate, return_index=True)[1])
    for idx, b in enumerate(boundaries[:-1]):
        avgs = np.mean(x1[boundaries[idx]:boundaries[idx+1]], axis=0) #mean across plates, not wavelength
        x2[boundaries[idx]:boundaries[idx+1]] = avgs
    #last section:
    avgs = np.mean(x1[boundaries[-1]:], axis=0) #mean across plates, not wavelength
    x2[boundaries[-1]:] = avgs

    print("check 7")
    
    x = np.memmap('x_mem.dat', dtype=np.float32, mode='w+', shape=x1.shape)
    for i in range(x.shape[0]):
        x[i] = np.subtract(x1[i], x2[i])

    print("check 8")

    #x unit conversion
    if boss:
        unit_factor = 7.392 * 10**-11
    else:
        unit_factor = 1.617 * 10**-10
    for i in range(x.shape[0]):
        x[i] *= unit_factor
    if mode == '1d' or mode == 'iris_1d':
        x = x.reshape(len(x), 1)
        
    #calculate alpha
    print("calculating alphas")
    
    xx = np.memmap('xx_mem.dat', dtype=np.float32, mode='w+', shape=x.shape)
    for i in range(x.shape[0]):
        xx[i] = np.multiply(x[i], x[i])

    print("check 9")

    yx = np.memmap('yx_mem.dat', dtype=np.float32, mode='w+', shape=y.shape)
    for i in range(y.shape[0]):
        yx[i] = np.multiply(y[i], x[i])

    print("check 10")
    
    yxsig = np.memmap('yxsig_mem.dat', dtype=np.float32, mode='w+', shape=yx.shape)
    for i in range(yx.shape[0]):
        yxsig[i] = np.multiply(yx[i], ivar[i])

    print("check 11")
    
    xxsig = np.memmap('xxsig_mem.dat', dtype=np.float32, mode='w+', shape=y.shape)
    for i in range(xx.shape[0]):
        xxsig[i] = np.multiply(xx[i], ivar[i])

    print("check 12")

    sums1 = np.zeros(yxsig.shape[1])
    for i in range(yxsig.shape[1]):
        sums1[i] = np.sum(yxsig[:,i])
    sums2 = np.zeros(xxsig.shape[1])
    for i in range(xxsig.shape[1]):
        sums2[i] = np.sum(xxsig[:,i])
    alphas = np.divide(sums1, sums2)
    alpha_std = np.sqrt(1/sums2)

    print("check 13")

    print("Alphas:") #TEMP
    print(alphas) #TEMP
    print(alpha_std) #TEMP
    
    #save alphas:
    if save:
        np.save('../alphas_and_stds/alphas_1d_boss.npy', alphas)
        np.save('../alphas_and_stds/alphas_1d_boss_stds.npy', alpha_std)
        print("alphas saved")

    #MAKE THIS INTO ITS OWN FUNCTION
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
        #plot alpha vs. wavelength
        plt.scatter(binned_lambdas, binned_alphas, s=5, c=color)
        plt.scatter(binned_lambdas, binned_std, s=5, c='r') #change this back to c=color
    else: 
        #plot alpha vs. wavelength
        plt.scatter(wavelength[50:-100], alphas[50:-100], s=5, c=color)
        plt.scatter(wavelength[50:-100], alpha_std[50:-100], s=5, c='r') #change back to c=color

    plt.xlabel("Wavelength")
    plt.ylabel("Alpha")

    return alphas, alpha_std, wavelength


#load data for BOSS
if boss:

    filenames = ['skyfibers_lam0.fits', 'skyfibers_lam1.fits', 'skyfibers_lam2.fits', 'skyfibers_lam3.fits', 'skyfibers_lam4.fits', \
                 'skyfibers_lam5.fits', 'skyfibers_lam6.fits', 'skyfibers_lam7.fits', 'skyfibers_lam8.fits', 'skyfibers_lam9.fits']

    #create unique plate identifier
    hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[0])
    plate = hdulist[6].data
    mjd = hdulist[5].data
    plate = 10000*plate + mjd%10000 #plate is now a unique identifier
    #RA and dec are hdulist[3] and hdulist[4]

    #get i100
    i100_old = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/SFD_i100_at_BOSS_locations.txt")[:,2].astype('float32')
    if mode == '2d':
        i100 = np.load("/Volumes/TOSHIBA EXT/Dust_Overflow/i100_tao_boss.npy", mmap_mode='r') #type: float64
    elif mode == 'iris':
        i100 = np.load("/Volumes/TOSHIBA EXT/Dust_Overflow/i100_tao_boss_iris.npy", mmap_mode='r').astype('float32')
    elif mode == 'iris_1d':
        i100 = np.load("/Users/blakechellew/Documents/DustProject/IRIS/iris_i100_at_boss.npy", mmap_mode='r').astype('float32')
    elif mode == '1d':
        i100 = i100_old
    else:
        print("Error: invalid mode")

#load data for SDSS
else:
    #load data
    fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
    #l = np.array(fiberinfo[:, 2])     # Galactic longitude, degrees
    #for i in range(len(l)): #convert l to range: -180 to 180
    #    if l[i] > 180:
    #        l[i] = l[i] - 360
    #b = np.array(fiberinfo[:, 3])     # Galactic latitude, degrees

    #get i100
    i100_old = np.array(fiberinfo[:, 4]).astype('float32')  # 100 micron intensity (MJy/sr)
    if mode == '2d':
        i100 = np.load("/Volumes/TOSHIBA EXT/Dust_Overflow/i100_tao.npy", mmap_mode='r').astype('float32')  #i100_tao.npy
    elif mode == 'iris':
        i100 = np.load("/Volumes/TOSHIBA EXT/Dust_Overflow/i100_iris_tao.npy", mmap_mode='r').astype('float32')
    elif mode == 'iris_1d':
        i100 = np.load("/Users/blakechellew/Documents/DustProject/npy_files/iris_i100_at_sfd.npy", mmap_mode='r').astype('float32')
    elif mode == '1d':
        i100 = i100_old
    else:
        print("Error: invalid mode")

    #get flambda and ivar
    hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
    plate = np.array(hdulist[0].data)
    wavelength = np.array(hdulist[1].data)  # Angstroms

    #create memmaps of flambda and ivar:
    flambda = np.memmap('flambda_mem.dat', dtype=np.float32, mode='w+', shape=hdulist[2].data.shape) # df/dlam, units of 1e-17 erg/s/cm^2/A
    flambda[:] = hdulist[2].data[:]
    ivar = np.memmap('ivar_mem.dat', dtype=np.float32, mode='w+', shape=hdulist[3].data.shape) # inverse variance, units of 1/flambda^2
    ivar[:] = hdulist[3].data[:]

#plot_alphas(i100, plate, flambda, ivar, 'k', bin=True)
#plt.show()



alphas_10 = np.zeros(0)
alpha_std_10 = np.zeros(0)
wavelength_10 = np.zeros(0)

if boss:
    num_files = 10
else:
    num_files = 1


for i in range(num_files):
    hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[i])
    flambda = hdulist[0].data #type: float64
    ivar = hdulist[1].data #type: float64
    wavelength = 10**( min(hdulist[2].data) + np.arange(flambda[0].shape[0])*1e-4 ).astype('float32')

    #process ivar:
    for i in range(ivar.shape[0]):
        ivar[i] *= (ivar[i] > 0) #correct the negative values
    #new masking
    for i in range(ivar.shape[0]):
        if i100_old[i] > 10:
            ivar[i] = 0
    #convert ivar to ivar of y
    for i in range(ivar.shape[0]):
        ivar[i] /= np.power(wavelength, 2)
    print("loaded data")

    alphas_i, alpha_std_i, wavelength_i = plot_alphas(i100, plate, flambda, ivar, 'k', bin=False)

    plt.clf() #TEST
    
    alphas_10 = np.append(alphas_10, alphas_i)
    
    alpha_std_10 = np.append(alpha_std_10, alpha_std_i)
    wavelength_10 = np.append(wavelength_10, wavelength_i)

np.save('../alphas_and_stds/alphas_1d_boss.npy', alphas_10)
np.save('../alphas_and_stds/alphas_1d_boss_stds.npy', alpha_std_10)
np.save('../alphas_and_stds/wavelength_boss_stds.npy', wavelength_10)
print("alphas saved")
    
bin = True
    
if bin:
    #binning:
    binned_lambdas = np.arange(3900, 9000, 50)
    binned_alphas = np.zeros(binned_lambdas.shape)
    binned_std = np.zeros(binned_lambdas.shape)
    for i, lmda in enumerate(binned_lambdas):
        indices = np.where((wavelength_10 > lmda-25) & (wavelength_10 < lmda+25))[0]
        relevant_alphas = alphas_10[indices]
        relevant_stds = alpha_std_10[indices]
        #weighted average:
        variance = np.power(relevant_stds, 2)
        numerator = np.sum(np.divide(relevant_alphas, variance))
        denominator = np.sum(np.divide(1, variance))
        avg1 = numerator / denominator
        avg2 = 1 / denominator
        binned_alphas[i] = avg1
        binned_std[i] = np.sqrt(avg2)
    #plot alpha vs. wavelength
    plt.scatter(binned_lambdas, binned_alphas, s=5, c='k')
    plt.scatter(binned_lambdas, binned_std, s=5, c='r')
else: 
    #plot alpha vs. wavelength
    plt.scatter(wavelength[50:-100], alphas[50:-100], s=5, c='k')
    plt.scatter(wavelength[50:-100], alpha_std[50:-100], s=5, c='r')

plt.show()



