#similar to reproduce_figs.py and reproduce_figs_boss.py
#but find alpha at just one wavelength (avoid waiting such a long time)

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from time import time
import sys
from math import floor #for calculating bin ranges

#command line args
if len(sys.argv) != 4:
    print("Usage: reproduce_figs.py [mode: 1d, 2d, iris, iris_1d] [boss: 0, 1] [save: 0, 1]")
    exit(0)

mode = sys.argv[1]
boss = int(sys.argv[2])
save = int(sys.argv[3])
    
#create scatterplot
def plot_alphas(alphas, alpha_std, wavelength, color1='k', color2='r', bin=False):

    if bin:
        #calculate bin ranges
        lambda_range = wavelength[-1] - wavelength[0]
        left_over = lambda_range - 50*floor(lambda_range / 50)
        
        binned_lambdas = np.arange(wavelength[0]+left_over/2, wavelength[-1], 50)
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
        plt.plot(binned_lambdas, binned_alphas, c=color1, drawstyle='steps')
        plt.plot(binned_lambdas, binned_std, c=color2, drawstyle='steps')
    else: 
        #plot alpha vs. wavelength
        plt.plot(wavelength[50:-100], alphas[50:-100], c=color1, drawstyle='steps')
        plt.plot(wavelength[50:-100], alpha_std[50:-100], c=color2, drawstyle='steps')

    plt.xlabel("Wavelength")
    plt.ylabel("Alpha")
    
#calculate x, y, alpha. then plot alpha vs. wavelength
def calc_alphas(i100, plate, flambda, ivar):

    #calculate y

    y = np.memmap('y_mem.dat', dtype=np.float32, mode='w+', shape=flambda.shape)
    for i in range(flambda.shape[0]):
        y[i]  = np.multiply(flambda[i], wavelength, dtype='float32')


    print("biggest and smallest y:")
    print(y[np.argsort(y)[:10]])
    print(y[np.argsort(y)[-10:]])
    print("ivar of those indices:")
    print(ivar[np.argsort(y)[:10]])
    print(ivar[np.argsort(y)[-10:]])
    


    plt.hist(i100, bins=50, range=(0,10))
    plt.show()
    #plt.hist(flambda, bins=50, range=(-10, 10))
    #plt.show()   
    #plt.hist(y, bins=50, range=(-50000, 50000))
    #plt.show()
    #plt.hist(ivar, bins=50)
    #plt.show()

    print("i100:", i100[3178])
    print("flam:", flambda[3178])
    print("ivar:", ivar[3178])
    print("y:", y[3178])

    exit()
        
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

    print("x1:", x1[3178])
    print("x2:", x2[3178])

    print("all x1:", x1)
    print("num larger x2:", len(np.unique(x2[x2>x2[3187]])))
    print("num smaller x2:", len(np.unique(x2[x2<x2[3187]])))

    #find std dev of elements on plate:
    x_on_plate = x1[plate == plate[3178]]
    print("all x1 on plate:", x_on_plate)
    idx_on_plate = np.argsort(x_on_plate)
    largest_on_plate = x_on_plate[idx_on_plate[-10:]]
    smallest_on_plate = x_on_plate[idx_on_plate[:10]]
    print("10 largest x1s:", largest_on_plate)
    print("10 smallest x1s:", smallest_on_plate)
    print("std on plate:", np.std(x_on_plate))
    

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

    print("x:", x[3178])

    plt.hist(x, bins=50, range=(-1000, 1000))
    plt.show()
    exit(0)
        
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

    #TEMP, uncomment this block, delete the one beneath
    '''
    sums1 = np.zeros(yxsig.shape[1]) 
    for i in range(yxsig.shape[1]):
        sums1[i] = np.sum(yxsig[:,i])
    sums2 = np.zeros(xxsig.shape[1])
    for i in range(xxsig.shape[1]):
        sums2[i] = np.sum(xxsig[:,i])
    alphas = np.divide(sums1, sums2)
    alpha_std = np.sqrt(1/sums2)
    '''

    print("yxsig:", yxsig)
    print("xxsig:", xxsig)
    print("outside of range:", yxsig[(yxsig<-.1) * (yxsig>.1)])
    print("outside of range:", xxsig[(xxsig<0) * (xxsig>.008)])

    #histograms of xxsig and yxsig
    '''
    plt.hist(yxsig, bins=50, range=(-.1, .1))
    plt.title("yxsig")
    plt.show()

    plt.hist(xxsig, bins=50, range=(0, .008))
    plt.title("xxsig")
    plt.show()
    #see if some plates are contributing more:
    #for each unique plates, get all the lambdas, see if any don't match
    unique_plate_idces = np.unique(plate, return_index=True)[1] #2512 unique plates

    yxsig_plate_means = []
    xxsig_plate_means = []
    for i in unique_plate_idces:
        p = hdulist[6].data[i]
        d = mjd[i]
        idces =  np.nonzero((hdulist[6].data == p) * (mjd == d))[0] #indices of rows with this p and d

        #yxsig and xxsig on this plate:
        yxsig_sub = yxsig[idces]
        yxsig_mean = np.mean(yxsig_sub)
        yxsig_plate_means.append(yxsig_mean)
        xxsig_sub = xxsig[idces]
        xxsig_mean = np.mean(xxsig_sub)
        xxsig_plate_means.append(xxsig_mean)
    plt.hist(yxsig_plate_means, bins=50, range=(-.1, .1))
    plt.title("yxsig plate means")
    plt.show()
    plt.hist(xxsig_plate_means, bins=50, range=(0, .002))
    plt.title("xxsig plate means")
    plt.show()
    ''' 
    
    #save:
    np.save("xxsig_iris_1d_82119.npy", xxsig)
    np.save("yxsig_iris_1d_82119.npy", yxsig)
    
    sums1 = np.sum(yxsig)
    sums2 = np.sum(xxsig)
    alphas = np.divide(sums1, sums2)
    alpha_std = np.sqrt(1/sums2)

    print("sums 1:", sums1)
    print("sums 2:", sums2)

    print("check 13")

    return alphas, alpha_std, wavelength


#load data for BOSS
if boss:
    
    filenames = ['skyfibers_lam0.fits', 'skyfibers_lam1.fits', 'skyfibers_lam2.fits', 'skyfibers_lam3.fits', 'skyfibers_lam4.fits', \
                 'skyfibers_lam5.fits', 'skyfibers_lam6.fits', 'skyfibers_lam7.fits', 'skyfibers_lam8.fits', 'skyfibers_lam9.fits']
    num_files = len(filenames)

    #get dimensions for selecting i100s:
    num_columns = np.zeros(num_files)
    for i in range(num_files):
        hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[i])
        num_columns[i] = hdulist[0].data.shape[1]
    #elements of num_columns will be tuples: (start_idx, num_elements)
    num_columns = [(np.sum(num_columns[:i], dtype=np.int32), int(n)) for (i, n) in enumerate(num_columns)]

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
    num_files = 1
    
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

#compute alphas separately for each file, then combine
alphas_10 = np.zeros(0)
alpha_std_10 = np.zeros(0)
wavelength_10 = np.zeros(0)

num_files = 1 #TEMP, remove

#preprocessing and calculate alphas for each file
for j in range(num_files):
    if boss:
        hdulist = fits.open("/Volumes/TOSHIBA EXT/Dust_Overflow/" + filenames[5]) #TEMP (5 -> j)
        flambda = hdulist[0].data[:, 149] #TEMP: remove last brackets  #type: float64
        ivar = hdulist[1].data[:,149] #TEMP: remove last brackets   #type: float64
        #wavelength = 10**( min(hdulist[2].data) + np.arange(flambda[0].shape[0])*1e-4 ).astype('float32')
        wavelength = 6302.31347656 #TEMP, remove and uncomment prev line

    #process ivar:
    for i in range(ivar.shape[0]):
        ivar[i] *= (ivar[i] > 0) #correct the negative values
    #new masking (threshold 10)
    for i in range(ivar.shape[0]):
        if i100_old[i] > 10: #10
            ivar[i] = 0

    ivar[3178] = 0 #TEMP
        #if plate[i] == 72576658 or plate[i] == 72586605 or plate[i] == 72596603 \
        #   or plate[i] == 72606679 or plate[i] == 72626659 or plate[i] == 72626683:
        #    ivar[i] = 0

    #convert ivar to ivar of y
    for i in range(ivar.shape[0]):
        ivar[i] /= np.power(wavelength, 2)
    print("loaded data")

    #get subset of i100:
    if boss and (mode == '2d' or mode == 'iris'):
        start_idx =  num_columns[j][0]
        num_elements = num_columns[j][1]
        i100_sub = i100[:, start_idx:start_idx+num_elements]
    else:
        i100_sub = i100

    #calculate alphas
    alphas_i, alpha_std_i, wavelength_i = calc_alphas(i100_sub, plate, flambda, ivar)
    alphas_10 = np.append(alphas_10, alphas_i)
    alpha_std_10 = np.append(alpha_std_10, alpha_std_i)
    wavelength_10 = np.append(wavelength_10, wavelength_i)

if save:
    np.save('../alphas_and_stds/alphas_boss_iris_1d_82019_5.npy', alphas_10)
    np.save('../alphas_and_stds/alphas_boss_iris_1d_stds_82019_5.npy', alpha_std_10)
    #np.save('../alphas_and_stds/wavelength_boss.npy', wavelength_10)
    print("alphas saved")

#plot
#plot_alphas(alphas_10, alpha_std_10, wavelength_10, bin=True) #TEMP: uncomment these lines
#plt.show()

print("alpha:", alphas_10) #TEMP: remove these lines
print("std:", alpha_std_10)



