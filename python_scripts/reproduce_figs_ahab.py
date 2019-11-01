#see reproduce_figs.py and reproduce_figs_boss.py
#this version uses a lot more RAM and is made to work on AHAB

#for bootstrapping:
# 


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from time import time
from time import sleep
import sys
from math import floor #for calculating bin ranges

#for tracking RAM:
import os
import psutil
import threading
import _thread #older version

#run continuously to make sure the program is not using too much RAM:
def check_memory():
    while True:
        process = psutil.Process(os.getpid())
        mem_usage = process.memory_info().rss / 1e9
        print("memory: ", mem_usage)  # in bytes
        if mem_usage > 28:
            print("exiting due to excessive RAM usage")
            _thread.interrupt_main()   #threading.currentThread()
        sleep(4)
        
thr1 = threading.Thread(target = check_memory, daemon=True)
thr1.start()

#command line args
if len(sys.argv) < 4 or len(sys.argv) > 6:
    print("Usage: reproduce_figs.py [mode: 1d, 2d, iris, iris_1d] [boss: 0, 1] [save: 0, savekey] [threshold=10] [bootstrap=0]")
    exit(0)
if len(sys.argv) >= 5:
    threshold = float(sys.argv[4]) #masking threshold
else:
    threshold = 10 #default
if len(sys.argv) == 6:
    bootstrap = int(sys.argv[5])
else:
    bootstrap = 0 #default

mode = sys.argv[1]
boss = int(sys.argv[2])
save = sys.argv[3]
    
#create scatterplot
def plot_alphas(alphas, alpha_std, wavelength, color1='k', color2='r', bin=False):

    if bin:
        #calculate bin ranges
        lambda_range = wavelength[-1] - wavelength[0]
        left_over = lambda_range - 50*floor(lambda_range / 50)
        binned_lambdas = np.arange(wavelength[0]+left_over/2, wavelength[-2], 50)

        binned_alphas = np.zeros(binned_lambdas.shape)
        binned_std = np.zeros(binned_lambdas.shape)
        for i, lmda in enumerate(binned_lambdas):
            indices = np.where((wavelength > lmda) & (wavelength < lmda+50))[0]
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
    
#calculate x, y, alpha, then plot alpha vs. wavelength
def calc_alphas(i100, plate, flambda, ivar):

    print("calculating x and y")
    
    #calculate y
    y = np.multiply(flambda, wavelength, dtype='float32')    
        
    #100 micron frequency
    lam100 = 100 * pow(10,-6) #100 microns
    c = 2.998 * pow(10, 8) #speed of light
    freq100 = c / lam100

    #copy i100:
    x1 = np.copy(i100)    
    
    #calculate x1 and x2
    x1 *= freq100   
    
    #avg x1 over plates (assuming grouped together)
    x2 = np.zeros(x1.shape, dtype=np.float32)    
    boundaries = np.sort(np.unique(plate, return_index=True)[1])
    for idx, b in enumerate(boundaries[:-1]):
        avgs = np.mean(x1[boundaries[idx]:boundaries[idx+1]], axis=0) #mean across plates, not wavelength
        x2[boundaries[idx]:boundaries[idx+1]] = avgs
    #last section:
    avgs = np.mean(x1[boundaries[-1]:], axis=0) #mean across plates, not wavelength
    x2[boundaries[-1]:] = avgs

    #calculate x
    x = np.subtract(x1, x2)

    #x unit conversion
    if boss:
        unit_factor = 7.384 * 10**-11
    else:
        unit_factor = 1.617 * 10**-10
    x *= unit_factor
    if mode == '1d' or mode == 'iris_1d':
        x = x.reshape(len(x), 1)
        
    print("calculating alphas")

    #calculate alpha
    xx = np.multiply(x, x)
    yx = np.multiply(y, x)
    yxsig = np.multiply(yx, ivar)
    xxsig = np.multiply(xx, ivar)    
    sums1 = np.sum(yxsig, axis=0)
    sums2 = np.sum(xxsig, axis=0)

    if bootstrap:
        return np.sum(yxsig, axis=0), np.sum(xxsig, axis=0)

    #check for division by 0
    for i in range(len(sums1)):
        if sums1[i] == 0 and sums2[i] == 0:
            sums1[i] = np.nan
            sums2[i] = np.nan
            print("sums1 and sums2 are both 0 at this wavelength")
        
    alphas = np.divide(sums1, sums2)
    alpha_std = np.sqrt(1/sums2)

    print("finished calculating alphas")
    return alphas, alpha_std, wavelength
    


#load data for BOSS
if boss:
    
    filenames = ['skyfibers_lam0.fits', 'skyfibers_lam1.fits', 'skyfibers_lam2.fits', 'skyfibers_lam3.fits', 'skyfibers_lam4.fits', \
                 'skyfibers_lam5.fits', 'skyfibers_lam6.fits', 'skyfibers_lam7.fits', 'skyfibers_lam8.fits', 'skyfibers_lam9.fits']
    num_files = len(filenames)

    #get dimensions for selecting i100s:
    num_columns = np.zeros(num_files)
    for i in range(num_files):
        hdulist = fits.open("../data/" + filenames[i])
        num_columns[i] = hdulist[0].data.shape[1]
    #elements of num_columns will be tuples: (start_idx, num_elements)
    num_columns = [(np.sum(num_columns[:i], dtype=np.int32), int(n)) for (i, n) in enumerate(num_columns)]

    #create unique plate identifier
    hdulist = fits.open("../data/" + filenames[0])
    plate = hdulist[6].data
    mjd = hdulist[5].data
    plate = 10000*plate + mjd%10000 #plate is now a unique identifier

    #get i100 (type: float64)
    i100_old = np.loadtxt("../data/SFD_i100_at_BOSS_locations.txt")[:,2]
    if mode == '2d':
        i100 = np.load("../data/i100_tao_boss.npy")
    elif mode == 'iris':
        i100 = np.load("../data/i100_tao_boss_iris.npy")
    elif mode == 'iris_1d':
        i100 = np.load("../data/iris_i100_at_boss.npy") 
    elif mode == '1d':
        i100 = i100_old
    else:
        print("Error: invalid mode")

#load data for SDSS
else:
    num_files = 1
    
    #load data
    fiberinfo = np.loadtxt('../data/fiberinfo_halpha.dat')
    #l = np.array(fiberinfo[:, 2])     # Galactic longitude, degrees
    #for i in range(len(l)): #convert l to range: -180 to 180
    #    if l[i] > 180:
    #        l[i] = l[i] - 360
    #b = np.array(fiberinfo[:, 3])     # Galactic latitude, degrees

    #get i100 (type: float64)
    i100_old = np.array(fiberinfo[:, 4])# 100 micron intensity (MJy/sr)
    if mode == '2d':
        i100 = np.load("../data/i100_tao.npy") #i100_tao.npy
    elif mode == 'iris':
        i100 = np.load("../data/i100_iris_tao.npy")
    elif mode == 'iris_1d':
        i100 = np.load("../data/iris_i100_at_sfd.npy")
    elif mode == '1d':
        i100 = i100_old
    else:
        print("Error: invalid mode")
        
    #get flambda and ivar
    hdulist = fits.open('../data/SDSS_allskyspec.fits')
    plate = np.array(hdulist[0].data)
    wavelength = np.array(hdulist[1].data)  # Angstroms

    #create memmaps of flambda and ivar:
    flambda = hdulist[2].data
    ivar = hdulist[3].data
    
#compute alphas separately for each file, then combine
alphas_10 = np.zeros(0)
alpha_std_10 = np.zeros(0)
wavelength_10 = np.zeros(0)

if bootstrap:
    unique_plates = np.unique(plate)
    print("number of unique plates:")

    #set up arrays for xxsig and yxsig. each unique plate has an entry at each wavelength
    yxsigs_bootstrap = np.zeros((len(unique_plates), 0))
    xxsigs_bootstrap = np.zeros((len(unique_plates), 0))

    #print("shapes before:")
    #print(i100_old.shape)
    #print(i100.shape)
    #print(plate.shape)

    #bootstrap_indices = np.random.choice([i for i in range(i100_old.shape[0]) if i != 3178], i100_old.shape[0])
    #bootstrap_indices = np.array([i for i in range(i100_old.shape[0]) if i != 3178])
    #bootstrap_indices = np.append(bootstrap_indices, bootstrap_indices[-1]) #to make it the right length                   
    #np.random.shuffle(bootstrap_indices)
    #plate_temp = plate[bootstrap_indices]
    #bootstrap_indices = bootstrap_indices[np.argsort(plate_temp)]
    #bootstrap_indices = np.array(range(i100_old.shape[0])) #test
    #bootstrap_indices[:10000] = bootstrap_indices[30000:40000]
    #print("bootstrap indices:", bootstrap_indices)

    #i100_old = i100_old[bootstrap_indices]
    #i100 = i100[bootstrap_indices]
    #plate = plate[bootstrap_indices]

    #print("shapes after:")
    #print(i100_old.shape)
    #print(i100.shape)
    #print(plate.shape)
    
#bootstrap code:
#btw stds prob irrelevant


#preprocessing and calculate alphas for each file
for j in range(num_files):
    if boss:
        hdulist = fits.open("../data/" + filenames[j])
        flambda = hdulist[0].data #type: float64
        ivar = hdulist[1].data #type: float64
        wavelength = 10**( min(hdulist[2].data) + np.arange(flambda[0].shape[0])*1e-4 ).astype('float32')

    #if bootstrap:
    #    print("shapes before:")
    #    print(flambda.shape)
    #    print(ivar.shape)
    #        
    #    flambda = flambda[bootstrap_indices]
    #    ivar = ivar[bootstrap_indices]
    #        
    #    print("shapes after:")
    #    print(flambda.shape)
    #    print(ivar.shape)

    #process ivar:
    ivar *= (ivar > 0)
    
    #masking
    ivar[i100_old > threshold] = 0
    #mask plates with large averages
    for p in np.unique(plate):
        if np.mean(i100_old[plate==p]) > threshold:
            ivar[plate==p] = 0
            print("masking whole plate")        
            
    if boss: # and not bootstrap: #TEMP might have to change this
        ivar[3178] = 0 #data at this location is bad
        
    #convert ivar to ivar of y
    ivar /= np.power(wavelength, 2)
    print("loaded data")

    #get subset of i100:
    if boss and (mode == '2d' or mode == 'iris'):
        start_idx =  num_columns[j][0]
        num_elements = num_columns[j][1]
        i100_sub = i100[:, start_idx:start_idx+num_elements]        
    else:
        i100_sub = i100

    if bootstrap:
        yxsig_partial = np.zeros(0)
        xxsig_partial = np.zeros(0)
        for p in unique_plates:
            i100_sub_p = i100_sub[plate==p]
            plate_p = plate[plate==p]
            flambda_p = flambda[plate==p]
            ivar_p = ivar[plate==p]
            print("shapes before calc_alphas:")
            print(i100_sub_p.shape)
            print(plate_p.shape)
            print(flambda_p.shape)
            print(ivar_p.shape)
            yxsig_p, xxsig_p = calc_alphas(i100_sub_p, plate_p, flambda_p, ivar_p)
            print(yxsig_p.shape)
            if len(yxsig_partial) == 0:
                yxsig_partial = yxsig_p.reshape(1, yxsig_p.shape[0])
                xxsig_partial = xxsig_p.reshape(1, xxsig_p.shape[0])
            else:
                yxsig_partial = np.append(yxsig_partial, yxsig_p.reshape(1, yxsig_p.shape[0]), axis=0)
                xxsig_partial = np.append(xxsig_partial, xxsig_p.reshape(1, xxsig_p.shape[0]), axis=0)
                print("shape:", yxsig_partial.shape)
        yxsigs_bootstrap = np.append(yxsigs_bootstrap, yxsig_partial, axis=1)
        xxsigs_bootstrap = np.append(xxsigs_bootstrap, xxsig_partial, axis=1)
        print("bigger shape:", yxsigs_bootstrap.shape)

    else:
        #calculate alphas
        alphas_i, alpha_std_i, wavelength_i = calc_alphas(i100_sub, plate, flambda, ivar)
        alphas_10 = np.append(alphas_10, alphas_i)
        alpha_std_10 = np.append(alpha_std_10, alpha_std_i)
        wavelength_10 = np.append(wavelength_10, wavelength_i)

if bootstrap:
    np.save('yx_bootstrap_' + save, yxsigs_bootstrap)
    np.save('xx_bootstrap_' + save, xxsigs_bootstrap)
    
if save != '0' and not bootstrap:
    #np.save('../alphas_and_stds/alphas_test.npy', alphas_10)
    np.save('../alphas_and_stds/alphas_' + save + '.npy', alphas_10)
    np.save('../alphas_and_stds/alpha_stds_' + save + '.npy', alpha_std_10)
    #np.save('../alphas_and_stds/wavelength_boss.npy', wavelength_10)
    print("alphas saved")

'''
#for if I want to bootstrap by fiber
if bootstrap:
    if os.path.isfile(save):
        bootstrap_alphas = np.load(save)  
        print("bootstrap alphas shape:", bootstrap_alphas.shape)
        np.save(save, np.append(bootstrap_alphas, alphas_10.reshape(1, len(alphas_10)), axis=0))
    else:
        np.save(save, alphas_10.reshape(1, len(alphas_10)))
'''



