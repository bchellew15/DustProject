from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from time import time
from time import sleep
import sys
from math import floor #for calculating bin ranges

##########################################
# SEE: generate_alphas_ahab.py
# This is a quick mod to run with different minimum thresholds.
##########################################

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
#see above for explanation
if len(sys.argv) != 8:
    print("Usage: generate_alphas_ahab.py [mode: 1d, 2d, iris_1d, iris_2d] [boss: 0, 1] [save: 0, savekey] [min_thresh=0] [threshold=10] [location=0, 1, 2] [bootstrap=0]")
    exit(0)
mode = sys.argv[1]
boss = int(sys.argv[2])
save = sys.argv[3]
min_thresh = float(sys.argv[4])
threshold = float(sys.argv[5])
location = int(sys.argv[6])
bootstrap = int(sys.argv[7])

#calculate x, y, alpha
def calc_alphas(i100, plate, flambda, ivar, boot=False):

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

    if boot:
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


# load data for BOSS
if boss:

    filenames = ['skyfibers_lam0.fits', 'skyfibers_lam1.fits', 'skyfibers_lam2.fits', 'skyfibers_lam3.fits', 'skyfibers_lam4.fits',
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

    fiber_id = hdulist[7].data

    #get i100 (type: float64)
    i100_old = np.loadtxt("../data/SFD_i100_at_BOSS_locations.txt")[:,2]
    if mode == '2d':
        i100 = np.load("../data/i100_tao_boss.npy")
    elif mode == 'iris_2d':
        i100 = np.load("../data/i100_tao_boss_iris.npy")
    elif mode == 'iris_1d':
        i100 = np.load("../data/iris_i100_at_boss.npy")
    elif mode == '1d':
        i100 = i100_old
    else:
        print("Error: invalid mode")

    plate = 10000*plate + mjd%10000 #plate is now a unique identifier

#load data for SDSS
else:
    num_files = 1

    #load data
    fiberinfo = np.loadtxt('../data/fiberinfo_halpha.dat')

    #get i100 (type: float64)
    i100_old = np.array(fiberinfo[:, 4])# 100 micron intensity (MJy/sr)
    if mode == '2d':
        i100 = np.load("../data/i100_tao.npy") #i100_tao.npy
    elif mode == 'iris_2d':
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

#preprocessing and calculate alphas for each file
for j in range(num_files):
    if boss:
        hdulist = fits.open("../data/" + filenames[j])
        flambda = hdulist[0].data #type: float64
        ivar = hdulist[1].data #type: float64
        wavelength = 10**( min(hdulist[2].data) + np.arange(flambda[0].shape[0])*1e-4 ).astype('float32')

    #process ivar:
    ivar *= (ivar > 0)
    #masking
    ivar[i100_old > threshold] = 0
    #mask plates with large averages
    for p in np.unique(plate):
        mean_i100_p = np.mean(i100_old[plate==p])
        if mean_i100_p > threshold:
            ivar[plate==p] = 0
            print("masking whole plate")
        # lower threshold, to check bias factor
        ivar[(plate == p) & (np.abs(i100_old - mean_i100_p) < min_thresh)] = 0

    if boss:
        ivar[3178] = 0 #data at this location is bad
        ivar[fiber_id == 840] = 0
        ivar[plate == np.unique(plate)[1509]] = 0
        ivar[plate == np.unique(plate)[1265]] = 0
        ivar[plate == np.unique(plate)[1786]] = 0
        ivar[(plate == np.unique(plate)[938]) * (fiber_id == 252)] = 0

    if location != 0:
        if boss:
            coords = np.loadtxt('BOSS_locations_galactic.txt')

            #find locations of problem plates:
            #print("plate 5896")
            #print(np.mean(coords[plate == 5896], axis = 0))
            #print("plate 5896")
            #print(np.mean(coords[plate == 5400], axis=0))
            #print("plate 5896")
            #print(np.mean(coords[plate == 6297], axis=0))
            #exit(0)

            l = coords[:,0]
            b = coords[:,1]
            if location == 1:
                ivar[b < 0] = 0 #north
            elif location == 2:
                ivar[b > 0] = 0 #south
        else:
            coords = np.loadtxt('infile.txt')
            l = coords[:,0]
            b = coords[:,1]
            if location == 1:
                ivar[b < 0] = 0 #north
            elif location == 2:
                ivar[b > 0] = 0 #south

    #convert ivar to ivar of y
    ivar /= np.power(wavelength, 2)
    print("loaded data")

    #get subset of i100:
    if boss and (mode == '2d' or mode == 'iris_2d'):
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
            yxsig_p, xxsig_p = calc_alphas(i100_sub_p, plate_p, flambda_p, ivar_p, boot=True)
            print(yxsig_p.shape)
            if len(yxsig_partial) == 0:
                yxsig_partial = yxsig_p.reshape(1, yxsig_p.shape[0])
                xxsig_partial = xxsig_p.reshape(1, xxsig_p.shape[0])
            else:
                yxsig_partial = np.append(yxsig_partial, yxsig_p.reshape(1, yxsig_p.shape[0]), axis=0)
                xxsig_partial = np.append(xxsig_partial, xxsig_p.reshape(1, xxsig_p.shape[0]), axis=0)
                print("current bootstrap shape:", yxsig_partial.shape)
        yxsigs_bootstrap = np.append(yxsigs_bootstrap, yxsig_partial, axis=1)
        xxsigs_bootstrap = np.append(xxsigs_bootstrap, xxsig_partial, axis=1)
        print("final bootstrap shape:", yxsigs_bootstrap.shape)

    else:
        #calculate alphas
        alphas_i, alpha_std_i, wavelength_i = calc_alphas(i100_sub, plate, flambda, ivar)
        alphas_10 = np.append(alphas_10, alphas_i)
        alpha_std_10 = np.append(alpha_std_10, alpha_std_i)
        wavelength_10 = np.append(wavelength_10, wavelength_i)

# print number of sky fibers
num_unmasked_fibers = np.count_nonzero(ivar[:, 0])
print("unmasked fibers:", num_unmasked_fibers)

print("saving alphas")

if bootstrap:
    np.save('yx_bootstrap_' + save, yxsigs_bootstrap)
    np.save('xx_bootstrap_' + save, xxsigs_bootstrap)
elif save != '0':

    save_path = '../alphas_and_stds/alphas_' + save + '.npy'
    save_num_fibers = '../alphas_and_stds/num_fibers_' + save + '.npy'
    import os.path
    if os.path.isfile(save_path):
        # save alphas
        alphas_from_file = np.load(save_path)
        alphas_10 = np.reshape(alphas_10, (1, len(alphas_10)))
        appended = np.concatenate((alphas_from_file, alphas_10))
        np.save(save_path, appended)

        # save num fibers
        num_fibers_arr = np.load(save_num_fibers)
        num_fibers_arr = np.append(num_fibers_arr, num_unmasked_fibers)
        np.save(save_num_fibers, num_fibers_arr)
    else:
        alphas_10 = np.reshape(alphas_10, (1, len(alphas_10)))
        np.save(save_path, alphas_10)
        np.save(save_num_fibers, [num_unmasked_fibers])

    # np.save('../alphas_and_stds/alpha_stds_' + save + '.npy', alpha_std_10)
    # np.save('../alphas_and_stds/wavelength_boss.npy', wavelength_10)
    print("alphas saved")



