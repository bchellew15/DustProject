from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from time import time
from time import sleep
import sys
from math import floor #for calculating bin ranges
# for masking ecliptic
import astropy.units as u
from astropy.coordinates import SkyCoord

#previously named reproduce_figs_ahab.py
#generate correlation spectra
#this version uses a lot of RAM (~15 GB)
#IMPORTANT: this code does NOT divide by flux conversion factor. For SDSS this is 1.38.
#IMPORTANT: if bootstrap is true, this script does NOT generate bootstrap samples.
#  Instead, it does the first part of the calculation and saves the result as files starting with xx_ and yy_.
#  The script bootstrap_ahab.py needs to be run on those to produce bootstrap samples.
#  The script script_looper.py is set up to do all of this together.
# bootstrap: masked plates still show up in output files as rows of all 0


'''
EXPLAIN PARAMETERS:

mode [1d, 2d, iris_1d, iris_2d]:
  1d refers to linear model where 100 micron intensity is proportional to optical intensity.
  2d refers to updated model that includes optical depth.
  '1d' and '2d' use 100 micron values from SFD, iris_1d and iris_2d use 100 micron values from IRIS
boss: 0 if optical intensity is from SDSS-II, 1 if BOSS
save: 0 to not save anything, or enter savekey to save alphas. if bootstrap = 1, 'save' is ignored.
threshold: helps determine which plates and fibers to mask. plates with avg 100 micron intensity above threshold
  are masked, and any other fibers with intensity above threshold.
location:
  0: uses sky fibers in all regions
  1: uses only fibers in the north
  2: uses fibers in the south 
bootstrap:
  0: generate alphas normally without bootstrapping
  1: save an intermediate step for later bootstrapping [handled by bootstrap_ahab.py]
'''

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
            print("exiting due to excessive RAM usage (" + str(mem_usage) + " GB")
            _thread.interrupt_main()   #threading.currentThread()
        sleep(4)
        
thr1 = threading.Thread(target = check_memory, daemon=True)
thr1.start()

#command line args
#see above for explanation
if len(sys.argv) != 8:
    print("Usage: generate_alphas_ahab.py [mode: 1d, 2d, iris_1d, iris_2d] [boss: 0, 1] [save: 0, savekey] [threshold=10] [location=0, 1, 2] [bootstrap=0] [get_correction=0]")
    exit(0)
mode = sys.argv[1]
boss = int(sys.argv[2])
save = sys.argv[3]
threshold = float(sys.argv[4])
location = int(sys.argv[5])
bootstrap = int(sys.argv[6])
get_correction = int(sys.argv[7])
mask_ecliptic = True
    
#calculate x, y, alpha
def calc_alphas(i100, plate, flambda, ivar, boot=False, i100_old=None):

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
    del x, y
    yxsig = np.multiply(yx, ivar)
    del yx
    xxsig = np.multiply(xx, ivar)
    del xx
    sums1 = np.sum(yxsig, axis=0)
    sums2 = np.sum(xxsig, axis=0)

    avg_correction = None
    avg_i100 = None
    if get_correction:
        correction_factor = i100 / i100_old
        # weighted_avg = (xxsig * correction_factor) / sums2
        # avg_correction = np.sum(weighted_avg, axis=0)
        avg_correction = np.mean(correction_factor, axis=0)
        weighted_i100 = (xxsig * i100_old) / sums2
        avg_i100 = np.sum(weighted_i100, axis=0)

        """
        print("SHAPES")
        print(weights_1col.shape)
        print(sums2_1col.shape)
        print(correction_factor.shape)
        print(weighted_avg.shape)
        print(avg_correction.shape)
        print(weighted_i100.shape)
        print(avg_i100.shape)
        """

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
    return alphas, alpha_std, wavelength, avg_correction, avg_i100


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
correction_10 = np.zeros(0)
i100_10 = np.zeros(0)

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
        if np.mean(i100_old[plate==p]) > threshold:
            ivar[plate==p] = 0
            print("masking whole plate")
            # i100_old[plate==p] = np.nan  # TEST avg i100

    # # TEST: calculate avg i100
    # i100_old[i100_old > threshold] = np.nan
    # # i100_old[ivar < 0] = np.nan
            
    if boss:
        # possible bad CCD columns: 40, 59, 60, 833, 839, and 840.
        ivar[fiber_id == 40] = 0
        ivar[fiber_id == 59] = 0
        ivar[fiber_id == 60] = 0
        ivar[fiber_id == 833] = 0
        ivar[fiber_id == 839] = 0
        ivar[fiber_id == 840] = 0
        # plate 36: the fiber at idx 3178 is bad. it was causing a discrepancy with SFD vs. IRIS results
        ivar[plate == np.unique(plate)[36]] = 0  # doesn't show up in jacknife bc already masked
        ivar[plate == np.unique(plate)[938]] = 0  # was masking fiber 252 which was wrong. now mask whole plate.
        ivar[plate == np.unique(plate)[1509]] = 0
        ivar[plate == np.unique(plate)[1265]] = 0
        ivar[plate == np.unique(plate)[1786]] = 0
        ivar[plate == np.unique(plate)[2141]] = 0
        ivar[plate == np.unique(plate)[2380]] = 0
        ivar[plate == np.unique(plate)[2383]] = 0
        ivar[plate == np.unique(plate)[2388]] = 0

        # # TEST:
        # i100_old[fiber_id == 40] = np.nan
        # i100_old[fiber_id == 59] = np.nan
        # i100_old[fiber_id == 60] = np.nan
        # i100_old[fiber_id == 833] = np.nan
        # i100_old[fiber_id == 839] = np.nan
        # i100_old[fiber_id == 840] = np.nan
        # # plate 36: the fiber at idx 3178 is bad. it was causing a discrepancy with SFD vs. IRIS results
        # i100_old[plate == np.unique(plate)[36]] = np.nan  # doesn't show up in jacknife bc already masked
        # i100_old[plate == np.unique(plate)[938]] = np.nan  # was masking fiber 252 which was wrong. now mask whole plate.
        # i100_old[plate == np.unique(plate)[1509]] = np.nan
        # i100_old[plate == np.unique(plate)[1265]] = np.nan
        # i100_old[plate == np.unique(plate)[1786]] = np.nan
        # i100_old[plate == np.unique(plate)[2141]] = np.nan
        # i100_old[plate == np.unique(plate)[2380]] = np.nan
        # i100_old[plate == np.unique(plate)[2383]] = np.nan
        # i100_old[plate == np.unique(plate)[2388]] = np.nan

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
            l_abs = np.copy(l)
            l_abs[l_abs > 180] -= 360  # make it range -180 to 180, for easier comparison
            l_abs = np.abs(l_abs)  # now only positive
            b_abs = np.abs(b)

            # diagnostics / histograms
            # plt.hist(l_abs, bins=100)
            # plt.title("Galactic Longitudes")
            # plt.show()
            # plt.hist(b, bins=100)
            # plt.title("Galactic Latitudes")
            # plt.show()
            # print("percentiles:")
            # print("longitude:")
            # print(np.percentile(np.abs(l_abs), [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]))
            # print("latitude:")
            # print(np.percentile(np.abs(b), [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]))
            # exit(0)

            if mask_ecliptic:
                c = SkyCoord(l=l * u.degree, b=b * u.degree, frame='galactic')
                ecliptic_coords = c.geocentricmeanecliptic
                ecliptic_longs = ecliptic_coords.lon.deg
                ecliptic_lats = ecliptic_coords.lat.deg

                for p in np.unique(plate):
                    if abs(np.mean(ecliptic_lats[plate == p])) < 10:
                        ivar[plate == p] = 0
                # skipping i100 masking
                # ivar[abs(ecliptic_lats) < 10] = 0  # to mask fibers instead of plates

                # # diagnostics:
                # # plt.hist(ecliptic_lats)
                # # plt.show()
                # num_within_10 = len(ecliptic_lats[abs(ecliptic_lats) < 10])
                # total_num = len(ecliptic_lats)
                # # north:
                # north_num_within_10 = len(ecliptic_lats[(abs(ecliptic_lats) < 10) & (b > 0)])
                # total_north = len(ecliptic_lats[b > 0])
                # south_num_within_10 = len(ecliptic_lats[(abs(ecliptic_lats) < 10) & (b < 0)])
                # total_south = len(ecliptic_lats[b < 0])
                # print("num within 10:", num_within_10)
                # print("fraction within 10:", num_within_10 / total_num)
                # print("north")
                # print("north within 10:", north_num_within_10)
                # print("fraction north within 10:", north_num_within_10 / total_north)
                # print("south within 10:", south_num_within_10)
                # print("fraction south within 10:", south_num_within_10 / total_south)

            if location == 1:
                ivar[b < 0] = 0 #north
                # i100_old[b < 0] = np.nan
            elif location == 2:
                ivar[b > 0] = 0 #south
                # i100_old[b > 0] = np.nan
            elif location == 3:
                # this is towards galactic center (l = 0)
                ivar[(l > 80) & (l < 280)] = 0
            elif location == 4:
                # this is away from galactic center (l = 180)
                ivar[l < 100] = 0
                ivar[l > 260] = 0
            elif location == 5:
                # latitude 30 to 51
                for p in np.unique(plate):
                    if np.mean(b_abs[plate == p]) < 30:
                        ivar[plate == p] = 0
                    if np.mean(b_abs[plate == p]) > 51:
                        ivar[plate == p] = 0
            elif location == 6:
                # latitude 51 to 71
                for p in np.unique(plate):
                    if np.mean(b_abs[plate == p]) < 51:
                        ivar[plate == p] = 0
                    if np.mean(b_abs[plate == p]) > 71:
                        ivar[plate == p] = 0
            elif location == 7:
                # latitude 30 to 41
                for p in np.unique(plate):
                    if np.mean(b_abs[plate == p]) < 30:
                        ivar[plate == p] = 0
                    if np.mean(b_abs[plate == p]) > 41:
                        ivar[plate == p] = 0
            elif location == 8:
                # latitude 41 to 51
                for p in np.unique(plate):
                    if np.mean(b_abs[plate == p]) < 41:
                        ivar[plate == p] = 0
                    if np.mean(b_abs[plate == p]) > 51:
                        ivar[plate == p] = 0
            elif location == 9:
                # latitude 51 to 59
                for p in np.unique(plate):
                    if np.mean(b_abs[plate == p]) < 51:
                        ivar[plate == p] = 0
                    if np.mean(b_abs[plate == p]) > 59:
                        ivar[plate == p] = 0
            elif location == 10:
                # latitude 59 to 71
                for p in np.unique(plate):
                    if np.mean(b_abs[plate == p]) < 59:
                        ivar[plate == p] = 0
                    if np.mean(b_abs[plate == p]) > 71:
                        ivar[plate == p] = 0
            elif location == 11:
                # longitude 36 to 116
                for p in np.unique(plate):
                    if np.mean(l_abs[plate == p]) < 36:
                        ivar[plate == p] = 0
                    if np.mean(l_abs[plate == p]) > 116:
                        ivar[plate == p] = 0
            elif location == 12:
                # longitude 116 to 171
                for p in np.unique(plate):
                    if np.mean(l_abs[plate == p]) < 116:
                        ivar[plate == p] = 0
                    if np.mean(l_abs[plate == p]) > 171:
                        ivar[plate == p] = 0
            elif location == 13:
                # longitude 36 to 84
                for p in np.unique(plate):
                    if np.mean(l_abs[plate == p]) < 36:
                        ivar[plate == p] = 0
                    if np.mean(l_abs[plate == p]) > 84:
                        ivar[plate == p] = 0
            elif location == 14:
                # longitude 84 to 116
                for p in np.unique(plate):
                    if np.mean(l_abs[plate == p]) < 84:
                        ivar[plate == p] = 0
                    if np.mean(l_abs[plate == p]) > 116:
                        ivar[plate == p] = 0
            elif location == 15:
                # longitude 116 to 145
                for p in np.unique(plate):
                    if np.mean(l_abs[plate == p]) < 116:
                        ivar[plate == p] = 0
                    if np.mean(l_abs[plate == p]) > 145:
                        ivar[plate == p] = 0
            elif location == 16:
                # longitude 145 to 171
                for p in np.unique(plate):
                    if np.mean(l_abs[plate == p]) < 145:
                        ivar[plate == p] = 0
                    if np.mean(l_abs[plate == p]) > 171:
                        ivar[plate == p] = 0
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
        i100_old_sub = None
        if get_correction:
            if mode == 'iris_2d':
                i100_old_sub = np.load("../data/iris_i100_at_boss.npy")[:, None]
            else:
                i100_old_sub = i100_old[:, None]
    else:
        i100_sub = i100[:, None]
        i100_old_sub = i100_old[:, None]
        if get_correction and (mode == '2d' or mode == 'iris_2d'):  # must be SDSS if it gets here
            if mode == 'iris_2d':
                i100_old_sub = np.load("../data/iris_i100_at_sfd.npy")[:, None]
            else:
                i100_old_sub = i100_old[:, None]

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
        alphas_i, alpha_std_i, wavelength_i, correction_i, i100_i = calc_alphas(i100_sub, plate, flambda, ivar, i100_old=i100_old_sub)
        alphas_10 = np.append(alphas_10, alphas_i)
        alpha_std_10 = np.append(alpha_std_10, alpha_std_i)
        wavelength_10 = np.append(wavelength_10, wavelength_i)
        correction_10 = np.append(correction_10, correction_i)
        i100_10 = np.append(i100_10, i100_i)
        
print("saving alphas")
    
if bootstrap:
    np.save('../data/yx_bootstrap_' + save, yxsigs_bootstrap)
    np.save('../data/xx_bootstrap_' + save, xxsigs_bootstrap)
elif save != '0':
    np.save('../alphas_and_stds/alphas_' + save + '.npy', alphas_10)
    np.save('../alphas_and_stds/alpha_stds_' + save + '.npy', alpha_std_10)
    #np.save('../alphas_and_stds/wavelength_boss.npy', wavelength_10)
    print("alphas saved")
if save != 0 and get_correction:
    np.save('../alphas_and_stds/correction_factor_' + save + '.npy', correction_10)
    # np.save('../alphas_and_stds/avg_i100_' + save + '.npy', i100_10)



