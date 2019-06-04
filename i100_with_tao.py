#multiply i100 by correction factor that takes optical depth into account
#i100 was calculated from E(B-V) values
#note: it does not use updated values of i100

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def convert_to_i100():
    #load data
    fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
    i100 = np.array(fiberinfo[:, 4])  # 100 micron intensity (MJy/sr)

    #convert i100 using optical depth:
    taos = np.load("/Users/blakechellew/Documents/DustProject/taos.npy")
    i100_factor = np.divide(1 - np.exp(np.negative(taos)), taos)
    i100 = i100.reshape(len(i100), 1)
    i100 = i100 * i100_factor

    np.save("/Users/blakechellew/Documents/DustProject/i100_tao.npy", i100)

#using small taos (0.01) to see limiting case
def convert_to_i100_test():
    #load data
    fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
    i100 = np.array(fiberinfo[:, 4])  # 100 micron intensity (MJy/sr)

    #convert i100 using optical depth:
    taos = np.load("/Users/blakechellew/Documents/DustProject/taos_small.npy")
    i100_factor = np.divide(1 - np.exp(np.negative(taos)), taos)
    i100 = i100.reshape(len(i100), 1)
    i100 = i100 * i100_factor

    np.save("/Users/blakechellew/Documents/DustProject/i100_tao_small.npy", i100)

#take i100 IRIS values and E(B-V) values, combine.
def convert_to_i100_iris():
    #load data
    i100 = np.load('iris_i100_at_sfd.npy')

    #convert i100 using optical depth:
    taos = np.load("/Users/blakechellew/Documents/DustProject/taos.npy")
    i100_factor = np.divide(1 - np.exp(np.negative(taos)), taos)
    i100 = i100.reshape(len(i100), 1)
    i100 = i100 * i100_factor

    np.save("/Users/blakechellew/Documents/DustProject/i100_iris_tao.npy", i100)

def load_data():
    #compare to old i100
    fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
    i100_old = np.array(fiberinfo[:, 4])  # 100 micron intensity (MJy/sr)
    hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
    wavelength = np.array(hdulist[1].data)  # Angstroms
    
    i100 = np.load("/Users/blakechellew/Documents/DustProject/i100_tao.npy")

    #want how much high and low wavelength values differ on average
    #high:
    high = i100[:,len(wavelength)-1]
    low = i100[:,0]
    high_avg = np.mean(np.subtract(high, i100_old))
    low_avg = np.mean(np.subtract(low, i100_old))
    print("high mean: ", high_avg)
    print("low mean: ", low_avg)

    '''
    #print i100 at each wavelength
    for idx, w in enumerate(wavelength):
        i100_new = i100[:,idx]
        print(i100_new)
    print(i100_old)
    '''
        
    
    
#convert_to_i100()
#load_data()
#convert_to_i100_test()

#convert_to_i100()
convert_to_i100_iris()
