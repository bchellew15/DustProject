#get tao (optical depth) as fn of position and wavelength
#using E(B-V) values calculated from the SFD map
#(see get_dustmaps.py)
#next step: combine with i100 values to get 2D i100 array

from ebv_to_tao import f99
import numpy as np
from astropy.io import fits

def get_tao():
    hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
    wavelength = np.array(hdulist[1].data)  # Angstroms

    ebv = np.load("/Users/blakechellew/Documents/DustProject/ebv.npy")

    taos = np.array([f99(wavelength, item) for item in ebv])
    print("taos: ", taos)
    print("shape: ", taos.shape)

    np.save("/Users/blakechellew/Documents/DustProject/taos.npy", taos)

#same but at BOSS locations
def get_boss_tao():

    #need to get i100 at all the wavelengths covered by BOSS: first to last + however many

    hdulist = fits.open('../BrandtFiles/skyfibers_nativelam.fits')
    flam = hdulist[0].data[0]

    #number of elements: (number of 10^-4 intervals beetween min and max starting val) + flam.shape[0]
    num_intervals = int(round((max(hdulist[2].data) - min(hdulist[2].data)) / 0.0001))
    lams = 10**( min(hdulist[2].data) + np.arange(num_intervals + flam.shape[0])*1e-4 )
    #lams has 4968 elements

    # save the full BOSS wavelength array
    # (so far I've been using a partial array, but need this for avg i100)
    np.save("/Users/blakechellew/Documents/DustProject/wavelength_boss_full.npy", lams)
    exit(0)
    
    ebv = np.load("/Users/blakechellew/Documents/DustProject/ebv_boss.npy")
    taos = np.array([f99(lams, item) for item in ebv])
    print("taos: ", taos)
    print("shape: ", taos.shape)
    np.save("/Users/blakechellew/Documents/DustProject/taos_boss.npy", taos)

def load_data():
    taos = np.load("/Users/blakechellew/Documents/DustProject/taos_boss.npy")
    print("taos: ", taos)
    print("Avg: ", np.mean(taos))

# load_data()
get_boss_tao()
#get_tao()




