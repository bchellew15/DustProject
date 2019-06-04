#get tao (optical depth) as fn of position and wavelength
#using E(B-V) values calculated from the SFD map
#see get_dustmaps.py

from ebv_to_tao import f99
import numpy as np
from astropy.io import fits

#see top comments
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
    hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
    wavelength = np.array(hdulist[1].data)  # Angstroms

    ebv = np.load("/Users/blakechellew/Documents/DustProject/ebv_boss.npy")
    taos = np.array([f99(wavelength, item) for item in ebv])
    print("taos: ", taos)
    print("shape: ", taos.shape)
    np.save("/Users/blakechellew/Documents/DustProject/taos_boss.npy", taos)

def load_data():
    taos = np.load("/Users/blakechellew/Documents/DustProject/taos.npy")
    for tao in taos:
        print(tao)
    print("taos: ", taos)
    print("Avg: ", np.mean(taos))

#load_data()
get_boss_tao()
#get_tao()




