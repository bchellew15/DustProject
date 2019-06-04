#get E(B-V) from SFD dustmaps at locations of the fibers from original paper
#save as .npy file
#this will be converted to tao as function of position and wavelength

import numpy as np
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery

def get_dustmap():
    #load coordinates of the fibers
    fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
    l = np.array(fiberinfo[:, 2])     # Galactic longitude, degrees
    b = np.array(fiberinfo[:, 3])     # Galactic latitude, degrees
    
    coords = SkyCoord(l, b, unit='deg', frame='galactic')
    sfd = SFDQuery()
    ebv = sfd(coords)
    
    print("E(B-V): ", ebv)
    print("shape: ", ebv.shape)
    print("avg: ", np.mean(ebv))

    #WARNING E(B-V) boss might be wrong because of this, should not have _boss
    np.save('/Users/blakechellew/Documents/DustProject/ebv_boss.npy', ebv)



def get_boss_dustmap():
    #now get dustmap from BOSS fibers

    with open("/Users/blakechellew/Documents/DustProject/BrandtFiles/sky_radec.dat") as f:
        lines = f.readlines()
        mjd = [line.split()[0].strip() for line in lines]
        plate = [line.split()[1].strip() for line in lines]
        fiber = [line.split()[2].strip() for line in lines]
        ra = [line.split()[3].strip() for line in lines]
        dec = [line.split()[4].strip() for line in lines]

    coords = SkyCoord(ra, dec, unit='deg', frame='fk5') #J2000 equatorial coordinates
    sfd = SFDQuery()
    ebv = sfd(coords)

    print("E(B-V): ", ebv)
    print("shape: ", ebv.shape)
    print("avg: ", np.mean(ebv))

    np.save('/Users/blakechellew/Documents/DustProject/ebv_boss.npy', ebv)
    #np.load(filepath)

def load_data():
    ebv = np.load('/Users/blakechellew/Documents/DustProject/ebv.npy')
    for e in ebv:
        print(e)
    print("shape: ", ebv.shape)
    print("avg: ", np.mean(ebv))

get_boss_dustmap()
#get_dustmap()
#load_data()




