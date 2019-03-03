import numpy as np
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery

#load coordinates of the fibers
fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
l = np.array(fiberinfo[:, 2])     # Galactic longitude, degrees
b = np.array(fiberinfo[:, 3])     # Galactic latitude, degrees

coords = SkyCoord(l, b, unit='deg', frame='galactic')
sfd = SFDQuery()
ebv = sfd(coords)

print("E(B-V): ", ebv)
print("shape: ", ebv.shape)

#np.save('filepath', ebv)
#np.load(filepath)

