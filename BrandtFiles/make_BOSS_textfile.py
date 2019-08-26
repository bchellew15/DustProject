#take sky_radec.dat and extract latitude and longitude
#output text file with just those

#1st coordinate column is 0 through 360, so it's longitude, then latitude.
#based on the Readme, this is same order as SFD code.

import numpy as np
from astropy.coordinates import SkyCoord


all_info = np.loadtxt("sky_radec.dat")
all_info = all_info[:,3:]

'''
print(np.min(all_info[:,0]))
print(np.max(all_info[:,0]))

print(np.min(all_info[:,1]))
print(np.max(all_info[:,1]))
'''


#no conversion:
#np.savetxt("BOSS_locations.txt", all_info)

#convert to galactic coordinates
c = SkyCoord(all_info[:,0], all_info[:,1], frame='fk5', unit="deg")

#print(c.fk5.ra[0])
#print(c.fk5.dec[0])
#print(c.galactic.l[0])
#print(c.galactic.l[0])


#make textfile:

longitudes = np.array(c.galactic.l.deg)
latitudes = np.array(c.galactic.b.deg)
longitudes = np.reshape(longitudes, (longitudes.shape[0], 1))
latitudes = np.reshape(latitudes, (latitudes.shape[0], 1))

all_info = np.concatenate((longitudes, latitudes), axis=1)

np.savetxt("BOSS_locations_galactic.txt", all_info)


'''
#convert the first few lines to various coordinate frames, for comparison
for i in range(10):
    #print(c.fk5.ra[i].deg, " ", c.fk5.dec[i].deg)
    print(c.fk4.ra[i].deg/15, " ", c.fk4.dec[i].deg/3600)
    #print(c.galactic.l.deg[i], " ", c.galactic.b.deg[i])
'''
    





