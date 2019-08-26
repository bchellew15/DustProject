#read the skyfibers.dat file, try to find out what units the coordinates are

import numpy as np
import matplotlib.pyplot as plt

coords = np.genfromtxt('skyfibers.dat', usecols=(5,6))

print(np.min(coords[:,0]))
print(np.max(coords[:,0]))

print(np.min(coords[:,1]))
print(np.max(coords[:,1]))

plt.hist(coords[:,0])
plt.show()

plt.hist(coords[:,1], bins=100, range=(0, 5))
plt.show()

#same for sky_radec:
coords2 = np.genfromtxt('sky_radec.dat', usecols=(3,4))

plt.hist(coords2[:,0])
plt.show()

plt.hist(coords2[:,1])
plt.show()
