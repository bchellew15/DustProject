#plot sky locations
#color coded according to density and weighted density (see Brandt paper)

import numpy as np
from matplotlib import pyplot as plt


#verify it works using SFD locations:
coords = np.loadtxt('/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/infile.txt')
longs = coords[:,0]
lats = coords[:,1]

ax = plt.axes(projection='mollweide')
ax.plot(longs, lats, 'k.', markersize=1)
ax.set_title('SFD locations (mollweide)')
ax.grid(True)
plt.show()

for i in range(len(longs)):
    if longs[i] > 180:
        longs[i] -= 360
ax = plt.axes(projection='mollweide')
ax.plot(longs, lats, 'k.', markersize=1)
ax.set_title('SFD locations (mollweide)')
ax.grid(True)
plt.show()
        
exit(0)

coords = np.loadtxt("/Users/blakechellew/Documents/DustProject/BrandtFiles/BOSS_locations_galactic.txt")
longs = coords[:,0]
lats = coords[:,1]

#plt.figure(figsize=(5, 4))

ax = plt.axes(projection='mollweide')
ax.plot(longs, lats, 'k.', markersize=1)
#ax.grid(True)
ax.set_title('Fiber Locations')
plt.show()

ax = plt.axes(projection='hammer')
ax.plot(longs, lats, 'k.', markersize=1)
#ax.grid(True)
ax.set_title('Fiber Locations')
plt.show()

ax = plt.axes(projection='lambert')
ax.plot(longs, lats, 'k.', markersize=1)
#ax.grid(True)
ax.set_title('Fiber Locations')
plt.show()
