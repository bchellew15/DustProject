'''
check SFD i100 values I got from C code against values I got from Brandt
I ran:
./dust_getval infile=infile.txt map=I100 interp=y verbose=y outfile=outfile.txt
where infile is list of SSDS sky survey coordinates.
'''

import numpy as np
import matplotlib.pyplot as plt

#from C code:
c_output = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/outfile.txt")[:,2]

#from Brandt:
fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
i100 = np.array(fiberinfo[:, 4])  # 100 micron intensity (MJy/sr)

plt.plot(i100, c_output, 'k.')
plt.plot(np.arange(0,200), np.arange(0,200), 'r')
plt.xlabel("SFD i100 (Dataset Values)")
plt.ylabel("SFD i100 (C code)")
plt.show()

#plt.hist(c_output, bins=10)


#check SFD against IRIS for BOSS:
c_output_boss = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/SFD_i100_at_BOSS_locations.txt")[:,2]
iris_boss = np.load("/Users/blakechellew/Documents/DustProject/IRIS/iris_i100_at_boss.npy")

#print(c_output_boss.shape)
#print(iris_boss.shape)


plt.plot(c_output_boss, iris_boss, 'k.')
plt.xlabel("SFD i100 at BOSS")
plt.ylabel("IRIS i100 at BOSS")
plt.plot([x for x in range(1000)], [x for x in range(1000)], 'b')
plt.show()


plt.hist(c_output_boss, bins=10, range=(0, 10))
plt.show()
plt.hist(iris_boss, bins=10, range=(0, 10))
plt.show()
