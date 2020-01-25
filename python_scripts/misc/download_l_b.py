#save SDSS coordinates in npy file (galactic)

import numpy as np

fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
l_b = np.array(fiberinfo[:, 2:4])    # Galactic longitude/latitude, degrees

#np.save("l_b_original.npy", l_b)


'''
#print out lat and long:
with open("infile.txt", 'w') as f:
    for i in range(l_b.shape[0]):
        f.write(str(l_b[i,0]) + " " + str(l_b[i,1]) + "\n")
'''
