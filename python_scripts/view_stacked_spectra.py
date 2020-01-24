#look at stacked spectra
#(an array of spectra where each row is a different spectrum)

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

np.set_printoptions(threshold=sys.maxsize)

#for the plot
xmin = 7800 #8200 #7900 #5888 #5750
xmax = 8200 #8500 #8100 #5893 #5900
ymin = -1 #0
ymax = 1 #0.6
#for the calculations
xmin_calc = 9369 #8298 #9302 #5888 #5840
xmax_calc = 9376 #8301 #9313 #5893 #5855



#jacknife: plot the first few
south = np.load("../alphas_and_stds/jacknife_alphas_south_111719.npy")
north = np.load("../alphas_and_stds/jacknife_alphas_north_111719.npy")
boss_wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')

'''
south_means = np.nanmean(south, axis=1)
north_means = np.nanmean(north, axis=1)
plt.hist(south_means, 50)
plt.title("south")
plt.show()
plt.hist(north_means, 50)
plt.title("north")
plt.show()

#south: .63 to .7 is normal range
print(np.where(south_means < .63)[0])
print(np.where(south_means > .7)[0])

print(south_means[1265])

#north: .34 to .38 is normal range
print(np.where(north_means < .34)[0])
print(np.where(north_means > .38)[0])

print(north_means[1265])

exit(0)
'''


indices = np.where((boss_wavelength>xmin_calc) * (boss_wavelength<xmax_calc))[0]
#for the dip around 8000:
indices = np.where((boss_wavelength>7921.5) * (boss_wavelength<7922.5) + (boss_wavelength>7978) * (boss_wavelength<7982) + (boss_wavelength>8208) * (boss_wavelength<8230))[0]

south_sum = np.sum(south[:, indices], axis=1)
north_sum = np.sum(north[:, indices], axis=1)

#print("sums")
#print(south_sum[:10])
#print(north_sum[:10])

print("diffs")
diffs = north_sum - south_sum
print(diffs[:10])

print(diffs[diffs<0.5])

#only one is <1, and the value is -0.3
bad_idx = np.where(diffs<0.5)[0]
print(bad_idx)

#plt.hist(diffs, bins=100)
#plt.show()



idx_removed = 1786
plt.plot(boss_wavelength, north[idx_removed], drawstyle='steps', label='north')
plt.plot(boss_wavelength, south[idx_removed], drawstyle='steps', label='south')
plt.legend(loc='upper center', frameon=False)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.show()


idx_removed = 1787
plt.plot(boss_wavelength, north[idx_removed], drawstyle='steps', label='north')
plt.plot(boss_wavelength, south[idx_removed], drawstyle='steps', label='south')
plt.legend(loc='upper center', frameon=False)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.show()


