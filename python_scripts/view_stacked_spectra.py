#look at stacked spectra
#(an array of spectra where each row is a different spectrum)

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from generate_plots import generate_binned_alphas

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

"""
# some plots
plt.plot(boss_wavelength, north[1000, :], 'k', label='regular')
plt.plot(boss_wavelength, north[1271, :], 'r', label='bad')
plt.legend()
plt.show()
exit(0)
"""


# view binned spectra:

# 1509 included even though not as significant
idx_array_north = [1509, 378, 1271, 211, 938, 1265, 2141]
idx_array_south = [22, 1081, 1786, 2380, 2381, 2383, 1393, 2011, 2388, 2416]

# for 6 sigma:
idx_array_north = [378, 1271, 2498, 211, 212, 905, 938, 1265, 1509, 2141, 2495, 2500, 2502]
idx_array_south = [22, 23, 389, 419, 440, 536, 667, 1065, 1081, 1106, 1684, 1690, 1739, 1786,
                   2229, 2380, 2381, 2383, 2388, 74, 1078, 1104, 1110, 1393, 1696, 1990, 1997,
                   2009, 2011, 2382, 2416]

# override to check just one
# idx_array_north = [2141]
# idx_array_south = []

for bad_idx in idx_array_north:

    wav_binned, alphas_binned_1, _ = generate_binned_alphas([north[bad_idx]], [np.ones(len(boss_wavelength))], boss_wavelength)
    wav_binned, alphas_binned_2, _ = generate_binned_alphas([north[1000]], [np.ones(len(boss_wavelength))], boss_wavelength)
    alphas_binned_1 = alphas_binned_1[0]
    alphas_binned_2 = alphas_binned_2[0]
    plt.plot(wav_binned, alphas_binned_1, 'k', label='good')
    plt.plot(wav_binned, alphas_binned_2, 'r', label='bad')
    plt.legend()
    plt.ylim(0, 1)
    plt.title("North idx " + str(bad_idx))
    plt.show()

for bad_idx in idx_array_south:
    wav_binned, alphas_binned_3, _ = generate_binned_alphas([south[bad_idx]], [np.ones(len(boss_wavelength))], boss_wavelength)
    wav_binned, alphas_binned_4, _ = generate_binned_alphas([south[1000]], [np.ones(len(boss_wavelength))], boss_wavelength)
    alphas_binned_3 = alphas_binned_3[0]
    alphas_binned_4 = alphas_binned_4[0]
    plt.plot(wav_binned, alphas_binned_3, 'k', label='good')
    plt.plot(wav_binned, alphas_binned_4, 'r', label='bad')
    plt.legend()
    plt.ylim(0, 1)
    plt.title("South idx " + str(bad_idx))
    plt.show()

# and unbinned:
# plt.plot(boss_wavelength, north[1509], 'r', label='bad?')
# plt.plot(boss_wavelength, north[1000], 'k', label='good?')
# plt.legend()
# plt.show()

exit(0)


# get average values in various wavelength ranges
start_wavs = [4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500]
end_wavs = [4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000]
north_avgs = np.zeros((south.shape[0], len(start_wavs)))
south_avgs = np.zeros((south.shape[0], len(start_wavs)))
for i in range(south.shape[0]):
    for j in range(len(start_wavs)):
        indices = np.where((boss_wavelength > start_wavs[j]) & (boss_wavelength < end_wavs[j]))[0]
        north_avgs[i, j] = np.nanmean(north[i][indices])
        south_avgs[i, j] = np.nanmean(south[i][indices])
# check if any are out of the ordinary
avg_north_avgs = np.mean(north_avgs, axis=0)
avg_south_avgs = np.mean(south_avgs, axis=0)
north_avgs_std = np.std(north_avgs, axis=0)
south_avgs_std = np.std(south_avgs, axis=0)

# north:
print("north too big:")
print(np.argwhere(north_avgs > avg_north_avgs + 6*north_avgs_std))
print("north too small:")
print(np.argwhere(north_avgs < avg_north_avgs - 6*north_avgs_std))

# south:
print("south too big:")
print(np.argwhere(south_avgs > avg_south_avgs + 6*south_avgs_std))
print("south too small:")
print(np.argwhere(south_avgs < avg_south_avgs - 6*south_avgs_std))

print(avg_north_avgs)
print(north_avgs_std)

exit(0)

# old stuff
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


