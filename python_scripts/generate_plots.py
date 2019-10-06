#generate plots of correlation spectra (overall and certain sections) for comparison
#this functionality was previously part of equiv_width_update.py

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys #for command line args
from math import floor #for binning range

#command line options
if len(sys.argv) != 4:
    print("Usage: equiv_width.py [boss: 0, 1] [save: 0, 1] [save_thresh: 0, 1]")
    exit(0)
boss = int(sys.argv[1])
save = int(sys.argv[2])
save_thresh = int(sys.argv[3])

# load wavelengths
if boss:
    wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')
else:
    hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
    wavelength = np.array(hdulist[1].data)

# load in npy files
# original, tao, tao AND iris, iris

#boss alphas:
alphas_boss = [np.load('../alphas_and_stds/alphas_boss_82219.npy'), np.load('../alphas_and_stds/alphas_boss_2d_82219.npy'), \
          np.load('../alphas_and_stds/alphas_boss_iris_82219.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_82219.npy')]
alpha_stds_boss = [np.load('../alphas_and_stds/alphas_boss_stds_82219.npy'), np.load('../alphas_and_stds/alphas_boss_2d_stds_82219.npy'), \
              np.load('../alphas_and_stds/alphas_boss_iris_stds_82219.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_stds_82219.npy')]

#sdss alphas
alphas_sdss = [np.load('../alphas_and_stds/alphas_1d.npy'), np.load('../alphas_and_stds/alphas_2d.npy'), \
          np.load('../alphas_and_stds/alphas_iris.npy'), np.load('../alphas_and_stds/alphas_iris_1d.npy')]
alpha_stds_sdss = [np.load('../alphas_and_stds/alphas_1d_stds.npy'), np.load('../alphas_and_stds/alphas_2d_stds.npy'), \
              np.load('../alphas_and_stds/alphas_iris_stds.npy'), np.load('../alphas_and_stds/alphas_iris_stds_1d.npy')]

#testing:
alphas_boss[0] = np.load('../alphas_and_stds/alphas_boss_91119_10.npy')
alphas_boss[1] = np.load('../alphas_and_stds/alphas_test.npy')

#flux conversion factor:
alphas_sdss = [a/1.38 for a in alphas_sdss]

print("wavelength and alphas:")
print(wavelength)
print(alphas_sdss[0])

if boss:
    alphas = alphas_boss
    alpha_stds = alpha_stds_boss
else:
    alphas = alphas_sdss
    alpha_stds = alpha_stds_sdss

num_arrays = len(alphas)


#plot unbinned spectra (wavelength ranges from paper)
def plot_emissions(alphas1, alphas2, alpha_std1, alpha_std2, label1, label2):
   plt.figure(figsize=(14, 6))

   #plot 4830 - 5040
   plt.subplot(1, 2, 1)
   plt.plot(wavelength, alphas1, c='k', drawstyle='steps', label=label1)
   plt.plot(wavelength, alpha_std1, c='k', drawstyle='steps')
   plt.plot(wavelength, alphas2, c='r', drawstyle='steps', label=label2)
   plt.plot(wavelength, alpha_std2, c='r', drawstyle='steps')

   plt.xlabel(r"Wavelength ($\AA$)")
   plt.ylabel(r"$\alpha_\lambda$")
   plt.legend(loc='upper center', frameon=False)
   plt.xlim(4830, 5040)
   plt.ylim(0, 0.6)
   xcoords = [4863, 4960, 5008]
   for xc in xcoords:
       plt.axvline(x=xc, color='k', linewidth=1, linestyle='--')
       
   #line from 03 continuum::
   plt.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
   #actual continuum for NII:
   plt.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')

   #plot 6530 - 6770 (original vs tao)
   plt.subplot(1, 2, 2)
   plt.plot(wavelength, alphas1, c='k', drawstyle='steps', label=label1)
   plt.plot(wavelength, alpha_std1, c='k', drawstyle='steps')
   plt.plot(wavelength, alphas2, c='r', drawstyle='steps', label=label2)
   plt.plot(wavelength, alpha_std2, c='r', drawstyle='steps')
   
   plt.xlabel(r"Wavelength ($\AA$)")
   plt.ylabel(r"$\alpha_\lambda$")
   plt.legend(loc='upper center', frameon=False)
   plt.xlim(6530, 6770)
   plt.ylim(0, 1.1)
   xcoords = [6550, 6565, 6585, 6718, 6733]
   for xc in xcoords:
       plt.axvline(x=xc, color='k', linewidth=1, linestyle='--')

   #line from 03 continuum::
   plt.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
   #actual continuum for NII:
   plt.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')
       
#plot unbinned spectra:
plot_emissions(alphas[0], alphas[1], alpha_stds[0], alpha_stds[1], "SFD", r"With $\tau$ Correction")
plt.show()
plot_emissions(alphas[0], alphas[3], alpha_stds[0], alpha_stds[3], "SFD", "With IRIS data")
plt.show()
plot_emissions(alphas[0], alphas[2], alpha_stds[0], alpha_stds[2], "SFD", r"With $\tau$ and IRIS" )
plt.show()
       
def generate_binned_alphas(alphas, alpha_stds, wavelength):
    #plot binned alpha vs wavelength (original)
    lambda_range = wavelength[-1] - wavelength[0]
    left_over = lambda_range - 50*floor(lambda_range / 50)  
    binned_lambdas = np.arange(wavelength[0]+left_over/2, wavelength[-2], 50) #[-2] to avoid going over the edge
    binned_alphas = []
    binned_stds = []

    for i in range(len(alphas)):

        binned_alpha_arr = np.zeros(binned_lambdas.shape)
        binned_std_arr = np.zeros(binned_lambdas.shape)
        for j, lmda in enumerate(binned_lambdas):
            indices = np.where((wavelength > lmda) & (wavelength < lmda+50))[0]
            relevant_alphas = alphas[i][indices]
            relevant_stds = alpha_stds[i][indices]
            #weighted average:
            variance = np.power(relevant_stds, 2)
            numerator = np.sum(np.divide(relevant_alphas, variance))
            denominator = np.sum(np.divide(1, variance))
            avg1 = numerator / denominator
            avg2 = 1 / denominator
            binned_alpha_arr[j] = avg1
            binned_std_arr[j] = np.sqrt(avg2)
        binned_alphas.append(binned_alpha_arr)
        binned_stds.append(binned_std_arr)

    return binned_lambdas, binned_alphas, binned_stds

binned_lambdas, binned_alphas, binned_stds = generate_binned_alphas(alphas, alpha_stds, wavelength)

plt.figure(figsize=(6, 14))
if boss:
    y_max = 0.55
    x_min = 3700
    x_max = 10000
else:
    y_max = 0.4
    x_min = 3850
    x_max = 9200

plt.subplot(3, 1, 1)    
#compare original to tao
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[1], c='r', drawstyle='steps', label=r'With $\tau$ Correction')
plt.plot(binned_lambdas, binned_stds[1], c='r', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(3, 1, 2)
#compare original to iris
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[3], c='r', drawstyle='steps', label='With IRIS Data')
plt.plot(binned_lambdas, binned_stds[3], c='r', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(3, 1, 3)
#compare original to tao AND iris
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[2], c='r', drawstyle='steps', label=r'With IRIS and $\tau$')
plt.plot(binned_lambdas, binned_stds[2], c='r', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.tight_layout()
if save:
    plt.savefig("../boss_binned_spectra_82319.png")
plt.show()

plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[2], c='r', drawstyle='steps', label=r'With IRIS and $\tau$')
plt.plot(binned_lambdas, binned_stds[2], c='r', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[1], c='b', drawstyle='steps', label='With tao corr.')
plt.plot(binned_lambdas, binned_stds[1], c='b', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[3], c='g', drawstyle='steps', label=r'with IRIS')
plt.plot(binned_lambdas, binned_stds[3], c='g', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)
plt.show()

'''
#alphas with various thresholds (SDSS, 1d, SFD)
alphas_thresh = [np.load('../alphas_and_stds/alphas_91019_5.npy'), np.load('../alphas_and_stds/alphas_91019_75.npy'), \
               np.load('../alphas_and_stds/alphas_91019_10.npy'), np.load('../alphas_and_stds/alphas_91019_125.npy'), \
                 np.load('../alphas_and_stds/alphas_91019_15.npy'), np.load('../alphas_and_stds/alphas_1d_92719_20.npy'), \
                 np.load('../alphas_and_stds/alphas_1d_92719_25.npy'), np.load('../alphas_and_stds/alphas_1d_92719_30.npy')]
alpha_stds_thresh = [np.load('../alphas_and_stds/alpha_stds_91019_5.npy'), np.load('../alphas_and_stds/alpha_stds_91019_75.npy'), \
               np.load('../alphas_and_stds/alpha_stds_91019_10.npy'), np.load('../alphas_and_stds/alpha_stds_91019_125.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_91019_15.npy'), np.load('../alphas_and_stds/alpha_stds_1d_92719_20.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_1d_92719_25.npy'), np.load('../alphas_and_stds/alpha_stds_1d_92719_30.npy')]
'''

#alphas with various thresholds (SDSS, 2d, SFD)
alphas_thresh = [np.load('../alphas_and_stds/alphas_2d_91119_5.npy'), np.load('../alphas_and_stds/alphas_2d_91119_75.npy'), \
               np.load('../alphas_and_stds/alphas_2d_91119_10.npy'), np.load('../alphas_and_stds/alphas_2d_91119_125.npy'), \
                 np.load('../alphas_and_stds/alphas_2d_91119_15.npy'), np.load('../alphas_and_stds/alphas_2d_92719_20.npy'), \
                 np.load('../alphas_and_stds/alphas_2d_92719_25.npy'), np.load('../alphas_and_stds/alphas_2d_92719_30.npy')]
alpha_stds_thresh = [np.load('../alphas_and_stds/alpha_stds_2d_91119_5.npy'), np.load('../alphas_and_stds/alpha_stds_2d_91119_75.npy'), \
               np.load('../alphas_and_stds/alpha_stds_2d_91119_10.npy'), np.load('../alphas_and_stds/alpha_stds_2d_91119_125.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_2d_91119_15.npy'), np.load('../alphas_and_stds/alpha_stds_2d_92719_20.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_2d_92719_25.npy'), np.load('../alphas_and_stds/alpha_stds_2d_92719_30.npy')]

'''
#alphas with various thresholds (BOSS, 1d, IRIS)
alphas_thresh = [np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_5.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_75.npy'), \
               np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_10.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_125.npy'), \
                 np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_15.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_92719_20.npy'), \
                 np.load('../alphas_and_stds/alphas_boss_iris_1d_92719_25.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_92719_30.npy')]
alpha_stds_thresh = [np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_5.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_75.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_10.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_125.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_15.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_20.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_25.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_30.npy')]

#alphas with various thresholds (BOSS, 2d, IRIS)
alphas_thresh = [np.load('../alphas_and_stds/alphas_boss_iris_91119_5.npy'), np.load('../alphas_and_stds/alphas_boss_iris_91119_75.npy'), \
               np.load('../alphas_and_stds/alphas_boss_iris_91119_10.npy'), np.load('../alphas_and_stds/alphas_boss_iris_91119_125.npy'), \
                 np.load('../alphas_and_stds/alphas_boss_iris_91119_15.npy'), np.load('../alphas_and_stds/alphas_boss_iris_92719_20.npy'), \
                 np.load('../alphas_and_stds/alphas_boss_iris_92719_25.npy'), np.load('../alphas_and_stds/alphas_boss_iris_92719_30.npy')]
alpha_stds_thresh = [np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_5.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_75.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_10.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_125.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_15.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_20.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_25.npy'), \
                     np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_30.npy')]
'''
#plot several thresholds side by side:

binned_lambdas, binned_alphas, binned_stds = generate_binned_alphas(alphas_thresh, alpha_stds_thresh, wavelength)


#5 plots for SFD:
plt.figure(figsize=(14, 3.5))
y_max = 0.45 #0.52 for boss, 0.45 for sdss
x_min = 3700
x_max = 10000

plt.subplot(1, 4, 1)
plt.plot(binned_lambdas, binned_alphas[2], c='k', drawstyle='steps', label='threshold = 10')
plt.plot(binned_lambdas, binned_stds[2], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 4, 2)
plt.plot(binned_lambdas, binned_alphas[4], c='k', drawstyle='steps', label='threshold = 15')
plt.plot(binned_lambdas, binned_stds[4], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 4, 3)
plt.plot(binned_lambdas, binned_alphas[5], c='k', drawstyle='steps', label='threshold = 20')
plt.plot(binned_lambdas, binned_stds[5], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 4, 4)
plt.plot(binned_lambdas, binned_alphas[6], c='k', drawstyle='steps', label='threshold = 25')
plt.plot(binned_lambdas, binned_stds[6], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.tight_layout()
if save_thresh:
    plt.savefig("../sdss_thresholds_2d_92719.png")
plt.show()

