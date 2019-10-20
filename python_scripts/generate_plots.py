#generate plots of correlation spectra (overall and certain sections) for comparison
#this functionality was previously part of equiv_width_update.py

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys #for command line args
from math import floor #for binning range

#command line options
if len(sys.argv) != 5:
    print("Usage: equiv_width.py [boss: 0, 1] [save: 0, 1] [save_thresh: 0, 1] [bootstrap: 0, 1]")
    exit(0)
boss = int(sys.argv[1])
save = int(sys.argv[2])
save_thresh = int(sys.argv[3])
bootstrap = int(sys.argv[4])

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

#sdss alphas (need to update to correct masking)
alphas_sdss = [np.load('../alphas_and_stds/alphas_91019_10.npy'), np.load('../alphas_and_stds/alphas_2d.npy'), \
          np.load('../alphas_and_stds/alphas_iris.npy'), np.load('../alphas_and_stds/alphas_iris_1d.npy')]
alpha_stds_sdss = [np.load('../alphas_and_stds/alphas_1d_stds.npy'), np.load('../alphas_and_stds/alphas_2d_stds.npy'), \
              np.load('../alphas_and_stds/alphas_iris_stds.npy'), np.load('../alphas_and_stds/alphas_iris_stds_1d.npy')]

bootstrap_alphas = np.load('../alphas_and_stds/bootstrap_alphas_sdss_1d_101819.npy')
print("number of bootstrap samples:")
print(bootstrap_alphas.shape)

#testing:
alphas_boss[0] = np.load('../alphas_and_stds/alphas_boss_91119_10.npy')
alphas_boss[1] = np.load('../alphas_and_stds/alphas_test.npy')
alphas_sdss[1] = bootstrap_alphas[-1]

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


#calculate bootstrap upper / lower bounds
#sdss 1d: alphas_boss_91119_10
#bootstrap samples: bootstrap_alphas_sdss_1d_101819.npy
bootstrap_lower_bound = np.percentile(bootstrap_alphas, 5, axis=0) / 1.38
bootstrap_upper_bound = np.percentile(bootstrap_alphas, 95, axis=0) / 1.38
#bootstrap estimate of std dev, for comparison:
bootstrap_std = (np.percentile(bootstrap_alphas, 84, axis=0) - np.percentile(bootstrap_alphas, 16, axis=0)) / 1.38




#plot unbinned spectra (wavelength ranges from paper)
def plot_emissions(alphas1, alphas2, alpha_std1, alpha_std2, label1, label2):
   plt.figure(figsize=(14, 6))

   #plot 4830 - 5040
   plt.subplot(1, 2, 1)
   plt.plot(wavelength, alphas1, c='k', drawstyle='steps', label=label1)
   plt.plot(wavelength, alpha_std1, c='k', drawstyle='steps')

   if not bootstrap:
       plt.plot(wavelength, alphas2, c='r', drawstyle='steps', label=label2)
       plt.plot(wavelength, alpha_std2, c='r', drawstyle='steps')

   if bootstrap:
       plt.fill_between(wavelength, bootstrap_lower_bound, bootstrap_upper_bound, alpha=0.5, step='pre')
       #plt.plot(wavelength, bootstrap_lower_bound, c='m', drawstyle='steps')
       #plt.plot(wavelength, bootstrap_upper_bound, c='m', drawstyle='steps')
       plt.plot(wavelength, bootstrap_std, c='m', drawstyle='steps')

   plt.xlabel(r"Wavelength ($\AA$)")
   plt.ylabel(r"$\alpha_\lambda$")
   plt.legend(loc='upper center', frameon=False)
   plt.xlim(4830, 5040)
   plt.ylim(0, 0.6)
   xcoords = [4863, 4960, 5008]
   for xc in xcoords:
       plt.axvline(x=xc, color='k', linewidth=1, linestyle='--')
       
   #line from 03 continuum::
   #plt.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
   #actual continuum for NII:
   #plt.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')

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
   #plt.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
   #actual continuum for NII:
   #plt.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')
       
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
if bootstrap:
    _, bootstrap_binned_lower, bootstrap_binned_stds = generate_binned_alphas([bootstrap_lower_bound], [bootstrap_std], wavelength)
    _, bootstrap_binned_upper, _ = generate_binned_alphas([bootstrap_upper_bound], [bootstrap_std], wavelength)
    bootstrap_binned_lower = bootstrap_binned_lower[0]
    bootstrap_binned_upper = bootstrap_binned_upper[0]
    bootstrap_binned_stds = bootstrap_binned_stds[0]
    

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
plt.fill_between(binned_lambdas, bootstrap_binned_lower, bootstrap_binned_upper, alpha=0.5, step='pre')
plt.plot(binned_lambdas, bootstrap_binned_stds, c='m', drawstyle='steps')
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



if boss:
    #alphas with various thresholds (BOSS, 1d, IRIS)
    alphas_thresh_1d = [np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_5.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_75.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_10.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_125.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_15.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_92719_20.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_92719_25.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_92719_30.npy')]
    alpha_stds_thresh_1d = [np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_5.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_75.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_10.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_125.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_15.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_20.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_25.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_30.npy')]
    #alphas with various thresholds (BOSS, 2d, IRIS)
    alphas_thresh_2d = [np.load('../alphas_and_stds/alphas_boss_iris_91119_5.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_91119_75.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_91119_10.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_91119_125.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_91119_15.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_92719_20.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_92719_25.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_92719_30.npy')]
    alpha_stds_thresh_2d = [np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_5.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_75.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_10.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_125.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_15.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_20.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_25.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_30.npy')]
else:
    #alphas with various thresholds (SDSS, 1d, SFD)
    alphas_thresh_1d = [np.load('../alphas_and_stds/alphas_91019_5.npy'), \
                        np.load('../alphas_and_stds/alphas_91019_75.npy'), \
                        np.load('../alphas_and_stds/alphas_91019_10.npy'), \
                        np.load('../alphas_and_stds/alphas_91019_125.npy'), \
                        np.load('../alphas_and_stds/alphas_91019_15.npy'), \
                        np.load('../alphas_and_stds/alphas_1d_92719_20.npy'), \
                        np.load('../alphas_and_stds/alphas_1d_92719_25.npy'), \
                        np.load('../alphas_and_stds/alphas_1d_92719_30.npy')]
    alpha_stds_thresh_1d = [np.load('../alphas_and_stds/alpha_stds_91019_5.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_91019_75.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_91019_10.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_91019_125.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_91019_15.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_1d_92719_20.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_1d_92719_25.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_1d_92719_30.npy')]
    #alphas with various thresholds (SDSS, 2d, SFD)
    alphas_thresh_2d = [np.load('../alphas_and_stds/alphas_2d_91119_5.npy'), \
                        np.load('../alphas_and_stds/alphas_2d_91119_75.npy'), \
                        np.load('../alphas_and_stds/alphas_2d_91119_10.npy'), \
                        np.load('../alphas_and_stds/alphas_2d_91119_125.npy'), \
                        np.load('../alphas_and_stds/alphas_2d_91119_15.npy'), \
                        np.load('../alphas_and_stds/alphas_2d_92719_20.npy'), \
                        np.load('../alphas_and_stds/alphas_2d_92719_25.npy'), \
                        np.load('../alphas_and_stds/alphas_2d_92719_30.npy')]
    alpha_stds_thresh_2d = [np.load('../alphas_and_stds/alpha_stds_2d_91119_5.npy'),
                            np.load('../alphas_and_stds/alpha_stds_2d_91119_75.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_2d_91119_10.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_2d_91119_125.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_2d_91119_15.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_2d_92719_20.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_2d_92719_25.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_2d_92719_30.npy')]


#plot several thresholds side by side:

binned_lambdas, binned_alphas_1d, binned_stds_1d = generate_binned_alphas(alphas_thresh_1d, alpha_stds_thresh_1d, wavelength)
binned_lambdas, binned_alphas_2d, binned_stds_2d = generate_binned_alphas(alphas_thresh_2d, alpha_stds_thresh_2d, wavelength)

fig = plt.figure(figsize=(8, 4), dpi=200)

print(binned_alphas_1d[2])

if boss:
    x_min = 3700
    x_max = 10100
    y_max = .52
else:
    x_min = 3800
    x_max = 9200
    y_max = .41

ax = fig.add_subplot(121)
plt.text(0.02, 0.98, 'Original\nModel', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10, fontweight='bold')

plt.plot(binned_lambdas, binned_alphas_1d[2], c='k', drawstyle='steps', label='threshold 10')
plt.plot(binned_lambdas, binned_stds_1d[2], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)

plt.plot(binned_lambdas, binned_alphas_1d[4], c='b', drawstyle='steps', linestyle='--', label='threshold 15')
plt.plot(binned_lambdas, binned_stds_1d[4], c='b', drawstyle='steps', linestyle='--')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)

plt.plot(binned_lambdas, binned_alphas_1d[5], c='m', drawstyle='steps', linestyle='dashdot', label='threshold 20')
plt.plot(binned_lambdas, binned_stds_1d[5], c='m', drawstyle='steps', linestyle='dashdot')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)

plt.plot(binned_lambdas, binned_alphas_1d[6], c='g', drawstyle='steps', linestyle=':', label='threshold 25')
plt.plot(binned_lambdas, binned_stds_1d[6], c='g', drawstyle='steps', linestyle=':')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)

if boss:
    plt.legend(frameon=False, loc='upper right')
else:
    plt.legend(frameon=False, loc='lower center')

plt.text(0, 100, 'Original\nModel')

ax = fig.add_subplot(122)
plt.text(0.02, 0.98, 'Tao\nModel', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10, fontweight='bold')

plt.plot(binned_lambdas, binned_alphas_2d[2], c='k', drawstyle='steps', label='threshold 10')
plt.plot(binned_lambdas, binned_stds_2d[2], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)

plt.plot(binned_lambdas, binned_alphas_2d[4], c='b', drawstyle='steps', linestyle='--', label='threshold 15')
plt.plot(binned_lambdas, binned_stds_2d[4], c='b', drawstyle='steps', linestyle='--')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)

plt.plot(binned_lambdas, binned_alphas_2d[5], c='m', drawstyle='steps', linestyle='dashdot', label='threshold 20')
plt.plot(binned_lambdas, binned_stds_2d[5], c='m', drawstyle='steps', linestyle='dashdot')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)

plt.plot(binned_lambdas, binned_alphas_2d[6], c='g', drawstyle='steps', linestyle=':', label='threshold 25')
plt.plot(binned_lambdas, binned_stds_2d[6], c='g', drawstyle='steps', linestyle=':')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)

if boss:
    plt.legend(frameon=False, loc='upper right')
else:
    plt.legend(frameon=False, loc='lower center')

plt.text(0, 100, 'Tao Model')

plt.tight_layout()

if save_thresh:
    if boss:
        plt.savefig("../boss_thresholds_2panel_10719.png")
    else:
        plt.savefig("../sdss_thresholds_2panel_10719.png")
plt.show()

