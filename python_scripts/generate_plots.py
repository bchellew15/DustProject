#generate plots of correlation spectra (overall and certain sections) for comparison
#this functionality was previously part of equiv_width_update.py

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys #for command line args
from math import floor #for binning range

#command line options
if len(sys.argv) != 3:
    print("Usage: equiv_width.py [boss: 0, 1] [save: 0, 1]")
    exit(0)
boss = int(sys.argv[1])
save = int(sys.argv[2])

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

alphas_sdss[0] = [0.18972603, 0.18287671, 0.1630137, 0.14657534, 0.09863014, 0.1609589, 0.18150685, 0.1890411, 0.16986301, 0.17328767, \
                  0.20547945, 0.20205479, 0.20205479, 0.23835616, 0.26438356, 0.24383562, 0.20136986, 0.16986301, 0.1390411, 0.13630137, 0.17876712, \
                  0.1739726, 0.28150685, 0.53835616, 0.61438356, 0.49452055, 0.31780822, 0.17465753, 0.19041096, 0.16643836, 0.16849315, \
                  0.16164384, 0.19383562, 0.19520548, 0.1760274, 0.18356164, 0.24109589, 0.33630137, 0.44315068, 0.4, 0.25, 0.20616438, \
                  0.1760274, 0.16780822, 0.17123288, 0.16712329, 0.17671233, 0.17260274, 0.17328767, 0.20205479, 0.17671233, 0.17191781, \
                  0.17808219, 0.20479452, 0.19383562, 0.1739726, 0.19383562, 0.19863014, 0.18630137, 0.19109589, 0.19246575, 0.18013699, \
                  0.17945205, 0.17671233, 0.16917808, 0.17260274, 0.16369863, 0.17671233, 0.18013699, 0.19726027, 0.16643836, 0.18561644, \
                  0.1739726, 0.17123288, 0.16712329, 0.18219178, 0.18630137, 0.17945205, 0.18287671, 0.17808219, 0.1760274, 0.19383562, \
                  0.16506849, 0.14931507, 0.16506849, 0.16575342, 0.17671233, 0.18630137, 0.18356164, 0.17534247, 0.18561644, 0.1869863, \
                  0.1609589, 0.17808219, 0.19109589, 0.18561644, 0.13424658, 0.17191781, 0.17945205, 0.18561644, 0.17534247, 0.1869863, \
                  0.17945205, 0.17739726, 0.17808219, 0.18630137, 0.18356164, 0.19726027, 0.20342466, 0.17739726, 0.15890411, 0.17260274, \
                  0.17945205, 0.17534247, 0.18082192, 0.16849315, 0.18219178, 0.20205479, 0.19794521, 0.16712329, 0.16849315, 0.18493151, \
                  0.17808219, 0.19520548, 0.17876712, 0.18835616, 0.19246575, 0.26712329, 0.37671233, 0.40410959, 0.3239726, 0.21849315, \
                  0.17465753, 0.19041096, 0.19383562, 0.18287671, 0.20273973, 0.26643836, 0.33972603, 0.34383562, 0.26712329, 0.19863014, \
                  0.17945205, 0.1869863, 0.19657534, 0.17260274, 0.18493151, 0.16369863, 0.17876712, 0.20273973, 0.20890411, 0.1869863, \
                  0.18013699, 0.16849315, 0.16643836, 0.18493151, 0.19041096, 0.17260274, 0.17739726, 0.17328767, 0.15958904, 0.18287671, \
                  0.19657534, 0.18561644, 0.19178082] #from Brandt paper
alphas_sdss[0] = np.pad(alphas_sdss[0], pad_width=((2449, 1386)), mode='constant')
alphas_sdss[3] = np.load('../alphas_and_stds/alphas_sdss_83019.npy') /1.38

print("lengths:")
print(len(alphas_sdss[0]))
print(len(alphas_sdss[3]))
#index 2487 is the NII peak
#it is index 37 in the new array
#pad with 2450 zeros at beginning
#lengths: 164 vs 4000 -> 3836 difference
# -> pad with 1386 zeros at the end


#other alphas:
alphas_misc = [np.load('../alphas_and_stds/alphas_boss_82319_5.npy'), np.load('../alphas_and_stds/alphas_boss_82319_75.npy'), \
               np.load('../alphas_and_stds/alphas_boss_82219.npy'), np.load('../alphas_and_stds/alphas_boss_82319_125.npy'), \
               np.load('../alphas_and_stds/alphas_boss_82319_15.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_82519_5.npy'), \
               np.load('../alphas_and_stds/alphas_boss_iris_1d_82519_75.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_82219.npy'), \
               np.load('../alphas_and_stds/alphas_boss_iris_1d_82519_125.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_82519_15.npy')]
alpha_stds_misc = [np.load('../alphas_and_stds/alphas_boss_stds_82319_5.npy'), np.load('../alphas_and_stds/alphas_boss_stds_82319_75.npy'), \
                   np.load('../alphas_and_stds/alphas_boss_stds_82219.npy'), np.load('../alphas_and_stds/alphas_boss_stds_82319_125.npy'), \
                   np.load('../alphas_and_stds/alphas_boss_stds_82319_15.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_stds_82519_5.npy'), \
                   np.load('../alphas_and_stds/alphas_boss_iris_1d_stds_82519_75.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_stds_82219.npy'), \
                   np.load('../alphas_and_stds/alphas_boss_iris_1d_stds_82519_125.npy'), np.load('../alphas_and_stds/alphas_boss_iris_1d_stds_82519_15.npy')]

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
#plot several thresholds side by side:

binned_lambdas, binned_alphas, binned_stds = generate_binned_alphas(alphas_misc, alpha_stds_misc, wavelength)


#5 plots for SFD:
plt.figure(figsize=(14, 6))
y_max = 0.55
x_min = 3700
x_max = 10000

plt.subplot(1, 5, 1)
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='threshold = 5')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 5, 2)
plt.plot(binned_lambdas, binned_alphas[1], c='k', drawstyle='steps', label='threshold = 7.5')
plt.plot(binned_lambdas, binned_stds[1], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 5, 3)
plt.plot(binned_lambdas, binned_alphas[2], c='k', drawstyle='steps', label='threshold = 10')
plt.plot(binned_lambdas, binned_stds[2], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 5, 4)
plt.plot(binned_lambdas, binned_alphas[3], c='k', drawstyle='steps', label='threshold = 12.5')
plt.plot(binned_lambdas, binned_stds[3], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 5, 5)
plt.plot(binned_lambdas, binned_alphas[4], c='k', drawstyle='steps', label='threshold = 15')
plt.plot(binned_lambdas, binned_stds[4], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.tight_layout()
#plt.savefig("../boss_spectra_overall.png")
plt.show()


#5 plots for IRIS:
plt.figure(figsize=(14, 6))
y_max = 0.55
x_min = 3700
x_max = 10000

plt.subplot(1, 5, 1)
plt.plot(binned_lambdas, binned_alphas[5], c='k', drawstyle='steps', label='threshold = 5')
plt.plot(binned_lambdas, binned_stds[5], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 5, 2)
plt.plot(binned_lambdas, binned_alphas[6], c='k', drawstyle='steps', label='threshold = 10')
plt.plot(binned_lambdas, binned_stds[6], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 5, 3)
plt.plot(binned_lambdas, binned_alphas[7], c='k', drawstyle='steps', label='threshold = 10')
plt.plot(binned_lambdas, binned_stds[7], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 5, 4)
plt.plot(binned_lambdas, binned_alphas[8], c='k', drawstyle='steps', label='threshold = 10')
plt.plot(binned_lambdas, binned_stds[8], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(1, 5, 5)
plt.plot(binned_lambdas, binned_alphas[9], c='k', drawstyle='steps', label='threshold = 15')
plt.plot(binned_lambdas, binned_stds[9], c='k', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.tight_layout()
#plt.savefig("../boss_spectra_overall.png")
plt.show()
'''

