#generate plots of correlation spectra (overall and certain sections) for comparison
#this functionality was previously part of equiv_width_update.py

#flux conversion factor is handled here, not in reproduce_figs

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys #for command line args
from math import floor #for binning range
from scipy.optimize import curve_fit #for checking if bootstrap histogram is gaussian

#command line options
if len(sys.argv) != 5:
    print("Usage: equiv_width.py [boss: 0, 1] [save: 0, 1] [save_thresh: 0, 1] [bootstrap: 0, 1]")
    exit(0)
boss = int(sys.argv[1])
save = int(sys.argv[2])
save_thresh = int(sys.argv[3])
bootstrap = int(sys.argv[4])

# load wavelengths
wavelength_boss = np.load('../alphas_and_stds/wavelength_boss.npy')
hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
wavelength_sdss = np.array(hdulist[1].data)
if boss:
    wavelength = wavelength_boss
else:
    wavelength = wavelength_sdss

sdss_fluxfactor = 1.38
boss_fluxfactor = 1.38

# load in npy files
# original, tao, tao AND iris, iris

#boss alphas:
alphas_boss = [np.load('../alphas_and_stds/alphas_boss_102019.npy'), \
               np.load('../alphas_and_stds/alphas_boss_2d_102019.npy'), \
               np.load('../alphas_and_stds/alphas_boss_iris_91119_10.npy'), \
               np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_10.npy')]
alpha_stds_boss = [np.load('../alphas_and_stds/alpha_stds_boss_102019.npy'), \
                   np.load('../alphas_and_stds/alpha_stds_boss_2d_102019.npy'), \
                   np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_10.npy'), \
                   np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_10.npy')]
if bootstrap:
    bootstrap_alphas_boss = [np.load('../alphas_and_stds/bootstrap_alphas_boss_1d_102019.npy'), \
                             np.load('../alphas_and_stds/bootstrap_alphas_boss_2d_102019.npy'), \
                             np.load('../alphas_and_stds/bootstrap_alphas_boss_iris_2d_102019.npy'), \
                             np.load('../alphas_and_stds/bootstrap_alphas_boss_iris_1d_102019.npy')]
    bootstrap_lower_boss = [np.percentile(b, 16, axis=0) / boss_fluxfactor for b in bootstrap_alphas_boss] #90 percent confidence interval
    bootstrap_upper_boss = [np.percentile(b, 84, axis=0) / boss_fluxfactor for b in bootstrap_alphas_boss]
    bootstrap_stds_boss= [(np.percentile(b, 84, axis=0) - np.percentile(b, 16, axis=0)) / 2 for b in bootstrap_alphas_boss]

#sdss alphas (need to update to correct masking)
alphas_sdss = [np.load('../alphas_and_stds/alphas_91019_10.npy'), \
               np.load('../alphas_and_stds/alphas_2d_91119_10.npy'), \
               np.load('../alphas_and_stds/alphas_sdss_iris_2d_102019.npy'), \
               np.load('../alphas_and_stds/alphas_sdss_iris_1d_102019.npy'), \
               np.load('../alphas_and_stds/alphas_boss_102019.npy')]
alpha_stds_sdss = [np.load('../alphas_and_stds/alpha_stds_91019_10.npy'), \
                   np.load('../alphas_and_stds/alpha_stds_2d_91119_10.npy'), \
                   np.load('../alphas_and_stds/alpha_stds_sdss_iris_2d_102019.npy'), \
                   np.load('../alphas_and_stds/alpha_stds_sdss_iris_1d_102019.npy'), \
                   np.load('../alphas_and_stds/alpha_stds_boss_102019.npy')]
if bootstrap:
    bootstrap_alphas_sdss = [np.load('../alphas_and_stds/bootstrap_alphas_sdss_1d_101819.npy'), \
                             np.load('../alphas_and_stds/bootstrap_alphas_sdss_2d_102019.npy'), \
                             np.load('../alphas_and_stds/bootstrap_alphas_sdss_iris_2d_102019.npy'), \
                             np.load('../alphas_and_stds/bootstrap_alphas_sdss_iris_1d_102019.npy'), \
                             np.load('../alphas_and_stds/bootstrap_alphas_boss_1d_102019.npy')]

    '''
    #histogram (verify that this is the MLE gaussian)
    data_hist = bootstrap_alphas_sdss[0][2473]     #[int(bootstrap_alphas_sdss[0].shape[1]/2)]
    data_hist = data_hist[(data_hist>0)*(data_hist<1)]

    avg_hist = np.mean(data_hist)
    var_hist = np.var(data_hist)
    pdf_x = np.linspace(0, 0.4, 100)
    pdf_y = 1.0/np.sqrt(2*np.pi*var_hist)*np.exp(-0.5*(pdf_x-avg_hist)**2/var_hist)

    h = plt.hist(data_hist, bins=50, range=(0, 0.4)) #,normed=True)
    plt.plot(pdf_x, np.max(h[0])/np.max(pdf_y)*pdf_y, 'k--')
    plt.show()

    #plt.hist(bootstrap_alphas_sdss[0][:int(bootstrap_alphas_sdss[0].shape[1]/2)], bins=50)

    plt.show()
    exit(0)
    '''
    

    bootstrap_lower_sdss = [np.percentile(b, 16, axis=0) / sdss_fluxfactor for b in bootstrap_alphas_sdss] #90 percent confidence interval
    bootstrap_upper_sdss = [np.percentile(b, 84, axis=0) / sdss_fluxfactor for b in bootstrap_alphas_sdss]
    bootstrap_stds_sdss = [(np.percentile(b, 84, axis=0) - np.percentile(b, 16, axis=0)) / 2 for b in bootstrap_alphas_sdss]
    
#flux conversion factor:
alphas_sdss = [a/sdss_fluxfactor for a in alphas_sdss]
alphas_boss = [a/sdss_fluxfactor for a in alphas_boss] #test

if boss:
    alphas = alphas_boss
    alpha_stds = alpha_stds_boss
    if bootstrap:
        bootstrap_lower = bootstrap_lower_boss
        bootstrap_upper = bootstrap_upper_boss
        bootstrap_stds = bootstrap_stds_boss
else:
    alphas = alphas_sdss
    alpha_stds = alpha_stds_sdss
    if bootstrap:
        bootstrap_lower = bootstrap_lower_sdss
        bootstrap_upper = bootstrap_upper_sdss
        bootstrap_stds = bootstrap_stds_sdss

num_arrays = len(alphas)


#plot unbinned spectra (wavelength ranges from paper)
def plot_emissions(idx1, idx2, label1, label2):
   plt.figure(figsize=(14, 6))

   #plot 4830 - 5040
   plt.subplot(1, 2, 1)
   plt.plot(wavelength, alphas[idx1], c='k', drawstyle='steps', label=label1)
   plt.plot(wavelength, alpha_stds[idx1], c='k', drawstyle='steps')
   plt.plot(wavelength, alphas[idx2], c='r', drawstyle='steps', label=label2)
   plt.plot(wavelength, alpha_stds[idx2], c='r', drawstyle='steps')

   if bootstrap:
       plt.fill_between(wavelength, bootstrap_lower[idx1], bootstrap_upper[idx1], linewidth=0.0, color='k', alpha=0.5, step='pre')
       plt.plot(wavelength, bootstrap_stds[idx1], c='m', drawstyle='steps')
       plt.fill_between(wavelength, bootstrap_lower[idx2], bootstrap_upper[idx2], linewidth=0.0, color='r', alpha=0.5, step='pre')
       plt.plot(wavelength, bootstrap_stds[idx2], c='m', drawstyle='steps')

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
   plt.plot(wavelength, alphas[idx1], c='k', drawstyle='steps', label=label1)
   plt.plot(wavelength, alpha_stds[idx1], c='k', drawstyle='steps')
   plt.plot(wavelength, alphas[idx2], c='r', drawstyle='steps', label=label2)
   plt.plot(wavelength, alpha_stds[idx2], c='r', drawstyle='steps')

   if bootstrap:
       plt.fill_between(wavelength, bootstrap_lower[idx1], bootstrap_upper[idx1], linewidth=0.0, color='k', alpha=0.5, step='pre')
       plt.plot(wavelength, bootstrap_stds[idx1], c='m', drawstyle='steps')
       plt.fill_between(wavelength, bootstrap_lower[idx2], bootstrap_upper[idx2], linewidth=0.0, color='r', alpha=0.5, step='pre')
       plt.plot(wavelength, bootstrap_stds[idx2], c='m', drawstyle='steps')
   
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
plot_emissions(0, 1, "SFD", r"With $\tau$ Correction")
plt.savefig('test.pdf')
plt.show()
plot_emissions(0, 3, "SFD", "With IRIS data")
plt.show()
plot_emissions(0, 2, "SFD", r"With $\tau$ and IRIS" )
plt.show()
       
def generate_binned_alphas(alphas, alpha_stds, wavelength_all, wavelength=None):
    #plot binned alpha vs wavelength (original)
    #wavelength_all is the one that determines the binned lambdas

    if wavelength is None:
        wavelength = wavelength_all
    
    lambda_range = wavelength_all[-1] - wavelength_all[0]
    left_over = lambda_range - 50*floor(lambda_range / 50)  
    binned_lambdas = np.arange(wavelength_all[1]+left_over/2, wavelength_all[-1], 50) #[1] to avoid going over left edge
    binned_alphas = []
    binned_stds = []

    #mask emission lines
    emission_line_mask = np.zeros(len(wavelength), dtype=int)
    emission_lines = [3727, 4863, 4960, 5008, 5877, 6550, 6565, 6585, 6718, 6733]
    for line in emission_lines:
        peak_idx = np.argmin(np.abs(wavelength-line))
        emission_line_mask[peak_idx-2:peak_idx+3] = 1
    #if boss, 3727
    
    for i in range(len(alphas)):

        binned_alpha_arr = np.zeros(binned_lambdas.shape)
        binned_std_arr = np.zeros(binned_lambdas.shape)
        for j, lmda in enumerate(binned_lambdas):
            indices = np.where((wavelength > lmda-50) & (wavelength < lmda) & np.logical_not(emission_line_mask))[0] #test
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
if not boss: #calculate binned spectrum for 1d boss
    binned_lambdas_boss, binned_alphas_boss, binned_stds_boss = generate_binned_alphas([alphas[-1]], [alpha_stds[-1]], wavelength, wavelength_boss)

if bootstrap:
    _, bootstrap_binned_lower, bootstrap_binned_stds = generate_binned_alphas(bootstrap_lower, bootstrap_stds, wavelength)
    _, bootstrap_binned_upper, _ = generate_binned_alphas(bootstrap_upper, bootstrap_stds, wavelength)
    if not boss:
        #first arg is arvitrary for binned_stds.
        _, _, bootstrap_binned_stds_boss = generate_binned_alphas([alphas[-1]], [bootstrap_stds[-1]], wavelength, wavelength_boss)
        _, bootstrap_binned_lower_boss, _ = generate_binned_alphas([bootstrap_lower_boss[0]], [bootstrap_stds[-1]], wavelength, wavelength_boss)
        _, bootstrap_binned_upper_boss, _ = generate_binned_alphas([bootstrap_upper_boss[0]], [bootstrap_stds[-1]], wavelength, wavelength_boss)
        bootstrap_binned_lower_boss = bootstrap_binned_lower_boss[0]
        bootstrap_binned_upper_boss = bootstrap_binned_upper_boss[0]


if boss:
    y_max = 0.55
    x_min = 3700
    x_max = 10000
else:
    y_max = 0.3
    x_min = 3850
    x_max = 9200


#2-panel plot, BOSS compared to SDSS (both 1d, etc.)

if not boss and bootstrap:
    plt.figure(figsize=(14, 6))

    plt.subplot(1, 2, 1)

    plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
    plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
    plt.plot(binned_lambdas, bootstrap_binned_stds[0], c='m', drawstyle='steps')
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[0], bootstrap_binned_upper[0], linewidth=0.0, color='k', alpha=0.2, step='pre')

    plt.legend(frameon=False)
    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"$\alpha_\lambda$")
    plt.xlim(x_min, x_max)
    plt.ylim(0, y_max)

    plt.subplot(1, 2, 2)

    plt.plot(binned_lambdas_boss, binned_alphas_boss[0], c='k', drawstyle='steps', label='BOSS')
    plt.plot(binned_lambdas_boss, bootstrap_binned_stds_boss[0], c='k', drawstyle='steps')
    plt.plot(binned_lambdas_boss, binned_stds_boss[0], c='m', drawstyle='steps')
    plt.fill_between(binned_lambdas, bootstrap_binned_lower_boss, bootstrap_binned_upper_boss, linewidth=0.0, color='k', alpha=0.2, step='pre')

    temp = y_max
    y_max = 0.3
    plt.legend(frameon=False)
    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"$\alpha_\lambda$")
    plt.xlim(x_min, x_max)
    plt.ylim(0, y_max)

    plt.tight_layout()
    plt.show()
    y_max = temp


    

    
'''
#vertical plot: 3 mods
    
plt.figure(figsize=(6, 14))
    
plt.subplot(3, 1, 1)    
#compare original to tao
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[1], c='r', drawstyle='steps', label=r'With $\tau$ Correction')
plt.plot(binned_lambdas, binned_stds[1], c='r', drawstyle='steps')
if bootstrap:
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[0], bootstrap_binned_upper[0], linewidth=0.0, color='k', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[0], c='m', drawstyle='steps')
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[1], bootstrap_binned_upper[1], linewidth=0.0, color='r', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[1], c='m', drawstyle='steps')
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
if bootstrap:
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[0], bootstrap_binned_upper[0], linewidth=0.0, color='k', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[0], c='m', drawstyle='steps')
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[3], bootstrap_binned_upper[3], linewidth=0.0, color='r', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[3], c='m', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.subplot(3, 1, 3)
#compare original to tao AND iris
#(or for SDSS, compare to BOSS)
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[2], c='r', drawstyle='steps', label=r'With IRIS and $\tau$')
plt.plot(binned_lambdas, binned_stds[2], c='r', drawstyle='steps')
if bootstrap:
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[0], bootstrap_binned_upper[0], linewidth=0.0, color='k', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[0], c='m', drawstyle='steps')
    plt.fill_between(binned_lambdas, bootstrap_binned_lower[2], bootstrap_binned_upper[2], linewidth=0.0, color='r', alpha=0.5, step='pre')
    plt.plot(binned_lambdas, bootstrap_binned_stds[2], c='m', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)

plt.tight_layout()
if save:
    plt.savefig("../bootstrap_sdss_binned_102919.png")
plt.show()
'''


'''
#plot original with all 3 modifications, all on one plot
#new model, IRIS, and BOSS

plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_alphas[1], c='r', drawstyle='steps', label=r'$\tau$ Correction')
plt.plot(binned_lambdas, binned_alphas[3], c='b', drawstyle='steps', label='IRIS')

if bootstrap:
    plt.plot(binned_lambdas, bootstrap_binned_stds[0], c='k', drawstyle='steps')
    plt.plot(binned_lambdas, bootstrap_binned_stds[1], c='r', drawstyle='steps')
    plt.plot(binned_lambdas, bootstrap_binned_stds[3], c='b', drawstyle='steps')
else:
    plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
    plt.plot(binned_lambdas, binned_stds[1], c='r', drawstyle='steps')
    plt.plot(binned_lambdas, binned_stds[3], c='b', drawstyle='steps')

temp = y_max
y_max = 0.35     
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(x_min, x_max)
plt.ylim(0, y_max)
plt.legend(frameon=False)
plt.show()
y_max = temp
'''


#threshold plots:

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
    #THIS CHUNK IS TEMP
    #plotting order: index: 2, 4, 5, 6
    alphas_thresh_1d = [np.load('../alphas_and_stds/alphas_boss_iris_1d_111119_3.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_111119_4.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_5.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_111119_6.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_111119_8.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_10.npy')]
    #also temp:
    alphas_thresh_1d = [np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_10.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_15.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_111119_175.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_92719_20.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_111119_225.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_92719_25.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_111119_275.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_92719_30.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_1d_111119_325.npy')]
    alpha_stds_thresh_1d = [np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_5.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_75.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_10.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_125.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_15.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_20.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_25.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_30.npy')]
    #TEMP:
    alpha_stds_thresh_1d = [np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_111119_3.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_111119_4.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_5.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_111119_6.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_111119_8.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_10.npy')]
    #also temp:
    alpha_stds_thresh_1d = [np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_10.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_15.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_111119_175.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_20.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_111119_225.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_25.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_111119_275.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_92719_30.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_111119_325.npy')]
    #alphas with various thresholds (BOSS, 2d, IRIS)
    alphas_thresh_2d = [np.load('../alphas_and_stds/alphas_boss_iris_91119_5.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_91119_75.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_91119_10.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_91119_125.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_91119_15.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_92719_20.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_92719_25.npy'), \
                     np.load('../alphas_and_stds/alphas_boss_iris_92719_30.npy')]
    #TEMP:
    alphas_thresh_2d = [np.load('../alphas_and_stds/alphas_boss_iris_111119_3.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_111119_4.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_91119_5.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_111119_6.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_111119_8.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_91119_10.npy')]
    #also TEMP:
    alphas_thresh_2d = [np.load('../alphas_and_stds/alphas_boss_iris_91119_10.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_91119_15.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_111119_175.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_92719_20.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_111119_225.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_92719_25.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_111119_275.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_92719_30.npy'), \
                        np.load('../alphas_and_stds/alphas_boss_iris_111119_325.npy')]
    alpha_stds_thresh_2d = [np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_5.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_75.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_10.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_125.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_15.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_20.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_25.npy'), \
                         np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_30.npy')]
    #TEMP:
    alpha_stds_thresh_2d = [np.load('../alphas_and_stds/alpha_stds_boss_iris_111119_3.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_111119_4.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_5.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_111119_6.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_111119_8.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_10.npy')]
    #also TEMP:
    alpha_stds_thresh_2d = [np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_10.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_15.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_111119_175.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_20.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_111119_225.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_25.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_111119_275.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_92719_30.npy'), \
                            np.load('../alphas_and_stds/alpha_stds_boss_iris_111119_325.npy')]
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
    x_min = 3850
    x_max = 9200
    y_max = .41

ax = fig.add_subplot(121)
plt.text(0.02, 0.98, 'Original\nModel', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10, fontweight='bold')

'''
colors = ['k', 'b', 'g', 'r', 'm', 'y']
labels = ['3', '4', '5', '6', '8', '10']
for i in range(6):
    plt.plot(binned_lambdas, binned_alphas_1d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i])
    plt.plot(binned_lambdas, binned_stds_1d[i], c=colors[i], drawstyle='steps')
    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"$\alpha_\lambda$")
    plt.xlim(x_min, x_max)
    plt.ylim(0, y_max)
'''

colors = ['k', 'b', 'g', 'r', 'm', 'y', 'brown', 'cyan', 'pink']
labels = ['10', '15', '17.5', '20', '22.5', '25', '27.5', '30', '32.5']
for i in range(9):
    plt.plot(binned_lambdas, binned_alphas_1d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i])
    plt.plot(binned_lambdas, binned_stds_1d[i], c=colors[i], drawstyle='steps')
    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"$\alpha_\lambda$")
    plt.xlim(x_min, x_max)
    plt.ylim(0, y_max)

if boss:
    leg = plt.legend(frameon=False, loc='upper right')
    plt.setp(leg.texts, family='monospace')
else:
    leg = plt.legend(frameon=False, loc='lower center')
    plt.setp(leg.texts, family='monospace')

plt.text(0, 100, 'Original\nModel')

ax = fig.add_subplot(122)
plt.text(0.02, 0.98, 'Tao\nModel', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10, fontweight='bold')

'''
colors = ['k', 'b', 'g', 'r', 'm', 'y']
labels = ['3', '4', '5', '6', '8', '10']
for i in range(6):
    plt.plot(binned_lambdas, binned_alphas_2d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i])
    plt.plot(binned_lambdas, binned_stds_2d[i], c=colors[i], drawstyle='steps')
    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"$\alpha_\lambda$")
    plt.xlim(x_min, x_max)
    plt.ylim(0, y_max)
'''

colors = ['k', 'b', 'g', 'r', 'm', 'y', 'brown', 'cyan', 'pink']
labels = ['10', '15', '17.5', '20', '22.5', '25', '27.5', '30', '32.5']
for i in range(9):
    plt.plot(binned_lambdas, binned_alphas_2d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i])
    plt.plot(binned_lambdas, binned_stds_2d[i], c=colors[i], drawstyle='steps')
    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"$\alpha_\lambda$")
    plt.xlim(x_min, x_max)
    plt.ylim(0, y_max)

if boss:
    leg = plt.legend(frameon=False, loc='upper right')
    plt.setp(leg.texts, family='monospace')
else:
    leg = plt.legend(frameon=False, loc='lower center')
    plt.setp(leg.texts, family='monospace')
    

plt.text(0, 100, 'Tao Model')

plt.tight_layout()

if save_thresh:
    if boss:
        plt.savefig("../boss_thresholds_2panel_110919.png")
    else:
        plt.savefig("../sdss_thresholds_2panel_110919.png")
plt.show()


