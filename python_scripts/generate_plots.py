# run this after generate_alphas_ahab.py. If want to use bootstrapping, also have to run bootstrap_ahab.py first.
# generate plots of correlation spectra (overall and certain sections)
# need to run 2x to cover all figures for the paper:
#   boss=0 and bootstrap=1
#   boss=1 and bootstrap=1
# flux conversion factor is handled here, not in reproduce_figs

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from astropy.io import fits
import sys #for command line args
from math import floor #for binning range

'''
COMMAND LINE OPTIONS:

boss: 1 to use boss data, 0 to use sdss-II data
save: enter savekey. figures for the paper will get filenames based on the key. otherwise figures are not saved.
bootstrap: 1 to plot with bootstrap errors (and sometimes envelopes)
loadkey: alphas will be loaded based on this key
'''

#command line options
if len(sys.argv) != 5:
    print("Usage: generate_plots.py [boss: 0, 1] [save: 0, savekey] [bootstrap: 0, 1] [loadkey]")
    exit(0)
boss = int(sys.argv[1])
save = sys.argv[2]
bootstrap = int(sys.argv[3])
loadkey = sys.argv[4]

# global plotting options
matplotlib.rcParams['axes.labelsize'] = 'large'
matplotlib.rcParams['xtick.labelsize'] = 'large'
matplotlib.rcParams['ytick.labelsize'] = 'large'
matplotlib.rcParams['legend.fontsize'] = 'large'
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['ytick.left'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['ytick.major.size'] = 6.0
matplotlib.rcParams['ytick.minor.size'] = 3.0
matplotlib.rcParams['xtick.minor.size'] = 3.0
matplotlib.rcParams['xtick.major.size'] = 6.0
#matplotlib.rcParams['font.size'] = 'large'

ahab = True

if ahab:
    alpha_direc = '../alphas_and_stds/'
    alpha_direc_boot = '../data/'  #'../alphas_and_stds'
    hdulist_direc = '../data/'  #'/Users/blakechellew/Documents/DustProject/BrandtFiles/'
else:
    alpha_direc = '../alphas_and_stds/'
    alpha_direc_boot = '../alphas_and_stds'
    hdulist_direc = '/Users/blakechellew/Documents/DustProject/BrandtFiles/'

# load wavelengths
wavelength_boss = np.load('../alphas_and_stds/wavelength_boss.npy')
hdulist = fits.open(hdulist_direc + 'SDSS_allskyspec.fits')
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
alphas_boss = [np.load(alpha_direc + 'alphas_boss_iris_2d_north_' + loadkey + '.npy'), \
               np.load(alpha_direc + 'alphas_boss_iris_2d_south_' + loadkey + '.npy'), \
               np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_10.npy')]
alpha_stds_boss = [np.load(alpha_direc + 'alpha_stds_boss_iris_2d_north_' + loadkey + '.npy'), \
                   np.load(alpha_direc + 'alpha_stds_boss_iris_2d_south_' + loadkey + '.npy'), \
                   np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')]

#sdss alphas
alphas_sdss = [np.load(alpha_direc + 'alphas_sdss_1d_' + loadkey + '.npy'), \
               np.load(alpha_direc + 'alphas_sdss_2d_' + loadkey + '.npy'), \
               np.load(alpha_direc + 'alphas_sdss_iris_2d_' + loadkey + '.npy'), \
               np.load(alpha_direc + 'alphas_sdss_iris_1d_' + loadkey + '.npy'), \
               np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_10.npy')] #don't change the last element
alpha_stds_sdss = [np.load(alpha_direc + 'alpha_stds_sdss_1d_' + loadkey + '.npy'), \
                   np.load(alpha_direc + 'alpha_stds_sdss_2d_' + loadkey + '.npy'), \
                   np.load(alpha_direc + 'alpha_stds_sdss_iris_2d_' + loadkey + '.npy'), \
                   np.load(alpha_direc + 'alpha_stds_sdss_iris_1d_' + loadkey + '.npy'), \
                   np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')]

if bootstrap:
    bootstrap_alphas_boss = [np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_north_' + loadkey + '.npy'), \
                             np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_south_' + loadkey + '.npy'), \
                             np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_10.npy')]
    bootstrap_alpha_stds_boss = [np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_north_' + loadkey + '.npy'), \
                                 np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_south_' + loadkey + '.npy'), \
                                 np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')]

    bootstrap_alphas_sdss = [np.load(alpha_direc_boot + 'bootstrap_alphas_sdss_1d_' + loadkey + '.npy'), \
                             np.load(alpha_direc_boot + 'bootstrap_alphas_sdss_2d_' + loadkey + '.npy'), \
                             np.load(alpha_direc_boot + 'bootstrap_alphas_sdss_iris_2d_' + loadkey + '.npy'), \
                             np.load(alpha_direc_boot + 'bootstrap_alphas_sdss_iris_1d_' + loadkey + '.npy'), \
                             np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_10.npy')]
    bootstrap_alpha_stds_sdss = [np.load(alpha_direc_boot + 'bootstrap_alpha_stds_sdss_1d_' + loadkey + '.npy'), \
                                 np.load(alpha_direc_boot + 'bootstrap_alpha_stds_sdss_2d_' + loadkey + '.npy'), \
                                 np.load(alpha_direc_boot + 'bootstrap_alpha_stds_sdss_iris_2d_' + loadkey + '.npy'), \
                                 np.load(alpha_direc_boot + 'bootstrap_alpha_stds_sdss_iris_1d_' + loadkey + '.npy'), \
                                 np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')]
    
    '''
    #temporary: truncate bootstrap arrays to make it faster
    num_samples = 1000
    bootstrap_alphas_boss = [b[:num_samples, :] for b in bootstrap_alphas_boss] #temp
    bootstrap_alpha_stds_boss = [b[:num_samples, :] for b in bootstrap_alpha_stds_boss] #temp
    bootstrap_alphas_sdss = [b[:num_samples, :] for b in bootstrap_alphas_sdss] #temp
    bootstrap_alpha_stds_sdss = [b[:num_samples, :] for b in bootstrap_alpha_stds_sdss] #temp
    '''
    
    if boss:
        bootstrap_lower = [np.nanpercentile(b, 16, axis=0) / boss_fluxfactor for b in bootstrap_alphas_boss]
        bootstrap_upper = [np.nanpercentile(b, 84, axis=0) / boss_fluxfactor for b in bootstrap_alphas_boss]
        bootstrap_stds = [(bootstrap_upper[i] - bootstrap_lower[i])/2 for i in range(len(bootstrap_lower))]
    else:
        bootstrap_lower = [np.nanpercentile(b, 16, axis=0) / sdss_fluxfactor for b in bootstrap_alphas_sdss]
        bootstrap_upper = [np.nanpercentile(b, 84, axis=0) / sdss_fluxfactor for b in bootstrap_alphas_sdss]
        bootstrap_stds = [(bootstrap_upper[i] - bootstrap_lower[i])/2 for i in range(len(bootstrap_lower))]
    
    
#flux conversion factor:
alphas_sdss = [a/sdss_fluxfactor for a in alphas_sdss]
alphas_boss = [a/boss_fluxfactor for a in alphas_boss]
alpha_stds_sdss = [a/sdss_fluxfactor for a in alpha_stds_sdss]
alpha_stds_boss = [a/boss_fluxfactor for a in alpha_stds_boss]

if boss:
    alphas = alphas_boss
    alpha_stds = alpha_stds_boss
    if bootstrap:
        bootstrap_alphas = bootstrap_alphas_boss
        bootstrap_alpha_stds = bootstrap_alpha_stds_boss
else:
    alphas = alphas_sdss
    alpha_stds = alpha_stds_sdss
    if bootstrap:
        bootstrap_alphas = bootstrap_alphas_sdss
        bootstrap_alpha_stds = bootstrap_alpha_stds_sdss
        
num_arrays = len(alphas)

#plot unbinned spectra (wavelength ranges: 4830-5040 and 6530-6770)
def plot_emissions(alpha_indices, labels, colors):
    plt.figure(figsize=(12, 5))
   
    #plot 4830 - 5040
    ax1 = plt.subplot(1, 2, 1)
    for i, idx in enumerate(alpha_indices):
        ax1.plot(wavelength, alphas[idx], c=colors[i], drawstyle='steps', label=labels[i])
        if bootstrap:
            #ax1.fill_between(wavelength, bootstrap_lower[idx], bootstrap_upper[idx], linewidth=0.0, color=colors[i], alpha=0.5, step='pre')
            ax1.plot(wavelength, bootstrap_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        else:
            ax1.plot(wavelength, alpha_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')

    ax1.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    ax1.set_ylabel(r"$\alpha_\lambda$")
    ax1.legend(loc='upper center', frameon=False)
    ax1.set_xlim(4830, 5040)
    ax1.set_ylim(0, 0.99)
    xcoords = [4863, 4960, 5008]
    for xc in xcoords:
        ax1.axvline(x=xc, color='k', linewidth=1, linestyle='--')
    ax1.text(4853, 0.4, r"H$\beta$")
    ax1.text(4943, 0.4, "O[III]")
    ax1.text(4991, 0.4, "O[III]")

    ax1.xaxis.set_major_locator(MultipleLocator(50))
    ax1.xaxis.set_minor_locator(MultipleLocator(10))
    ax1.yaxis.set_major_locator(MultipleLocator(0.2))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.04))

    #line from 03 continuum::
    #ax1.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
    #actual continuum for NII:
    #ax1.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')

    #plot 6530 - 6770 (original vs tao)
    ax2 = plt.subplot(1, 2, 2)
    for i, idx in enumerate(alpha_indices):
        ax2.plot(wavelength, alphas[idx], c=colors[i], drawstyle='steps', label=labels[i])
        if bootstrap:
            #ax2.fill_between(wavelength, bootstrap_lower[idx], bootstrap_upper[idx], linewidth=0.0, color=colors[i], alpha=0.5, step='pre')
            ax2.plot(wavelength, bootstrap_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        else:
            ax2.plot(wavelength, alpha_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')

    ax2.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    ax2.set_ylabel(r"$\alpha_\lambda$")
    ax2.legend(loc='upper center', frameon=False)
    ax2.set_xlim(6530, 6770)
    ax2.set_ylim(0, 0.99)
    xcoords = [6550, 6565, 6585, 6718, 6733]
    for xc in xcoords:
        ax2.axvline(x=xc, color='k', linewidth=1, linestyle='--')
    ax2.text(6535, 0.4, "N[II]")
    ax2.text(6568, 0.9, r"H$\alpha$")
    ax2.text(6590, 0.6, "N[II]")
    ax2.text(6700, 0.6, "S[II]")
    ax2.text(6738, 0.5, "S[II]")

    ax2.xaxis.set_major_locator(MultipleLocator(50))
    ax2.xaxis.set_minor_locator(MultipleLocator(10))
    ax2.yaxis.set_major_locator(MultipleLocator(0.2))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.04))

    #line from 03 continuum::
    #ax2.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
    #actual continuum for NII:
    #ax2.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')


if boss:
    plot_emissions([0, 1, 2], ["North", "South", "Full Sky"], ['#004488', '#BB5566', '#DDAA33'])
    if save != '0' and bootstrap:
        plt.savefig('../paper_figures/unbinned_' + loadkey + '.pdf', bbox_inches='tight')
        plt.clf()
    else:
        plt.show()
        

'''
#unbinned plots, SFD 1d vs 3 others
plot_emissions([0, 1], ["SFD", r"With $\tau$ Correction"], ['k', 'r'])
plt.show()
plot_emissions([0, 3], ["SFD", "With IRIS data"], ['k', 'r'])
plt.show()
plot_emissions([0, 2], ["SFD", r"With $\tau$ and IRIS"], ['k', 'r'])
plt.show()
'''
       
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
    emission_lines = [4863, 4960, 5008, 5877, 6550, 6565, 6585, 6718, 6733]
    if boss:
        emission_lines.insert(0, 3727)
    for line in emission_lines:
        peak_idx = np.argmin(np.abs(wavelength-line))
        emission_line_mask[peak_idx-3:peak_idx+4] = 1
    
    for i in range(len(alphas)):

        binned_alpha_arr = np.zeros(binned_lambdas.shape)
        binned_std_arr = np.zeros(binned_lambdas.shape)
        if alphas[i].ndim > 1:
            binned_alpha_arr = np.zeros((alphas[i].shape[0], binned_lambdas.shape[0]))
            binned_std_arr = np.zeros((alphas[i].shape[0], binned_lambdas.shape[0]))
        for j, lmda in enumerate(binned_lambdas):
            indices = np.where((wavelength > lmda-50) & (wavelength < lmda) & np.logical_not(emission_line_mask))[0] #test
            if alphas[i].ndim > 1:
                relevant_alphas = alphas[i][:,indices]
                relevant_stds = alpha_stds[i][:,indices]
            else:
                relevant_alphas = alphas[i][indices]
                relevant_stds = alpha_stds[i][indices]
                
            #weighted average:
            variance = np.power(relevant_stds, 2)
            if alphas[i].ndim > 1:
                numerator = np.sum(np.divide(relevant_alphas, variance), axis=1)
                denominator = np.sum(np.divide(1, variance), axis=1)
            else:
                numerator = np.sum(np.divide(relevant_alphas, variance))
                denominator = np.sum(np.divide(1, variance))
            avg1 = numerator / denominator
            avg2 = 1 / denominator
            if alphas[i].ndim > 1:
                binned_alpha_arr[:,j] = avg1
                binned_std_arr[:,j] = np.sqrt(avg2)
            else:
                binned_alpha_arr[j] = avg1
                binned_std_arr[j] = np.sqrt(avg2)
                
        binned_alphas.append(binned_alpha_arr)
        binned_stds.append(binned_std_arr)

    return binned_lambdas, binned_alphas, binned_stds


binned_lambdas, binned_alphas, binned_stds = generate_binned_alphas(alphas, alpha_stds, wavelength)
if not boss: #calculate binned spectrum for 1d boss
    binned_lambdas_boss, binned_alphas_boss, binned_stds_boss = generate_binned_alphas([alphas[-1]], [alpha_stds[-1]], wavelength, wavelength_boss)

    
if bootstrap:
    #bin all the bootstrap spectra
    #use one spectrum to get binned wavelength
    _, bootstrap_binned_alphas, _ = generate_binned_alphas(bootstrap_alphas, bootstrap_alpha_stds, wavelength)

    #now look at percentiles:
    bootstrap_binned_lower = [np.nanpercentile(b, 16, axis=0) / boss_fluxfactor for b in bootstrap_binned_alphas] #68 percent confidence interval
    bootstrap_binned_upper = [np.nanpercentile(b, 84, axis=0) / boss_fluxfactor for b in bootstrap_binned_alphas]
    bootstrap_binned_stds = [(bootstrap_binned_upper[i] - bootstrap_binned_lower[i])/2 for i in range(len(bootstrap_binned_lower))]

    if not boss:
        _, bootstrap_binned_alphas_boss, _ = generate_binned_alphas([bootstrap_alphas[-1]], [bootstrap_alpha_stds[-1]], wavelength, wavelength_boss)
        bootstrap_binned_lower_boss = [np.nanpercentile(b, 16, axis=0) / boss_fluxfactor for b in bootstrap_binned_alphas_boss] #68 percent confidence interval
        bootstrap_binned_upper_boss = [np.nanpercentile(b, 84, axis=0) / boss_fluxfactor for b in bootstrap_binned_alphas_boss]
        bootstrap_binned_stds_boss = [(bootstrap_binned_upper_boss[i] - bootstrap_binned_lower_boss[i])/2 for i in range(len(bootstrap_binned_lower_boss))]


if boss:
    y_max = 0.3
    x_min = 3700
    x_max = 10000
else:
    y_max = 0.3
    x_min = 3850
    x_max = 9200

#2-panel plot, BOSS compared to SDSS (both 1d, etc.)
if not boss and bootstrap:
    plt.figure(figsize=(12, 5))

    ax1 = plt.subplot(1, 2, 1)

    ax1.plot(binned_lambdas, binned_alphas[2], c='k', drawstyle='steps', label='SDSS')
    ax1.plot(binned_lambdas, binned_stds[2], c='k', drawstyle='steps', linestyle='--')
    ax1.plot(binned_lambdas, bootstrap_binned_stds[2], c='m', drawstyle='steps', linestyle='--')
    ax1.fill_between(binned_lambdas, bootstrap_binned_lower[2], bootstrap_binned_upper[2], linewidth=0.0, color='k', alpha=0.2, step='pre')

    ax1.legend(frameon=False)
    ax1.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    ax1.set_ylabel(r"$\alpha_\lambda$")
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(0, y_max)

    ax1.xaxis.set_major_locator(MultipleLocator(1000))
    ax1.xaxis.set_minor_locator(MultipleLocator(200))
    ax1.yaxis.set_major_locator(MultipleLocator(0.05))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.01))

    ax2 = plt.subplot(1, 2, 2)

    ax2.plot(binned_lambdas_boss, binned_alphas_boss[0], c='k', drawstyle='steps', label='BOSS')
    ax2.plot(binned_lambdas_boss, binned_stds_boss[0], c='k', drawstyle='steps', linestyle='--')
    ax2.plot(binned_lambdas_boss, bootstrap_binned_stds_boss[0], c='m', drawstyle='steps', linestyle='--')
    ax2.fill_between(binned_lambdas_boss, bootstrap_binned_lower_boss[0], bootstrap_binned_upper_boss[0], linewidth=0.0, color='k', alpha=0.2, step='pre')

    temp = y_max
    y_max = 0.29
    ax2.legend(frameon=False)
    ax2.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    ax2.set_ylabel(r"$\alpha_\lambda$")
    ax2.set_xlim(x_min, x_max)
    ax2.set_ylim(0, y_max)

    ax2.xaxis.set_major_locator(MultipleLocator(1000))
    ax2.xaxis.set_minor_locator(MultipleLocator(200))
    ax2.yaxis.set_major_locator(MultipleLocator(0.05))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.01))

    # plt.tight_layout()

    if save != '0':
        plt.savefig('../paper_figures/compare_boss_sdss_' + loadkey + '.pdf', bbox_inches='tight')
        plt.clf()
    else:
        plt.show()
    y_max = temp

#plot binned alphas
#takes alphas already binned. use generate_binned_alphas
def plot_binned(alpha_indices, colors, labels, envelope=False):
    fig, ax = plt.subplots()

    for i, idx in enumerate(alpha_indices):
        ax.plot(binned_lambdas, binned_alphas[idx], c=colors[i], drawstyle='steps', label=labels[i])
        if bootstrap:
            ax.plot(binned_lambdas, bootstrap_binned_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        else:
            ax.plot(binned_lambdas, binned_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        if envelope:
            ax.fill_between(binned_lambdas, bootstrap_binned_lower[idx], bootstrap_binned_upper[idx], linewidth=0.0, color=colors[i], alpha=0.2, step='pre')
            
    y_max = 0.35
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(200))
    ax.yaxis.set_minor_locator(MultipleLocator(0.02))
    #ax.figure(figsize=(6, 5))
    ax.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    ax.set_ylabel(r"$\alpha_\lambda$")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, y_max)
    ax.legend(frameon=False)
    y_max = 0.3


if boss and bootstrap:
    #calculate difference between average north and average south
    #also compare a color index
    #using bootstrapping to get the std devs

    #color index:
    wavelength_deltas = np.array([wavelength[i+1] - wavelength[i] for i in range(len(wavelength)-1)])
    wavelength_deltas = np.append(wavelength_deltas, wavelength_deltas[-1])
    idx_4000 = np.argmin(np.abs(wavelength-4000))
    idx_5000 = np.argmin(np.abs(wavelength-5000))
    idx_6000 = np.argmin(np.abs(wavelength-6000))
    idx_7000 = np.argmin(np.abs(wavelength-7000))
    idx_8000 = np.argmin(np.abs(wavelength-8000))
    idx_9000 = np.argmin(np.abs(wavelength-9000))
    
    #north:
    variance = np.power(bootstrap_stds[0], 2)
    rel_alphas = alphas[0]
    numerator = np.nansum(np.divide(rel_alphas, variance))
    denominator = np.nansum(np.divide(1, variance))
    avg_north = numerator / denominator
    avg_north_var = np.nanvar(avg_north)
    avg_north = np.nanmean(avg_north)

    integrand = np.multiply(rel_alphas, wavelength_deltas)
    integrand_std = np.multiply(variance, wavelength_deltas**2)
    integral_4_5 = np.nansum(integrand[idx_4000:idx_5000])
    integral_6_7 = np.nansum(integrand[idx_6000:idx_7000])
    integral_8_9 = np.nansum(integrand[idx_8000:idx_9000])
    integral_4_5_std = np.sqrt(np.nansum(integrand_std[idx_4000:idx_5000]))
    integral_6_7_std = np.sqrt(np.nansum(integrand_std[idx_6000:idx_7000]))
    integral_8_9_std = np.sqrt(np.nansum(integrand_std[idx_8000:idx_9000]))
    north_color_idx = integral_8_9 - integral_4_5
    north_color_idx_avg = np.mean(north_color_idx)
    north_color_idx_var = np.nanvar(north_color_idx)
    north_ere = 2 * integral_6_7 / (integral_4_5 + integral_8_9)
    north_ere_frac_1 = np.sqrt(integral_4_5_std**2 + integral_8_9_std**2)/(integral_4_5 + integral_8_9)
    north_ere_frac_2 = integral_6_7_std / integral_6_7
    north_ere_std = 2 * np.sqrt(north_ere_frac_1**2 + north_ere_frac_2**2)

    #south:
    variance = np.power(bootstrap_stds[1], 2)
    rel_alphas = alphas[1]
    numerator = np.nansum(np.divide(rel_alphas, variance))
    denominator = np.nansum(np.divide(1, variance))
    avg_south = numerator / denominator
    avg_south_var = np.nanvar(avg_south)
    avg_south = np.nanmean(avg_south)

    integrand = np.multiply(rel_alphas, wavelength_deltas)
    integrand_std = np.multiply(variance, wavelength_deltas ** 2)
    integral_4_5 = np.nansum(integrand[idx_4000:idx_5000])
    integral_6_7 = np.nansum(integrand[idx_6000:idx_7000])
    integral_8_9 = np.nansum(integrand[idx_8000:idx_9000])
    integral_4_5_std = np.sqrt(np.nansum(integrand_std[idx_4000:idx_5000]))
    integral_6_7_std = np.sqrt(np.nansum(integrand_std[idx_6000:idx_7000]))
    integral_8_9_std = np.sqrt(np.nansum(integrand_std[idx_8000:idx_9000]))
    south_color_idx = integral_8_9 - integral_4_5
    south_color_idx_avg = np.mean(south_color_idx)
    south_color_idx_var = np.nanvar(south_color_idx)
    south_ere = 2 * integral_6_7 / (integral_4_5 + integral_8_9)
    south_ere_frac_1 = np.sqrt(integral_4_5_std ** 2 + integral_8_9_std ** 2) / (integral_4_5 + integral_8_9)
    south_ere_frac_2 = integral_6_7_std / integral_6_7
    south_ere_std = 2 * np.sqrt(south_ere_frac_1 ** 2 + south_ere_frac_2 ** 2)
    
    #print results
    print("north avg:", avg_north)
    print("south avg:", avg_south)
    print("north idx:", north_color_idx_avg)
    print("south idx:", south_color_idx_avg)
    print("north ere:", north_ere, north_ere_std)
    print("south ere:", south_ere, south_ere_std)

    #calculate significance:
    #estimate variance:
    var_diff = avg_north_var + avg_south_var
    z = (avg_north - avg_south) / np.sqrt(var_diff)
    print("z value for avg:")
    print(z)

    var_diff = north_color_idx_var + south_color_idx_var
    z = (north_color_idx_avg - south_color_idx_avg) / np.sqrt(var_diff)
    print("z value for idx:")
    print(z)

    var_diff = north_ere_std**2 + south_ere_std**2
    z = (south_ere - north_ere) / np.sqrt(var_diff)
    print("z value for ere", z)
    
    
#plot original with 2 modifications, all on one plot
#new model, IRIS
if not boss:
    plot_binned([0, 1, 3], ['#004488', '#BB5566', '#DDAA33'], ['SFD', 'Nonlinear Model', 'IRIS'])
    if save != '0' and bootstrap:
        plt.savefig('../paper_figures/all_3_mods_' + loadkey + '.pdf', bbox_inches='tight')
        plt.clf()
    else:
        plt.show()

#plot north and south on same plot (boss)
if boss:
    envelope = bootstrap
    plot_binned([0, 1], ['#004488', '#BB5566'], ['North', 'South'], envelope=envelope)
    if save != '0' and bootstrap:
        plt.savefig('../paper_figures/boss_north_south_' + loadkey + '.pdf', bbox_inches='tight')
        plt.clf()
    else:
        plt.show()

#threshold plots:
if boss:
    #alphas with various thresholds (BOSS, 1d, IRIS)
    alphas_thresh_1d = [np.load(alpha_direc + 'alphas_boss_iris_1d_' + loadkey + '_10.npy'), \
                        np.load(alpha_direc + 'alphas_boss_iris_1d_' + loadkey + '_15.npy'), \
                        np.load(alpha_direc + 'alphas_boss_iris_1d_' + loadkey + '_20.npy'), \
                        np.load(alpha_direc + 'alphas_boss_iris_1d_' + loadkey + '_25.npy'), \
                        np.load(alpha_direc + 'alphas_boss_iris_1d_' + loadkey + '_30.npy')]
    alpha_stds_thresh_1d = [np.load(alpha_direc + 'alpha_stds_boss_iris_1d_' + loadkey + '_10.npy'), \
                            np.load(alpha_direc + 'alpha_stds_boss_iris_1d_' + loadkey + '_15.npy'), \
                            np.load(alpha_direc + 'alpha_stds_boss_iris_1d_' + loadkey + '_20.npy'), \
                            np.load(alpha_direc + 'alpha_stds_boss_iris_1d_' + loadkey + '_25.npy'), \
                            np.load(alpha_direc + 'alpha_stds_boss_iris_1d_' + loadkey + '_30.npy')]
    #alphas with various thresholds (BOSS, 2d, IRIS)
    alphas_thresh_2d = [np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_10.npy'), \
                        np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_15.npy'), \
                        np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_20.npy'), \
                        np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_25.npy'), \
                        np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_30.npy')]
    alpha_stds_thresh_2d = [np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_10.npy'), \
                            np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_15.npy'), \
                            np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_20.npy'), \
                            np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_25.npy'), \
                            np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_30.npy')]

    if bootstrap:
        #bootstrap alphas and stds
        bootstrap_alphas_thresh_1d = [np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_10.npy'), \
                                      np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_15.npy'), \
                                      np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_20.npy'), \
                                      np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_25.npy'), \
                                      np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_30.npy')]
        bootstrap_alpha_stds_thresh_1d = [np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_10.npy'), \
                                          np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_15.npy'), \
                                          np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_20.npy'), \
                                          np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_25.npy'), \
                                          np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_30.npy')]
        bootstrap_alphas_thresh_2d = [np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_10.npy'), \
                                      np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_15.npy'), \
                                      np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_20.npy'), \
                                      np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_25.npy'), \
                                      np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_30.npy')]
        bootstrap_alpha_stds_thresh_2d = [np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_10.npy'), \
                                          np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_15.npy'), \
                                          np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_20.npy'), \
                                          np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_25.npy'), \
                                          np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_30.npy')]

    #flux factor corrections:
    alphas_thresh_1d = [a/boss_fluxfactor for a in alphas_thresh_1d]
    alphas_stds_thresh_1d = [a/boss_fluxfactor for a in alpha_stds_thresh_1d]
    alphas_thresh_2d = [a/boss_fluxfactor for a in alphas_thresh_2d]
    alphas_stds_thresh_2d = [a/boss_fluxfactor for a in alpha_stds_thresh_2d]

    #binning
    binned_lambdas, binned_alphas_1d, binned_stds_1d = generate_binned_alphas(alphas_thresh_1d, alpha_stds_thresh_1d, wavelength)
    binned_lambdas, binned_alphas_2d, binned_stds_2d = generate_binned_alphas(alphas_thresh_2d, alpha_stds_thresh_2d, wavelength)
    if bootstrap:
        _, bootstrap_binned_alphas_thresh_1d, _ = generate_binned_alphas(bootstrap_alphas_thresh_1d, bootstrap_alpha_stds_thresh_1d, wavelength)
        _, bootstrap_binned_alphas_thresh_2d, _ = generate_binned_alphas(bootstrap_alphas_thresh_2d, bootstrap_alpha_stds_thresh_2d, wavelength)

    if bootstrap:
        #bootstrap percentiles
        bootstrap_binned_lower_thresh_1d = [np.nanpercentile(b, 16, axis=0) / boss_fluxfactor for b in bootstrap_binned_alphas_thresh_1d] #68 percent confidence interval
        bootstrap_binned_upper_thresh_1d = [np.nanpercentile(b, 84, axis=0) / boss_fluxfactor for b in bootstrap_binned_alphas_thresh_1d]
        bootstrap_binned_stds_thresh_1d = [(bootstrap_binned_upper_thresh_1d[i] - bootstrap_binned_lower_thresh_1d[i])/2 for i in range(len(bootstrap_binned_lower_thresh_1d))]
        bootstrap_binned_lower_thresh_2d = [np.nanpercentile(b, 16, axis=0) / boss_fluxfactor for b in bootstrap_binned_alphas_thresh_2d] #68 percent confidence interval
        bootstrap_binned_upper_thresh_2d = [np.nanpercentile(b, 84, axis=0) / boss_fluxfactor for b in bootstrap_binned_alphas_thresh_2d]
        bootstrap_binned_stds_thresh_2d = [(bootstrap_binned_upper_thresh_2d[i] - bootstrap_binned_lower_thresh_2d[i])/2 for i in range(len(bootstrap_binned_lower_thresh_2d))]
    
    fig = plt.figure(figsize=(12, 5), dpi=200)
    x_min = 3700
    x_max = 10100
    y_max = .3

    ax1 = fig.add_subplot(121)
    ax1.text(0.98, 0.98, 'Linear\nModel', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes, fontsize='large')

    #"bright" color scheme from https://personal.sron.nl/~pault/
    # (order: blue, red, green, yellow, cyan)
    colors = ['#4477AA', '#CCBB44', '#66CCEE', '#EE6677', '#228833']
    labels = ['10', '15', '20', '25', '30']
    for i in range(len(labels)):
        ax1.plot(binned_lambdas, binned_alphas_1d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i])
        if bootstrap:
            ax1.plot(binned_lambdas, bootstrap_binned_stds_thresh_1d[i], c = colors[i], drawstyle='steps', linestyle='--')
            #ax1.fill_between(binned_lambdas, bootstrap_binned_lower_thresh_1d[i], bootstrap_binned_upper_thresh_1d[i], linewidth=0.0, color=colors[i], alpha=0.5, step='pre') #TEMP
        else:
            ax1.plot(binned_lambdas, binned_stds_1d[i], c=colors[i], drawstyle='steps', linestyle='--')
        ax1.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
        ax1.set_ylabel(r"$\alpha_\lambda$")
        ax1.set_xlim(x_min, x_max)
        ax1.set_ylim(0, y_max)

    leg = ax1.legend(frameon=False, loc='upper left', ncol=2)

    ax1.xaxis.set_major_locator(MultipleLocator(1000))
    ax1.xaxis.set_minor_locator(MultipleLocator(200))
    ax1.yaxis.set_major_locator(MultipleLocator(0.05))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.01))

    ax2 = fig.add_subplot(122)
    ax2.text(0.98, 0.98, 'Nonlinear\nModel', horizontalalignment='right', verticalalignment='top', transform=ax2.transAxes, fontsize='large')

    for i in range(len(labels)):
        ax2.plot(binned_lambdas, binned_alphas_2d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i])
        if bootstrap:
            ax2.plot(binned_lambdas, bootstrap_binned_stds_thresh_2d[i], c = colors[i], drawstyle='steps', linestyle='--')
            #ax2.fill_between(binned_lambdas, bootstrap_binned_lower_thresh_2d[i], bootstrap_binned_upper_thresh_2d[i], linewidth=0.0, color=colors[i], alpha=0.5, step='pre') #TEMP
        else:
            ax2.plot(binned_lambdas, binned_stds_2d[i], c=colors[i], drawstyle='steps', linestyle='--')
        ax2.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
        ax2.set_ylabel(r"$\alpha_\lambda$")
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(0, y_max)

    leg = ax2.legend(frameon=False, loc='lower center', ncol=2)

    ax2.xaxis.set_major_locator(MultipleLocator(1000))
    ax2.xaxis.set_minor_locator(MultipleLocator(200))
    ax2.yaxis.set_major_locator(MultipleLocator(0.05))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.01))

    if save != '0' and bootstrap:
        plt.savefig('../paper_figures/boss_thresholds_2panel_' + loadkey + '.pdf', bbox_inches='tight')
        plt.clf()
    else:
        plt.show()

else:
    pass #don't want to plot SDSS thresholds right now


    '''
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

    binned_lambdas, binned_alphas_1d, binned_stds_1d = generate_binned_alphas(alphas_thresh_1d, alpha_stds_thresh_1d, wavelength)
    binned_lambdas, binned_alphas_2d, binned_stds_2d = generate_binned_alphas(alphas_thresh_2d, alpha_stds_thresh_2d, wavelength)
    fig = plt.figure(figsize=(8, 4), dpi=200)
    x_min = 3850
    x_max = 9200
    y_max = .41

    ax = fig.add_subplot(121)
    plt.text(0.02, 0.98, 'Original\nModel', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10, fontweight='bold')

    colors = ['k', 'b', 'g', 'r', 'm', 'y', 'brown', 'cyan']
    labels = ['5', '7.5', '10', '12.5', '15', '20', '25', '30']
    for i in range(len(labels)):
        plt.plot(binned_lambdas, binned_alphas_1d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i])
        plt.plot(binned_lambdas, binned_stds_1d[i], c=colors[i], drawstyle='steps')
        plt.xlabel(r"Wavelength ($\AA$)")
        plt.ylabel(r"$\alpha_\lambda$")
        plt.xlim(x_min, x_max)
        plt.ylim(0, y_max)

    leg = plt.legend(frameon=False, loc='lower center')
    plt.setp(leg.texts, family='monospace')

    plt.text(0, 100, 'Original\nModel')

    ax = fig.add_subplot(122)
    plt.text(0.02, 0.98, 'Tao\nModel', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10, fontweight='bold')

    for i in range(len(labels)):
        plt.plot(binned_lambdas, binned_alphas_2d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i])
        plt.plot(binned_lambdas, binned_stds_2d[i], c=colors[i], drawstyle='steps')
        plt.xlabel(r"Wavelength ($\AA$)")
        plt.ylabel(r"$\alpha_\lambda$")
        plt.xlim(x_min, x_max)
        plt.ylim(0, y_max)

    leg = plt.legend(frameon=False, loc='lower center')
    plt.setp(leg.texts, family='monospace')
    plt.text(0, 100, 'Tao Model')
    plt.tight_layout()

    plt.show()
    '''






