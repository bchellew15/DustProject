# flux conversion factor is handled here, not in reproduce_figs

# should allow more things to be passed from cmd line

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys  # for command line args
from math import floor  # for binning range

'''
COMMAND LINE OPTIONS:

boss: 1 to use boss data, 0 to use sdss-II data
save: enter savekey. figures for the paper will get filenames based on the key. otherwise figures are not saved.
bootstrap: 1 to plot with bootstrap errors (and sometimes envelopes)
loadkey: alphas will be loaded based on this key
'''

# command line options
if len(sys.argv) != 6:
    print("Usage: generate_plots.py [boss: 0, 1] [save: 0, savekey] [bootstrap: 0, 1] [loadkey] [binned: 0, 1] ")
    exit(0)
boss = int(sys.argv[1])
save = sys.argv[2]
bootstrap = int(sys.argv[3])
loadkey = sys.argv[4]
binned = int(sys.argv[5])

alpha_direc = '../alphas_and_stds/'
alpha_direc_boot = '../data/'  # '../alphas_and_stds'
alpha_direc_boot = '../alphas_and_stds'
hdulist_direc = '../data/'  # '/Users/blakechellew/Documents/DustProject/BrandtFiles/'
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

# boss alphas:
alphas = [
    np.load('../alphas_and_stds/partial_jacknife_alphas_boss_iris_1d_012720_30.npy')[0],
    np.load('../alphas_and_stds/partial_jacknife_alphas_boss_iris_1d_012720_30.npy')[1],
    np.load('../alphas_and_stds/partial_jacknife_alphas_boss_iris_1d_012720_30.npy')[2],
    np.load('../alphas_and_stds/partial_jacknife_alphas_boss_iris_1d_012720_30.npy')[3],
    np.load('../alphas_and_stds/partial_jacknife_alphas_boss_iris_1d_012720_30.npy')[4]
    # np.load(alpha_direc + 'alphas_092220.npy'),
    # np.load(alpha_direc + 'alphas_above_50_092220.npy'),
    # np.load(alpha_direc + 'alphas_above_50_A_092220.npy')
    # np.load(alpha_direc + 'alphas_above_50_B_092220.npy'),
    # np.load(alpha_direc + 'alphas_above_50_C_092220.npy'),
]
alpha_stds = [
    np.load(alpha_direc + 'alpha_stds_boss_91119_10.npy'),
    np.load(alpha_direc + 'alpha_stds_boss_91119_10.npy'),
    np.load(alpha_direc + 'alpha_stds_boss_91119_10.npy'),
    np.load(alpha_direc + 'alpha_stds_boss_91119_10.npy'),
    np.load(alpha_direc + 'alpha_stds_boss_91119_10.npy')
    # np.load(alpha_direc + 'alpha_stds_092220.npy'),
    # np.load(alpha_direc + 'alpha_stds_above_50_092220.npy'),
    # np.load(alpha_direc + 'alpha_stds_above_50_A_092220.npy')
    # np.load(alpha_direc + 'alpha_stds_above_50_B_092220.npy'),
    # np.load(alpha_direc + 'alpha_stds_above_50_C_092220.npy'),
]

# flux conversion factor:
if boss:
    alphas = [a / boss_fluxfactor for a in alphas]
    alpha_stds = [a / boss_fluxfactor for a in alpha_stds]
else:
    alphas = [a / sdss_fluxfactor for a in alphas]
    alpha_stds = [a / sdss_fluxfactor for a in alpha_stds]

num_arrays = len(alphas)


# plot unbinned spectra (wavelength ranges: 4830-5040 and 6530-6770)
def plot_emissions(alpha_indices, labels, colors):
    plt.figure(figsize=(12, 5))

    # plot 4830 - 5040
    plt.subplot(1, 2, 1)
    for i, idx in enumerate(alpha_indices):
        plt.plot(wavelength, alphas[idx], c=colors[i], drawstyle='steps', label=labels[i])
        if bootstrap:
            # plt.fill_between(wavelength, bootstrap_lower[idx], bootstrap_upper[idx], linewidth=0.0, color=colors[i], alpha=0.5, step='pre')
            plt.plot(wavelength, bootstrap_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        else:
            plt.plot(wavelength, alpha_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')

    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"$\alpha_\lambda$")
    plt.legend(loc='upper center', frameon=False)
    plt.xlim(4830, 5040)
    plt.ylim(0, 1)
    xcoords = [4863, 4960, 5008]
    for xc in xcoords:
        plt.axvline(x=xc, color='k', linewidth=1, linestyle='--')
    plt.text(4853, 0.4, r"H$\beta$")
    plt.text(4943, 0.4, "O[III]")
    plt.text(4991, 0.4, "O[III]")

    # line from 03 continuum::
    # plt.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
    # actual continuum for NII:
    # plt.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')

    # plot 6530 - 6770 (original vs tao)
    plt.subplot(1, 2, 2)
    for i, idx in enumerate(alpha_indices):
        plt.plot(wavelength, alphas[idx], c=colors[i], drawstyle='steps', label=labels[i])
        if bootstrap:
            # plt.fill_between(wavelength, bootstrap_lower[idx], bootstrap_upper[idx], linewidth=0.0, color=colors[i], alpha=0.5, step='pre')
            plt.plot(wavelength, bootstrap_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        else:
            plt.plot(wavelength, alpha_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')

    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"$\alpha_\lambda$")
    plt.legend(loc='upper center', frameon=False)
    plt.xlim(6530, 6770)
    plt.ylim(0, 1)
    xcoords = [6550, 6565, 6585, 6718, 6733]
    for xc in xcoords:
        plt.axvline(x=xc, color='k', linewidth=1, linestyle='--')
    plt.text(6535, 0.4, "N[II]")
    plt.text(6568, 0.9, r"H$\alpha$")
    plt.text(6590, 0.6, "N[II]")
    plt.text(6700, 0.6, "S[II]")
    plt.text(6738, 0.5, "S[II]")

    # line from 03 continuum::
    # plt.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
    # actual continuum for NII:
    # plt.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')

    plt.tight_layout()

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
    # plot binned alpha vs wavelength (original)
    # wavelength_all is the one that determines the binned lambdas

    if wavelength is None:
        wavelength = wavelength_all

    lambda_range = wavelength_all[-1] - wavelength_all[0]
    left_over = lambda_range - 50 * floor(lambda_range / 50)
    binned_lambdas = np.arange(wavelength_all[1] + left_over / 2, wavelength_all[-1],
                               50)  # [1] to avoid going over left edge
    binned_alphas = []
    binned_stds = []

    # mask emission lines
    emission_line_mask = np.zeros(len(wavelength), dtype=int)
    emission_lines = [4863, 4960, 5008, 5877, 6550, 6565, 6585, 6718, 6733]
    if boss:
        emission_lines.insert(0, 3727)
    for line in emission_lines:
        peak_idx = np.argmin(np.abs(wavelength - line))
        emission_line_mask[peak_idx - 3:peak_idx + 4] = 1

    for i in range(len(alphas)):

        binned_alpha_arr = np.zeros(binned_lambdas.shape)
        binned_std_arr = np.zeros(binned_lambdas.shape)
        if alphas[i].ndim > 1:
            binned_alpha_arr = np.zeros((alphas[i].shape[0], binned_lambdas.shape[0]))
            binned_std_arr = np.zeros((alphas[i].shape[0], binned_lambdas.shape[0]))
        for j, lmda in enumerate(binned_lambdas):
            indices = np.where((wavelength > lmda - 50) & (wavelength < lmda) & np.logical_not(emission_line_mask))[
                0]  # test
            if alphas[i].ndim > 1:
                relevant_alphas = alphas[i][:, indices]
                relevant_stds = alpha_stds[i][:, indices]
            else:
                relevant_alphas = alphas[i][indices]
                relevant_stds = alpha_stds[i][indices]

            # weighted average:
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
                binned_alpha_arr[:, j] = avg1
                binned_std_arr[:, j] = np.sqrt(avg2)
            else:
                binned_alpha_arr[j] = avg1
                binned_std_arr[j] = np.sqrt(avg2)

        binned_alphas.append(binned_alpha_arr)
        binned_stds.append(binned_std_arr)

    return binned_lambdas, binned_alphas, binned_stds

#bin alphas
if binned:
    binned_lambdas, binned_alphas, binned_stds = generate_binned_alphas(alphas, alpha_stds, wavelength)

if boss:
    y_max = 0.3
    x_min = 3700
    x_max = 10000
else:
    y_max = 0.3
    x_min = 3850
    x_max = 9200

# plot binned alphas
# takes alphas already binned. use generate_binned_alphas
def plot_binned(alpha_indices, colors, labels, envelope=False):
    for i, idx in enumerate(alpha_indices):
        plt.plot(binned_lambdas, binned_alphas[idx], c=colors[i], drawstyle='steps', label=labels[i])
        if bootstrap:
            plt.plot(binned_lambdas, bootstrap_binned_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        else:
            plt.plot(binned_lambdas, binned_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        if envelope:
            plt.fill_between(binned_lambdas, bootstrap_binned_lower[idx], bootstrap_binned_upper[idx], linewidth=0.0,
                             color=colors[i], alpha=0.2, step='pre')

    y_max = 0.35 #why manually adjusting this?
    plt.xlabel(r"Wavelength ($\AA$)")
    plt.ylabel(r"$\alpha_\lambda$")
    plt.xlim(x_min, x_max)
    plt.ylim(0, y_max)
    plt.legend(frameon=False)
    plt.show()
    y_max = 0.3


indices = [0, 1, 2, 3, 4]
labels = ['2389', '2390', '2384', '2395', '2386']
colors = ['k', 'r', 'g', 'b', 'pink']

# plot north and south on same plot (boss)
if binned:
    plot_binned(indices, colors, labels)
    if save != '0':
        plt.savefig('../paper_figures/' + save + '.pdf')
    else:
        plt.show()
else:
    plot_emissions(indices, labels, colors)
    if save != '0':
        plt.savefig('../paper_figures/' + save + '.pdf')
    else:
        plt.show()
