# Mostly copy / paste from generate_plots_ahab,
# and modified to make plots with ecliptic masked.

import matplotlib
from matplotlib import pyplot
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib import gridspec
import numpy as np
from astropy.io import fits
import sys  # for command line args
import pickle

'''
COMMAND LINE OPTIONS:

boss: 1 to use boss data, 0 to use sdss-II data
save: enter savekey. figures for the paper will get filenames based on the key. otherwise figures are not saved.
bootstrap: 1 to plot with bootstrap errors (and sometimes envelopes)
loadkey: alphas will be loaded based on this key
'''

def generate_binned_alphas(alphas, alpha_stds, wavelength_all, wavelength=None, boss=False, bin_offset=0, bin_width=50):
    # plot binned alpha vs wavelength (original)
    # "wavelength_all" is the one that determines the binned lambdas.
    # "wavelength" should match the alphas.

    if wavelength is None:
        wavelength = wavelength_all

    # find the bottom edge for the bins
    min_wav = wavelength_all[0]
    # add 50 to get inside the range, and another 50 bc these are right edges
    bin_start = min_wav - min_wav % bin_width + 2*bin_width + bin_offset
    binned_lambdas = np.arange(bin_start, wavelength_all[-1], bin_width)
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
            indices = np.where((wavelength > lmda - bin_width) & (wavelength < lmda) & np.logical_not(emission_line_mask))[
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

if __name__ == "__main__":

    #command line options
    if len(sys.argv) != 6:
        print("Usage: generate_plots.py [boss: 0, 1] [save: 0, savekey] [bootstrap: 0, 1] [show_plots: 0, 1] [loadkey]")
        exit(0)
    boss = int(sys.argv[1])
    save = sys.argv[2]
    bootstrap = int(sys.argv[3])
    show_plots = int(sys.argv[4])
    loadkey = sys.argv[5]

    # specify parameters
    alpha_direc = '../alphas_and_stds/'
    alpha_direc_boot = '../data/'  # '../alphas_and_stds'
    hdulist_direc = '../data/'  # '/Users/blakechellew/Documents/DustProject/BrandtFiles/'
    sdss_fluxfactor = 1.38
    boss_fluxfactor = 1.38
    if boss:
        fluxfactor = boss_fluxfactor
    else:
        fluxfactor = sdss_fluxfactor

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
    font_size = 12
    matplotlib.rcParams['axes.unicode_minus'] = False

    # load wavelengths
    wavelength_boss = np.load('../alphas_and_stds/wavelength_boss.npy')
    hdulist = fits.open(hdulist_direc + 'SDSS_allskyspec.fits')
    wavelength_sdss = np.array(hdulist[1].data)
    if boss:
        wavelength = wavelength_boss
    else:
        wavelength = wavelength_sdss

    # load alphas
    if boss:
        alphas = [np.load(alpha_direc + 'alphas_boss_iris_2d_north_' + loadkey + '.npy'),
                  np.load(alpha_direc + 'alphas_boss_iris_2d_south_' + loadkey + '.npy'),
                  np.load(alpha_direc + 'alphas_north_mask_ecliptic.npy'),
                  np.load(alpha_direc + 'alphas_south_mask_ecliptic.npy'),
                  np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_10.npy')]
        alpha_stds = [np.load(alpha_direc + 'alpha_stds_boss_iris_2d_north_' + loadkey + '.npy'),
                      np.load(alpha_direc + 'alpha_stds_boss_iris_2d_south_' + loadkey + '.npy'),
                      np.load(alpha_direc + 'alpha_stds_north_mask_ecliptic.npy'),
                      np.load(alpha_direc + 'alpha_stds_south_mask_ecliptic.npy'),
                      np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')]
    else:
        alphas = [np.load(alpha_direc + 'alphas_sdss_1d_' + loadkey + '.npy'),
                  np.load(alpha_direc + 'alphas_sdss_2d_' + loadkey + '.npy'),
                  np.load(alpha_direc + 'alphas_sdss_iris_2d_' + loadkey + '.npy'),
                  np.load(alpha_direc + 'alphas_sdss_iris_1d_' + loadkey + '.npy'),
                  np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_10.npy')] #don't change the last element
        alpha_stds = [np.load(alpha_direc + 'alpha_stds_sdss_1d_' + loadkey + '.npy'),
                      np.load(alpha_direc + 'alpha_stds_sdss_2d_' + loadkey + '.npy'),
                      np.load(alpha_direc + 'alpha_stds_sdss_iris_2d_' + loadkey + '.npy'),
                      np.load(alpha_direc + 'alpha_stds_sdss_iris_1d_' + loadkey + '.npy'),
                      np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')]
    # flux conversion factor:
    alphas = [a / fluxfactor for a in alphas]
    alpha_stds = [a / fluxfactor for a in alpha_stds]

    if bootstrap:
        if boss:
            with open(alpha_direc_boot + 'bootstrap_lower_boss' + loadkey + '.p', 'rb') as pf:
                bootstrap_lower = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_upper_boss' + loadkey + '.p', 'rb') as pf:
                bootstrap_upper = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_stds_boss' + loadkey + '.p', 'rb') as pf:
                bootstrap_stds = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_binned_lower_boss' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_lower = pickle.load(pf)
                # insert the ecliptic stuff
                bootstrap_binned_lower.insert(2, np.load('../data/bootstrap_binned_lower_north_mask_ecliptic.npy'))
                bootstrap_binned_lower.insert(3, np.load('../data/bootstrap_binned_lower_south_mask_ecliptic.npy'))
            with open(alpha_direc_boot + 'bootstrap_binned_upper_boss' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_upper = pickle.load(pf)
                # insert the ecliptic stuff
                bootstrap_binned_upper.insert(2, np.load('../data/bootstrap_binned_upper_north_mask_ecliptic.npy'))
                bootstrap_binned_upper.insert(3, np.load('../data/bootstrap_binned_upper_south_mask_ecliptic.npy'))
            with open(alpha_direc_boot + 'bootstrap_binned_stds_boss' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_stds = pickle.load(pf)
        else:
            with open(alpha_direc_boot + 'bootstrap_lower_sdss' + loadkey + '.p', 'rb') as pf:
                bootstrap_lower = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_upper_sdss' + loadkey + '.p', 'rb') as pf:
                bootstrap_upper = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_stds_sdss' + loadkey + '.p', 'rb') as pf:
                bootstrap_stds = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_binned_lower_sdss' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_lower = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_binned_upper_sdss' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_upper = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_binned_stds_sdss' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_stds = pickle.load(pf)
            # load the relevant boss alphas
            bootstrap_binned_lower_boss = np.load(alpha_direc_boot + 'bootstrap_binned_lower_boss_to_sdss' + loadkey + '.npy') / boss_fluxfactor
            bootstrap_binned_upper_boss = np.load(alpha_direc_boot + 'bootstrap_binned_upper_boss_to_sdss' + loadkey + '.npy') / boss_fluxfactor
            bootstrap_binned_stds_boss = np.load(alpha_direc_boot + 'bootstrap_binned_stds_boss_to_sdss' + loadkey + '.npy') / boss_fluxfactor
        # flux conversion factor
        bootstrap_lower = [b / fluxfactor for b in bootstrap_lower]
        bootstrap_upper = [b / fluxfactor for b in bootstrap_upper]
        bootstrap_stds = [b / fluxfactor for b in bootstrap_stds]
        bootstrap_binned_lower = [b / fluxfactor for b in bootstrap_binned_lower]
        bootstrap_binned_upper = [b / fluxfactor for b in bootstrap_binned_upper]
        bootstrap_binned_stds = [b / fluxfactor for b in bootstrap_binned_stds]

if __name__ == "__main__":
    # bin the regular spectra
    binned_lambdas, binned_alphas, binned_stds = generate_binned_alphas(alphas, alpha_stds, wavelength, boss=boss)
    if not boss: #calculate binned spectrum for boss, but using bins based on sdss
        binned_lambdas_boss, binned_alphas_boss, binned_stds_boss = generate_binned_alphas([alphas[-1]], [alpha_stds[-1]], wavelength, wavelength_boss, boss=boss)

# plot binned alphas
# takes alphas already binned. use generate_binned_alphas
def plot_binned(alpha_indices, colors, labels, envelope=False, thicks=None):
    fig, ax = plt.subplots(figsize=(6, 4.8))

    if thicks is None:
        thicks = [1.5 for l in labels]

    for i, idx in enumerate(alpha_indices):
        ax.plot(binned_lambdas, binned_alphas[idx], c=colors[i], drawstyle='steps', label=labels[i], linewidth=thicks[i])
        # if bootstrap:
        #     ax.plot(binned_lambdas, bootstrap_binned_stds[idx], c=colors[i], drawstyle='steps', linestyle='--', linewidth=thicks[i])
        # else:
        #     ax.plot(binned_lambdas, binned_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        if envelope:
            ax.fill_between(binned_lambdas, bootstrap_binned_lower[idx], bootstrap_binned_upper[idx], linewidth=0.0,
                            color=colors[i], alpha=0.2, step='pre')

    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(200))
    ax.yaxis.set_minor_locator(MultipleLocator(0.02))
    # ax.figure(figsize=(6, 5))
    ax.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    ax.set_ylabel(r"$\alpha_\lambda^{\prime}$ = $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
    if boss:
        ax.set_xlim(3700, 10000)
    else:
        ax.set_xlim(3850, 9200)
    ax.set_ylim(0, .4)
    ax.legend(frameon=False)
    return ax

if __name__ == "__main__":

    envelope = False

    ##############################
    # PLOT NORTH AND SOUTH BOSS
    ##############################
    if boss:
        ax = plot_binned([1, 0], ['#BB5566', '#004488'], ['South', 'North'], envelope=envelope, thicks=[1.5, 1.7])
        ax.text(0.5, 0.2, "Bootstrap Uncertainties", transform=ax.transAxes, fontsize=font_size, horizontalalignment='center')
        ax.arrow(0.5, 0.17, 0, -0.04, transform=ax.transAxes, width=.001, color='k', head_width=.01)
        if save != '0' and bootstrap:
            plt.savefig('../paper_figures/boss_north_south_' + save + '.pdf', bbox_inches='tight')
            plt.clf()
        elif show_plots:
            plt.show()

    ##############################
    # NORTH AND SOUTH BOSS WITH ECLIPTIC MASKED
    ##############################
    if boss:
        ax = plot_binned([3, 2], ['#BB5566', '#004488'], ['South (ecliptic masked)', 'North (ecliptic masked)'], envelope=envelope, thicks=[1.5, 1.7])
        ax.text(0.5, 0.2, "Bootstrap Uncertainties", transform=ax.transAxes, fontsize=font_size,
                horizontalalignment='center')
        ax.arrow(0.5, 0.17, 0, -0.04, transform=ax.transAxes, width=.001, color='k', head_width=.01)
        if save != '0' and bootstrap:
            plt.savefig('../paper_figures/boss_north_south_' + save + '.pdf', bbox_inches='tight')
            plt.clf()
        elif show_plots:
            plt.show()

    ##############################
    # SOUTH WITH / WITHOUT ECLIPTIC
    ##############################
    if boss:
        ax = plot_binned([1, 3], ['#BB5566', '#004488'], ['South', 'South (ecliptic masked)'], envelope=envelope, thicks=[1.5, 1.7])
        ax.text(0.5, 0.2, "Bootstrap Uncertainties", transform=ax.transAxes, fontsize=font_size,
                horizontalalignment='center')
        ax.arrow(0.5, 0.17, 0, -0.04, transform=ax.transAxes, width=.001, color='k', head_width=.01)
        if save != '0' and bootstrap:
            plt.savefig('../paper_figures/boss_north_south_' + save + '.pdf', bbox_inches='tight')
            plt.clf()
        elif show_plots:
            plt.show()

    ##############################
    # NORTH WITH / WITHOUT ECLIPTIC
    ##############################
    if boss:
        ax = plot_binned([0, 2], ['#BB5566', '#004488'], ['North', 'North (ecliptic masked)'], envelope=envelope, thicks=[1.5, 1.7])
        ax.text(0.5, 0.2, "Bootstrap Uncertainties", transform=ax.transAxes, fontsize=font_size,
                horizontalalignment='center')
        ax.arrow(0.5, 0.17, 0, -0.04, transform=ax.transAxes, width=.001, color='k', head_width=.01)
        if save != '0' and bootstrap:
            plt.savefig('../paper_figures/boss_north_south_' + save + '.pdf', bbox_inches='tight')
            plt.clf()
        elif show_plots:
            plt.show()