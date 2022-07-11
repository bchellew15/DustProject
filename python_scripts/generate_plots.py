# this should run on ahab.
# run this after generate_alphas_ahab.py. If want to use bootstrapping, also have to run bootstrap_ahab.py first.
# also need to run create_bootstrap_envelopes to get the bootstrap envelopes.

# generate plots of correlation spectra (overall and certain sections)
# need to run 2x to cover all figures for the paper:
#   boss=0 and bootstrap=1
#   boss=1 and bootstrap=1
# flux conversion factor is handled here, not in reproduce_figs

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
        alphas = [np.load(alpha_direc + 'alphas_boss_iris_2d_north_' + loadkey + '.npy'), \
                  np.load(alpha_direc + 'alphas_boss_iris_2d_south_' + loadkey + '.npy'), \
                  np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_10.npy')]
        alpha_stds = [np.load(alpha_direc + 'alpha_stds_boss_iris_2d_north_' + loadkey + '.npy'), \
                      np.load(alpha_direc + 'alpha_stds_boss_iris_2d_south_' + loadkey + '.npy'), \
                      np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')]
        correction_factors = [np.load('../alphas_and_stds/correction_factor_boss_iris_north_' + loadkey + '_smooth.npy'),
                              np.load('../alphas_and_stds/correction_factor_boss_iris_south_' + loadkey + '_smooth.npy'),
                              np.load('../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_10_smooth.npy')]
    else:
        alphas = [np.load(alpha_direc + 'alphas_sdss_1d_' + loadkey + '.npy'), \
                  np.load(alpha_direc + 'alphas_sdss_2d_' + loadkey + '.npy'), \
                  np.load(alpha_direc + 'alphas_sdss_iris_2d_' + loadkey + '.npy'), \
                  np.load(alpha_direc + 'alphas_sdss_iris_1d_' + loadkey + '.npy'), \
                  np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_10.npy')] #don't change the last element
        alpha_stds = [np.load(alpha_direc + 'alpha_stds_sdss_1d_' + loadkey + '.npy'), \
                      np.load(alpha_direc + 'alpha_stds_sdss_2d_' + loadkey + '.npy'), \
                      np.load(alpha_direc + 'alpha_stds_sdss_iris_2d_' + loadkey + '.npy'), \
                      np.load(alpha_direc + 'alpha_stds_sdss_iris_1d_' + loadkey + '.npy'), \
                      np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')]
        correction_factors = [np.ones(len(wavelength_sdss)),
                              np.load('../alphas_and_stds/correction_factor_sdss_sfd_smooth.npy'),
                              np.load('../alphas_and_stds/correction_factor_sdss_iris_smooth.npy'),
                              np.ones(len(wavelength_sdss)),
                              np.load('../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_10_smooth.npy')]
    # flux conversion factor and correction factor:
    alphas = [a / fluxfactor * corr for a, corr in zip(alphas, correction_factors)]
    alpha_stds = [a / fluxfactor * corr for a, corr in zip(alpha_stds, correction_factors)]

    # bin the correction factors
    _, binned_corrections, _ = generate_binned_alphas(correction_factors, 4*[np.ones(len(correction_factors[0]))] + [np.ones(len(correction_factors[-1]))],
                                                      wavelength, boss=boss)

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
            with open(alpha_direc_boot + 'bootstrap_binned_upper_boss' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_upper = pickle.load(pf)
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
            bootstrap_binned_lower_boss = np.load(alpha_direc_boot + 'bootstrap_binned_lower_boss_to_sdss' + loadkey + '.npy') / boss_fluxfactor * binned_corrections[-1]
            bootstrap_binned_upper_boss = np.load(alpha_direc_boot + 'bootstrap_binned_upper_boss_to_sdss' + loadkey + '.npy') / boss_fluxfactor * binned_corrections[-1]
            bootstrap_binned_stds_boss = np.load(alpha_direc_boot + 'bootstrap_binned_stds_boss_to_sdss' + loadkey + '.npy') / boss_fluxfactor * binned_corrections[-1]
        # flux conversion factor
        bootstrap_lower = [b / fluxfactor * corr for b, corr in zip(bootstrap_lower, correction_factors)]
        bootstrap_upper = [b / fluxfactor * corr for b, corr in zip(bootstrap_upper, correction_factors)]
        bootstrap_stds = [b / fluxfactor * corr for b, corr in zip(bootstrap_stds, correction_factors)]
        bootstrap_binned_lower = [b / fluxfactor * corr for b, corr in zip(bootstrap_binned_lower, binned_corrections)]
        bootstrap_binned_upper = [b / fluxfactor * corr for b, corr in zip(bootstrap_binned_upper, binned_corrections)]
        bootstrap_binned_stds = [b / fluxfactor * corr for b, corr in zip(bootstrap_binned_stds, binned_corrections)]

# plot unbinned spectra (wavelength ranges: 4830-5040 and 6530-6770)
def plot_emissions(alpha_indices, labels, colors, show_o3=False, thicks=None, zorders=None):
    fig = plt.figure(figsize=(12, 4.8)) # TEST (changed all from 5 to 4)

    if thicks is None:
        thicks = [1.5 for _ in labels]
    if zorders is None:
        zorders = [1 for _ in labels]

    num_plots = 2
    if show_o3:
        num_plots = 3

    spec = gridspec.GridSpec(ncols=num_plots, nrows=1, width_ratios=[1, 3, 3], wspace=0.05, hspace=0, left=0.1, right=0.9)

    # plot 4830 - 5040
    # ax1 = plt.subplot(1, num_plots, num_plots-1)
    ax1 = fig.add_subplot(spec[-2])
    for i, idx in enumerate(alpha_indices):
        ax1.plot(wavelength, alphas[idx], c=colors[i], drawstyle='steps-mid', label=labels[i], linewidth=thicks[i], zorder=zorders[i])
        if bootstrap:
            # ax1.fill_between(wavelength, bootstrap_lower[idx], bootstrap_upper[idx], linewidth=0.0, color=colors[i], alpha=0.5, step='pre')
            ax1.plot(wavelength, bootstrap_stds[idx], c=colors[i], drawstyle='steps-mid', linestyle='--', linewidth=thicks[i], zorder=zorders[i])
        else:
            ax1.plot(wavelength, alpha_stds[idx], c=colors[i], drawstyle='steps-mid', linestyle='--')

    if show_o3:
        # ax1.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
        ax1.set_yticklabels([])
        fig.supxlabel(r"Wavelength ($\mathrm{\AA}$)")
    else:
        ax1.set_ylabel(r"$\alpha_\lambda^{\prime}$ $\propto$ $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
    ax1.legend(loc='upper center', bbox_to_anchor=(0.45, 0.95), frameon=False)
    ax1.set_xlim(4805, 5045)
    ax1.set_ylim(-0.16, 0.96)
    xcoords = [4863, 4960, 5008]
    for xc in xcoords:
        ax1.axvline(x=xc, color='k', linewidth=1, linestyle='--')
    ax1.hlines(0, 4000, 6000, color='gray')
    ax1.text(4856, 0.4, r"H$\beta$", fontsize=font_size, backgroundcolor='white')
    ax1.text(4946, 0.4, "[O III]", fontsize=font_size, backgroundcolor='white')
    ax1.text(4994, 0.4, "[O III]", fontsize=font_size, backgroundcolor='white')

    ax1.xaxis.set_major_locator(MultipleLocator(50))
    ax1.xaxis.set_minor_locator(MultipleLocator(10))
    ax1.yaxis.set_major_locator(MultipleLocator(0.2))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.04))

    #line from 03 continuum::
    #ax1.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
    #actual continuum for NII:
    #ax1.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')

    #plot 6530 - 6770 (original vs tao)
    # ax2 = plt.subplot(1, num_plots, num_plots)
    ax2 = fig.add_subplot(spec[-1])
    for i, idx in enumerate(alpha_indices):
        ax2.plot(wavelength, alphas[idx], c=colors[i], drawstyle='steps-mid', label=labels[i], linewidth=thicks[i], zorder=zorders[i])
        if bootstrap:
            #ax2.fill_between(wavelength, bootstrap_lower[idx], bootstrap_upper[idx], linewidth=0.0, color=colors[i], alpha=0.5, step='pre')
            ax2.plot(wavelength, bootstrap_stds[idx], c=colors[i], drawstyle='steps-mid', linestyle='--', linewidth=thicks[i], zorder=zorders[i])
        else:
            ax2.plot(wavelength, alpha_stds[idx], c=colors[i], drawstyle='steps-mid', linestyle='--')

    if show_o3:
        ax2.set_yticklabels([])
    else:
        ax2.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
        ax2.set_ylabel(r"$\alpha_\lambda^{\prime}$ $\propto$ $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
    # ax2.legend(loc='upper center', frameon=False)
    ax2.set_xlim(6530, 6770)
    ax2.set_ylim(-0.16, 0.96)
    xcoords = [6550, 6565, 6585, 6718, 6733]
    for xc in xcoords:
        ax2.axvline(x=xc, color='k', linewidth=1, linestyle='--')
    ax2.hlines(0, 6000, 7000, color='gray')
    ax2.text(6535, 0.4, "[N II]", fontsize=font_size, backgroundcolor='white', bbox=dict(facecolor='none', edgecolor='none', pad=0, fc='white'))
    ax2.text(6567, 0.89, r"H$\alpha$", fontsize=font_size)
    ax2.text(6575, 0.72, "[N II]", fontsize=font_size, backgroundcolor='white')
    ax2.text(6705, 0.68, "[S II]", fontsize=font_size, backgroundcolor='white')
    ax2.text(6728, 0.55, "[S II]", fontsize=font_size, backgroundcolor='white')

    ax2.xaxis.set_major_locator(MultipleLocator(50))
    ax2.xaxis.set_minor_locator(MultipleLocator(10))
    ax2.yaxis.set_major_locator(MultipleLocator(0.2))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.04))

    if show_o3:
        # ax3 = plt.subplot(1, num_plots, 1)
        ax3 = fig.add_subplot(spec[0])
        for i, idx in enumerate(alpha_indices):
            ax3.plot(wavelength, alphas[idx], c=colors[i], drawstyle='steps-mid', label=labels[i], linewidth=thicks[i], zorder=zorders[i])
            if bootstrap:
                # ax2.fill_between(wavelength, bootstrap_lower[idx], bootstrap_upper[idx], linewidth=0.0, color=colors[i], alpha=0.5, step='pre')
                ax3.plot(wavelength, bootstrap_stds[idx], c=colors[i], drawstyle='steps-mid', linestyle='--', linewidth=thicks[i], zorder=zorders[i])
            else:
                ax3.plot(wavelength, alpha_stds[idx], c=colors[i], drawstyle='steps-mid', linestyle='--')

        ax3.set_ylabel(r"$\alpha_\lambda^{\prime}$ $\propto$ $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
        # ax3.legend(loc='upper center', frameon=False)
        ax3.set_xlim(3712, 3742)
        ax3.set_ylim(-0.16, 0.96)
        # ax3.set_ylim(-0.24, 0.99-0.24)
        # ax3.set_ylim(-0.15, 0.3)
        xcoords = [3726, 3729]
        for xc in xcoords:
            ax3.axvline(x=xc, color='k', linewidth=1, linestyle='--')
        ax3.hlines(0, 3700, 3800, color='gray')
        ax3.text(3718, 0.54, "[O II]", fontsize=font_size, backgroundcolor='white')
        ax3.text(3728, 0.64, "[O II]", fontsize=font_size, backgroundcolor='white')

        ax3.xaxis.set_major_locator(MultipleLocator(20))
        # ax3.xaxis.set_minor_locator(MultipleLocator(10))
        ax3.yaxis.set_major_locator(MultipleLocator(0.2))
        ax3.yaxis.set_minor_locator(MultipleLocator(0.04))


    #line from 03 continuum::
    #ax2.axhline(y=0.14898818311840933, color='r', linewidth=1, linestyle='--')
    #actual continuum for NII:
    #ax2.axhline(y=0.17930096676470586, color='r', linewidth=1, linestyle='--')


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
        if bootstrap:
            ax.plot(binned_lambdas, bootstrap_binned_stds[idx], c=colors[i], drawstyle='steps', linestyle='--', linewidth=thicks[i])
        else:
            ax.plot(binned_lambdas, binned_stds[idx], c=colors[i], drawstyle='steps', linestyle='--')
        if envelope:
            ax.fill_between(binned_lambdas, bootstrap_binned_lower[idx], bootstrap_binned_upper[idx], linewidth=0.0,
                            color=colors[i], alpha=0.2, step='pre')

    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(200))
    ax.yaxis.set_minor_locator(MultipleLocator(0.02))
    # ax.figure(figsize=(6, 5))
    ax.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    ax.set_ylabel(r"$\alpha_\lambda^{\prime}$ $\propto$ $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
    if boss:
        ax.set_xlim(3700, 10000)
    else:
        ax.set_xlim(3850, 9200)
    ax.set_ylim(0, .29)
    ax.legend(frameon=False)
    return ax

if __name__ == "__main__":

    ###############################
    # PLOT: ALPHAS VS BC03 WITH RAD. TRANSFER
    ###############################


    #################################
    # PLOT: BOSS NORTH VS SOUTH EMISSIONS
    #################################
    if boss:
        plot_emissions([1, 0], ["South", "North"], ['#BB5566', '#004488'], show_o3=True, thicks=[1.5, 1.7], zorders=[3, 2])
        # removed "Full Sky" at index 2 with color '#DDAA33'
        if save != '0' and bootstrap:
            print("saving as:", '../paper_figures/unbinned_' + save + '.pdf')
            plt.savefig('../paper_figures/unbinned_' + save + '.pdf', bbox_inches='tight')
            plt.clf()
        elif show_plots:
            plt.show()

    #################################
    # 2-PANEL BOSS VS SDSS (both 1d, etc.)
    #################################
    if not boss and bootstrap:
        x_min = 3850 # 3700
        x_max = 9200  # 10000

        plt.figure(figsize=(12, 4.8))

        ax1 = plt.subplot(1, 2, 1)

        a1, = ax1.plot(binned_lambdas, binned_alphas[2], c='k', drawstyle='steps', label='SDSS-II')
        ax1.fill_between(binned_lambdas, bootstrap_binned_lower[2], bootstrap_binned_upper[2], linewidth=0.0, color='k', alpha=0.2, step='pre')
        a2, = ax1.plot(binned_lambdas, bootstrap_binned_stds[2], c='k', drawstyle='steps', linestyle='--', label='Formal Uncertainties')
        a3, = ax1.plot(binned_lambdas, binned_stds[2], c='m', drawstyle='steps', linestyle='dotted', label='Bootstrap Uncertainties')

        legend1 = ax1.legend([a1], ['SDSS-II'], frameon=False, loc='upper left')
        legend2 = ax1.legend([a2, a3], ['Bootstrap Uncertainties', 'Formal Uncertainties'], frameon=False, loc='lower center', bbox_to_anchor=(0.5, 0.1))
        pyplot.gca().add_artist(legend1)
        ax1.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
        ax1.set_ylabel(r"$\alpha_\lambda^{\prime}$ $\propto$ $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
        ax1.set_xlim(x_min, x_max)
        ax1_ymin, ax1_ymax = 0, 0.29
        ax1.set_ylim(ax1_ymin, ax1_ymax)

        ax1.xaxis.set_major_locator(MultipleLocator(1000))
        ax1.xaxis.set_minor_locator(MultipleLocator(200))
        ax1.yaxis.set_major_locator(MultipleLocator(0.1))
        ax1.yaxis.set_minor_locator(MultipleLocator(0.02))

        ax2 = plt.subplot(1, 2, 2)

        a4, = ax2.plot(binned_lambdas_boss, binned_alphas_boss[0], c='k', drawstyle='steps', label='BOSS')
        ax2.fill_between(binned_lambdas_boss, bootstrap_binned_lower_boss, bootstrap_binned_upper_boss, linewidth=0.0, color='k', alpha=0.2, step='pre')
        a5, = ax2.plot(binned_lambdas_boss, bootstrap_binned_stds_boss, c='k', drawstyle='steps', linestyle='--', label='Formal Uncertainties')
        a6, = ax2.plot(binned_lambdas_boss, binned_stds_boss[0], c='m', drawstyle='steps', linestyle='dotted', label='Bootstrap Uncertainties')

        legend3 = ax2.legend([a4], ['BOSS'], frameon=False, loc='upper left')
        legend4 = ax2.legend([a5, a6], ['Bootstrap Uncertainties', 'Formal Uncertainties'], frameon=False, loc='lower center', bbox_to_anchor=(0.5, 0.1))
        pyplot.gca().add_artist(legend3)
        ax2.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
        ax2.set_ylabel(r"$\alpha_\lambda^{\prime}$ $\propto$ $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
        ax2.set_xlim(x_min, x_max)
        ax2_ymin, ax2_ymax = 0, 0.29
        ax2.set_ylim(ax2_ymin, ax2_ymax)

        ax2.xaxis.set_major_locator(MultipleLocator(1000))
        ax2.xaxis.set_minor_locator(MultipleLocator(200))
        ax2.yaxis.set_major_locator(MultipleLocator(0.1))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.02))

        if save != '0':
            plt.savefig('../paper_figures/compare_boss_sdss_' + save + '.pdf', bbox_inches='tight')
            plt.clf()
        elif show_plots:
            plt.show()

    ##############################
    # PLOT ORIGINAL AND 2 MODS
    # (new model, IRIS)
    ##############################
    if not boss:
        plot_binned([0, 1, 3], ['#004488', '#BB5566', '#DDAA33'], ['SFD', 'Nonlinear Model', 'IRIS'])
        if save != '0' and bootstrap:
            plt.savefig('../paper_figures/all_3_mods_' + save + '.pdf', bbox_inches='tight')
            plt.clf()
        elif show_plots:
            plt.show()

    ##############################
    # PLOT NORTH AND SOUTH BOSS
    ##############################
    if boss:
        envelope = bootstrap
        ax = plot_binned([1, 0], ['#BB5566', '#004488'], ['South', 'North'], envelope=envelope, thicks=[1.5, 1.7])
        ax.text(0.5, 0.2, "Bootstrap Uncertainties", transform=ax.transAxes, fontsize=font_size, horizontalalignment='center')
        ax.arrow(0.5, 0.17, 0, -0.04, transform=ax.transAxes, width=.001, color='k', head_width=.01)
        if save != '0' and bootstrap:
            plt.savefig('../paper_figures/boss_north_south_' + save + '.pdf', bbox_inches='tight')
            plt.clf()
        elif show_plots:
            plt.show()

    ##############################
    # THRESHOLD PLOTS
    ##############################
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
        correction_factors_thresh = [np.load('../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_10_smooth.npy'),
                                     np.load('../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_15_smooth.npy'),
                                     np.load('../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_20_smooth.npy'),
                                     np.load('../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_25_smooth.npy'),
                                     np.load('../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_30_smooth.npy')]
        # test
        plt.plot(wavelength, correction_factors_thresh[0], label='10')
        plt.plot(wavelength, correction_factors_thresh[1], label='15')
        plt.plot(wavelength, correction_factors_thresh[2], label='20')
        plt.plot(wavelength, correction_factors_thresh[3], label='25')
        plt.plot(wavelength, correction_factors_thresh[4], label='30')
        plt.legend()
        plt.show()

        # bin the correction factors
        _, binned_corrections_thresh, _ = generate_binned_alphas(correction_factors_thresh, 5*[np.ones(len(correction_factors_thresh[0]))], wavelength, boss=boss)

        if bootstrap:
            with open(alpha_direc_boot + 'bootstrap_binned_lower_thresh_1d' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_lower_thresh_1d = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_binned_upper_thresh_1d' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_upper_thresh_1d = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_binned_stds_thresh_1d' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_stds_thresh_1d = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_binned_lower_thresh_2d' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_lower_thresh_2d = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_binned_upper_thresh_2d' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_upper_thresh_2d = pickle.load(pf)
            with open(alpha_direc_boot + 'bootstrap_binned_stds_thresh_2d' + loadkey + '.p', 'rb') as pf:
                bootstrap_binned_stds_thresh_2d = pickle.load(pf)

        # flux factor corrections:
        alphas_thresh_1d = [a / boss_fluxfactor for a in alphas_thresh_1d]
        alphas_stds_thresh_1d = [a / boss_fluxfactor for a in alpha_stds_thresh_1d]
        alphas_thresh_2d = [a / boss_fluxfactor * corr for a, corr in zip(alphas_thresh_2d, correction_factors_thresh)]
        alphas_stds_thresh_2d = [a / boss_fluxfactor * corr for a, corr in zip(alpha_stds_thresh_2d, correction_factors_thresh)]
        if bootstrap:
            bootstrap_binned_lower_thresh_1d = [b / boss_fluxfactor for b in bootstrap_binned_lower_thresh_1d]
            bootstrap_binned_upper_thresh_1d = [b / boss_fluxfactor for b in bootstrap_binned_upper_thresh_1d]
            bootstrap_binned_stds_thresh_1d = [b / boss_fluxfactor for b in bootstrap_binned_stds_thresh_1d]
            bootstrap_binned_lower_thresh_2d = [b / boss_fluxfactor * corr for b, corr in zip(bootstrap_binned_lower_thresh_2d, binned_corrections_thresh)]
            bootstrap_binned_upper_thresh_2d = [b / boss_fluxfactor * corr for b, corr in zip(bootstrap_binned_upper_thresh_2d, binned_corrections_thresh)]
            bootstrap_binned_stds_thresh_2d = [b / boss_fluxfactor * corr for b, corr in zip(bootstrap_binned_stds_thresh_2d, binned_corrections_thresh)]

        #binning
        binned_lambdas, binned_alphas_1d, binned_stds_1d = generate_binned_alphas(alphas_thresh_1d, alpha_stds_thresh_1d, wavelength, boss=boss)
        binned_lambdas, binned_alphas_2d, binned_stds_2d = generate_binned_alphas(alphas_thresh_2d, alpha_stds_thresh_2d, wavelength, boss=boss)

        fig = plt.figure(figsize=(12, 4.8))
        x_min = 3700
        x_max = 10000
        y_max = 0.29

        ax1 = fig.add_subplot(121)
        ax1.set_title('Linear Model')
        # ax1.text(0.98, 0.98, 'Linear\nModel', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes, fontsize='large')

        #"bright" color scheme from https://personal.sron.nl/~pault/
        # (not in order: blue, yellow, cyan, red, green, purple)
        # colors = ['#4477AA', '#CCBB44', '#66CCEE', '#EE6677', '#228833', '#AA3377']
        colors = ['k', '#CCBB44', '#66CCEE', '#EE6677', '#228833', '#AA3377']
        labels = ['10', '15', '20', '25', '30']
        thicks = [2.4, 1.2, 1.8, 1.2, 1.8]
        styles = ['solid', 'solid', 'solid', 'solid', 'solid']
        for i in range(len(labels)):
            ax1.plot(binned_lambdas, binned_alphas_1d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i], linewidth=thicks[i], linestyle=styles[i])
            if bootstrap:
                ax1.plot(binned_lambdas, bootstrap_binned_stds_thresh_1d[i], c = colors[i], drawstyle='steps', linestyle='--')
                #ax1.fill_between(binned_lambdas, bootstrap_binned_lower_thresh_1d[i], bootstrap_binned_upper_thresh_1d[i], linewidth=0.0, color=colors[i], alpha=0.5, step='pre') #TEMP
            else:
                ax1.plot(binned_lambdas, binned_stds_1d[i], c=colors[i], drawstyle='steps', linestyle='--')
            ax1.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
            ax1.set_ylabel(r"$\alpha_\lambda$ $\propto$ $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
            ax1.set_xlim(x_min, x_max)
            ax1.set_ylim(0, y_max)

        leg = ax1.legend(frameon=False, loc='upper left', ncol=2)

        ax1.xaxis.set_major_locator(MultipleLocator(1000))
        ax1.xaxis.set_minor_locator(MultipleLocator(200))
        ax1.yaxis.set_major_locator(MultipleLocator(0.1))
        ax1.yaxis.set_minor_locator(MultipleLocator(0.02))

        ax2 = fig.add_subplot(122)
        ax2.set_title('Nonlinear Model')
        # ax2.text(0.98, 0.98, 'Nonlinear\nModel', horizontalalignment='right', verticalalignment='top', transform=ax2.transAxes, fontsize='large')

        for i in range(len(labels)):
            ax2.plot(binned_lambdas, binned_alphas_2d[i], c=colors[i], drawstyle='steps', label=r'I$_{100} < %s$' % labels[i], linewidth=thicks[i], linestyle=styles[i])
            if bootstrap:
                ax2.plot(binned_lambdas, bootstrap_binned_stds_thresh_2d[i], c = colors[i], drawstyle='steps', linestyle='--')
                #ax2.fill_between(binned_lambdas, bootstrap_binned_lower_thresh_2d[i], bootstrap_binned_upper_thresh_2d[i], linewidth=0.0, color=colors[i], alpha=0.5, step='pre') #TEMP
            else:
                ax2.plot(binned_lambdas, binned_stds_2d[i], c=colors[i], drawstyle='steps', linestyle='--')
            ax2.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
            ax2.set_ylabel(r"$\alpha_\lambda^{\prime}$ $\propto$ $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
            ax2.set_xlim(x_min, x_max)
            ax2.set_ylim(0, y_max)

        leg = ax2.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, 0.15), ncol=2)

        ax2.xaxis.set_major_locator(MultipleLocator(1000))
        ax2.xaxis.set_minor_locator(MultipleLocator(200))
        ax2.yaxis.set_major_locator(MultipleLocator(0.1))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.02))

        if save != '0' and bootstrap:
            plt.savefig('../paper_figures/boss_thresholds_2panel_' + save + '.pdf', bbox_inches='tight')
            plt.clf()
        elif show_plots:
            plt.show()

    #############################
    # CALCULATIONS COMPARING NORTH AND SOUTH
    ################################
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
        std_north = bootstrap_stds[0]
        rel_std_north = std_north[(wavelength > 5000) & (wavelength < 8000)]
        rel_variance_north = rel_std_north**2
        avg_north_std = np.nanmean(rel_std_north)
        north_alphas = alphas[0]
        rel_north_alphas = north_alphas[(wavelength > 5000) & (wavelength < 8000)]
        numerator = np.nansum(np.divide(rel_north_alphas, rel_variance_north))
        denominator = np.nansum(np.divide(1, rel_variance_north))
        avg_north = numerator / denominator
        print("avg north:", avg_north)
        print(avg_north.shape)
        print("avg north std:", avg_north_std)

        integrand = np.multiply(north_alphas, wavelength_deltas)
        integrand_std = np.multiply(std_north**2, wavelength_deltas**2)
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
        std_south = bootstrap_stds[1]
        rel_std_south = std_south[(wavelength > 5000) & (wavelength < 8000)]
        rel_variance_south = rel_std_south**2
        avg_south_std = np.nanmean(rel_std_south)
        south_alphas = alphas[1]
        rel_south_alphas = south_alphas[(wavelength > 5000) & (wavelength < 8000)]
        numerator = np.nansum(np.divide(rel_south_alphas, rel_variance_south))
        denominator = np.nansum(np.divide(1, rel_variance_south))
        avg_south = numerator / denominator
        print("avg south:", avg_south)
        print("avg south std:", avg_south_std)

        integrand = np.multiply(south_alphas, wavelength_deltas)
        integrand_std = np.multiply(std_south**2, wavelength_deltas ** 2)
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

        # overall BOSS:
        rel_std_boss = bootstrap_stds[2][(wavelength > 5000) & (wavelength < 8000)]
        rel_variance_boss = rel_std_boss ** 2
        avg_boss_std = np.nanmean(rel_std_boss)
        rel_boss_alphas = alphas[2][(wavelength > 5000) & (wavelength < 8000)]
        numerator = np.nansum(np.divide(rel_boss_alphas, rel_variance_boss))
        denominator = np.nansum(np.divide(1, rel_variance_boss))
        avg_boss = numerator / denominator
        print("avg BOSS:", avg_boss)
        print("avg BOSS std:", avg_boss_std)

        #print results
        print("north avg:", avg_north)
        print("south avg:", avg_south)
        print("north idx:", north_color_idx_avg)
        print("south idx:", south_color_idx_avg)
        print("north ere:", north_ere, north_ere_std)
        print("south ere:", south_ere, south_ere_std)

        #calculate significance:
        #estimate variance:
        var_diff = avg_north_std**2 + avg_south_std**2
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

    if not boss and bootstrap:
        # find average and uncertainty for SDSS
        rel_std_sdss = bootstrap_stds[2][(wavelength > 5000) & (wavelength < 8000)]
        rel_variance_sdss = rel_std_sdss ** 2
        avg_sdss_std = np.nanmean(rel_std_sdss)
        rel_sdss_alphas = alphas[2][(wavelength > 5000) & (wavelength < 8000)]
        numerator = np.nansum(np.divide(rel_sdss_alphas, rel_variance_sdss))
        denominator = np.nansum(np.divide(1, rel_variance_sdss))
        avg_sdss = numerator / denominator
        print("avg SDSS:", avg_sdss)
        print("avg SDSS std:", avg_sdss_std)


##################################################

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

'''
#unbinned plots, SFD 1d vs 3 others
plot_emissions([0, 1], ["SFD", r"With $\tau$ Correction"], ['k', 'r'])
plt.show()
plot_emissions([0, 3], ["SFD", "With IRIS data"], ['k', 'r'])
plt.show()
plot_emissions([0, 2], ["SFD", r"With $\tau$ and IRIS"], ['k', 'r'])
plt.show()
'''






