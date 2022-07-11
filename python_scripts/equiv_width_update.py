#calculations to find equivalent widths, line ratios, temperatures

# bootstrapping: already set up to handle multiple spectra

# ISSUES:
# 5008 looks like double peak
# error prop from H corrections? (can get error on 4000 A break)
# investigate the integration warnings (now I am ignoring)
# continuum: should use weighted mean?

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy import interpolate
from scipy import integrate
from scipy import optimize
import sys #for command line args
from numpy.linalg import lstsq #for fitting gaussian
from scipy.integrate import quad
from tabulate import tabulate #for tables

import warnings
warnings.simplefilter("ignore")

#command line options
if len(sys.argv) != 4:
    print("Usage: equiv_width.py [boss: 0, 1] [bootstrap: 0] [loadkey]")
    exit(0)
boss = int(sys.argv[1])
bootstrap = int(sys.argv[2])
loadkey = sys.argv[3]

alpha_direc = '../alphas_and_stds/'
hdulist_direc = '../data/'  #'/Users/blakechellew/Documents/DustProject/BrandtFiles/'
    
# load wavelengths
if boss:
    wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')
else:
    hdulist = fits.open(hdulist_direc + 'SDSS_allskyspec.fits')
    wavelength = np.array(hdulist[1].data)

#load in alphas
if bootstrap:
    # full sky
    alphas = np.load('../alphas_and_stds/bootstrap_alphas_boss_iris_2d_' + loadkey + '_10.npy')
    alpha_stds = np.load('../alphas_and_stds/bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')
    # north
    # alphas = np.load('../alphas_and_stds/bootstrap_alphas_boss_iris_2d_north_' + loadkey + '.npy')
    # alpha_stds = np.load('../alphas_and_stds/bootstrap_alpha_stds_boss_iris_2d_north_' + loadkey + '.npy')
    # south
    # alphas = np.load('../alphas_and_stds/bootstrap_alphas_boss_iris_2d_south_' + loadkey + '.npy')
    # alpha_stds = np.load('../alphas_and_stds/bootstrap_alpha_stds_boss_iris_2d_south_' + loadkey + '.npy')
else:
    if boss:
        alphas = [np.load(alpha_direc + 'alphas_boss_iris_2d_' + loadkey + '_10.npy'), \
                  np.load(alpha_direc + 'alphas_boss_iris_2d_north_' + loadkey + '.npy'), \
                  np.load(alpha_direc + 'alphas_boss_iris_2d_south_' + loadkey + '.npy')]
        alpha_stds = [np.load(alpha_direc + 'alpha_stds_boss_iris_2d_' + loadkey + '_10.npy'), \
                      np.load(alpha_direc + 'alpha_stds_boss_iris_2d_north_' + loadkey + '.npy'), \
                      np.load(alpha_direc + 'alpha_stds_boss_iris_2d_south_' + loadkey + '.npy')]
    else:
        alphas = [np.load(alpha_direc + 'alphas_sdss_1d_' + loadkey + '.npy')]
        alpha_stds  = [np.load(alpha_direc + 'alpha_stds_sdss_1d_' + loadkey + '.npy')]

peaks = [4863, 4960, 5008, 5877, 6550, 6565, 6585, 6718, 6733]
left_ranges = [(4829, 4893), (4905, 4977), (4987, 5026), (5827, 5887), (6470, 6554), (6555, 6573), (6574.6, 6700), (6600, 6722), (6722, 6785)]
right_ranges = [None, None, None, None, (6600, 6700), (6600, 6700), None, (6745, 6785), (6600, 6700)]

idx_6565 = 5
idx_4863 = 0
idx_6585 = 6
idx_6550 = 4
idx_6718 = 7
idx_6733 = 8
if boss:
    # second o2 line. same continuum.
    peaks.insert(0, 3729.9)
    left_ranges.insert(0, (3700, 3772))
    right_ranges.insert(0, None)
    # first o2 line
    peaks.insert(0, 3727.1)
    left_ranges.insert(0, (3700, 3772))
    right_ranges.insert(0, None)

    idx_6565 += 2
    idx_4863 += 2
    idx_6585 += 2
    idx_6550 += 2
    idx_6718 += 2
    idx_6733 += 2

# get the line width at a particular wavelength based on BOSS resolution
# SDSS site: from R = 1560 at 3700 Angstroms to R= 2270 at 6000 Angstroms (blue channel),
#   and from R = 1850 at 6000 Angstroms to R = 2650 at 9000 Angstroms (reds channel)
def wav_to_width(wav):
    if wav < 6000:
        R1 = 1560
        R2 = 2270
        wav1 = 3700
        wav2 = 6000
    elif wav >= 6000:
        R1 = 1850
        R2 = 2650
        wav1 = 6000
        wav2 = 9000
    range_fraction = (wav - wav1) / (wav2 - wav1)
    R = R1 + (R2 - R1) * range_fraction
    fwhm = wav / R
    width_a = fwhm / 2
    return width_a
# xvals = np.linspace(3700, 9000, 500)
# plt.plot(xvals, [wav_to_width(v) for v in xvals])
# plt.show()
# exit(0)

# only for the 3-parameter fit with width, bc it's nonlinear
# do a 2-parameter fit and then calculate / return chi^2
def chi_sqr_width(width, *args):
    rel_alphas, rel_lambdas, rel_sigmas, peak_l = args
    a1, a2, _, _ = linear_fit(rel_alphas, rel_lambdas, rel_sigmas, peak_l, width)
    chi_2 = np.sum((gaussian_func(rel_lambdas, a1, a2, peak_l, width) - rel_alphas)**2 / rel_sigmas**2)
    return chi_2

#gaussian for fitting
def gaussian_func(x, a1, a2, a3, a4):
        return a1 + a2* np.exp(-np.power(x-a3, 2)/(2*a4**2))

# helper for fitting gaussian
# argument order: continuum, height, center, width
# is is a set of wavelengths
def width_helper(x, a1, a2, a3, a4):
    return (gaussian_func(x, a1, a2, a3, a4) - a1) / a1

#overall peak fitting
def linear_fit(rel_alphas, rel_lambdas, rel_sigmas, a3, a4):

    rel_vars = np.power(rel_sigmas, 2)
    rel_ivars = 1 / rel_vars
    y = np.exp(-np.power(rel_lambdas-a3, 2)/(2*a4**2))
    A = [[np.sum(rel_ivars), np.sum(np.multiply(y, rel_ivars))], \
         [np.sum(np.multiply(y, rel_ivars)), np.sum(np.multiply(np.power(y, 2), rel_ivars))]]
    b = [np.sum(np.multiply(rel_alphas, rel_ivars)), np.sum(np.multiply(np.multiply(rel_alphas, y), rel_ivars))]
    
    a1, a2 = lstsq(A, b)[0]

    #compute errors:
    bp = np.concatenate((rel_ivars.reshape(1, len(rel_ivars)), np.multiply(y, rel_ivars).reshape(1, len(rel_ivars))), axis=0) #b prime (derivative)
    Abp = np.dot(np.linalg.inv(A), bp)
    Abp_sqr = np.power(Abp, 2)
    a_vars = np.dot(Abp_sqr, rel_vars)
    a_stds = np.sqrt(a_vars)

    return a1, a2, a_stds[0], a_stds[1]


# fit the o2 peak with two gaussians, fixed wavelength
def linear_fit_two_fixed(rel_alphas, rel_lambdas, rel_sigmas, a3, a4, a6, a7):
    rel_vars = np.power(rel_sigmas, 2)
    rel_ivars = 1 / rel_vars
    y = np.exp(-np.power(rel_lambdas - a3, 2) / (2 * a4 ** 2)) + 1.5 * np.exp(-np.power(rel_lambdas - a6, 2) / (2 * a7 ** 2))
    A = [[np.sum(rel_ivars), np.sum(np.multiply(y, rel_ivars))], \
         [np.sum(np.multiply(y, rel_ivars)), np.sum(np.multiply(np.power(y, 2), rel_ivars))]]
    b = [np.sum(np.multiply(rel_alphas, rel_ivars)), np.sum(np.multiply(np.multiply(rel_alphas, y), rel_ivars))]

    a1, a2 = lstsq(A, b)[0]

    # compute errors:
    bp = np.concatenate((rel_ivars.reshape(1, len(rel_ivars)), np.multiply(y, rel_ivars).reshape(1, len(rel_ivars))),
                        axis=0)  # b prime (derivative)
    Abp = np.dot(np.linalg.inv(A), bp)
    Abp_sqr = np.power(Abp, 2)
    a_vars = np.dot(Abp_sqr, rel_vars)
    a_stds = np.sqrt(a_vars)

    return a1, a2  # skip stds for the moment

# for the doublet
def linear_fit_two_peaks(rel_alphas, rel_lambdas, rel_sigmas, a3, a4, a6, a7):
    rel_vars = np.power(rel_sigmas, 2)
    rel_ivars = 1 / rel_vars
    y1 = np.exp(-np.power(rel_lambdas - a3, 2) / (2 * a4 ** 2))
    y2 = np.exp(-np.power(rel_lambdas - a6, 2) / (2 * a7 ** 2))
    A = [[np.sum(rel_ivars), np.sum(y1 * rel_ivars), np.sum(y2 * rel_ivars)],
         [np.sum(y1 * rel_ivars), np.sum(y1 * y1 * rel_ivars), np.sum(y1 * y2 * rel_ivars)],
         [np.sum(y2 * rel_ivars), np.sum(y1 * y2 * rel_ivars), np.sum(y2 * y2 * rel_ivars)]]
    b = [np.sum(rel_alphas * rel_ivars), np.sum(rel_alphas * y1 * rel_ivars), np.sum(rel_alphas * y2 * rel_ivars)]

    a1, a2, a5 = lstsq(A, b)[0]
    return a1, a2, a5

    # skipping uncertainties bc finding with bootstrap
    
#fitting a line + gaussian to peak, calculate equiv width
def equiv_width(peak_l, alphas, alpha_stds, range1, range2=None, two_peaks=False, peak_l_2=None):

    #get indices of continuum boundaries
    range1_idx1 = np.argmin(np.abs(wavelength-range1[0]))
    range1_idx2 = np.argmin(np.abs(wavelength-range1[1]))
    if range2 is not None:
        range2_idx1 = np.argmin(np.abs(wavelength-range2[0]))
        range2_idx2 = np.argmin(np.abs(wavelength-range2[1]))
    else:
        range2_idx1, range2_idx2 = (0, 0)
    
    range_indices = np.concatenate((np.arange(range1_idx1, range1_idx2), np.arange(range2_idx1, range2_idx2)))
    rel_alphas = alphas[range_indices]
    rel_sigmas = alpha_stds[range_indices]
    rel_lambdas = wavelength[range_indices]

    peak_width = wav_to_width(peak_l)
    # peak_width = 1.85 #1.85 width from H-alpha

    if not two_peaks:
        a1, a2, a1_std, a2_std = linear_fit(rel_alphas, rel_lambdas, rel_sigmas, peak_l, peak_width)

        # 3 param fit for H-alpha width (only if not bootstrapping bc would take long time)
        if not bootstrap:
            best_fit_width = optimize.minimize(chi_sqr_width, 4, args=(rel_alphas, rel_lambdas, rel_sigmas, peak_l)).x
            print("best fit width for wavelength" + str(peak_l))
            print(best_fit_width)
            # a1, a2, _, _ = linear_fit(rel_alphas, rel_lambdas, rel_sigmas, peak_l, best_fit_width)
            # plt.plot(rel_lambdas, rel_alphas, '.')
            # x = np.linspace(6530, 6590, 300)
            # plt.plot(x, gaussian_func(x, a1, a2, peak_l, best_fit_width))
            # plt.show()
            # plt.xlim(6530, 6590)

        #integrate (over 10 sigma)
        num_sigmas = 10
        width, _ = quad(width_helper, peak_l-num_sigmas*peak_width, peak_l + num_sigmas*peak_width, args=(a1, a2, peak_l, peak_width))
        #error:
        gauss_integral, _ = quad(width_helper, peak_l-num_sigmas*peak_width, peak_l + num_sigmas*peak_width, args=(1, 1, peak_l, peak_width))
        frac_err = np.sqrt(np.power(a1_std/a1, 2) + np.power(a2_std/a2, 2))
        err = frac_err*(a2/a1)*gauss_integral

        # plot the fitted gaussian:
        """
        x_range = np.arange(range1[0], range1[1], .01)
        y = np.exp(-np.power(x_range-peak_l, 2)/(2*peak_width**2))
        alpha_pred = a1+a2*y
        plt.plot(x_range, alpha_pred, '.')
        plt.plot(rel_lambdas, rel_alphas, 'r.')
        plt.show()
        """

        return width, err
    else:
        # a1, a2, a5 = linear_fit_two_peaks(rel_alphas, rel_lambdas, rel_sigmas, peak_l, peak_width, peak_l_2, peak_width)

        peak_width_1 = wav_to_width(peak_l)
        peak_width_2 = wav_to_width(peak_l_2)

        a1, a2 = linear_fit_two_fixed(rel_alphas, rel_lambdas, rel_sigmas, peak_l, peak_width_1, peak_l_2, peak_width_2)

        # integrate (over 10 sigma)
        num_sigmas = 10
        width1, _ = quad(width_helper, peak_l - num_sigmas * peak_width, peak_l + num_sigmas * peak_width,
                        args=(a1, a2, peak_l, peak_width))
        width2, _ = quad(width_helper, peak_l_2 - num_sigmas * peak_width, peak_l_2 + num_sigmas * peak_width,
                         args=(a1, 1.5 * a2, peak_l_2, peak_width))
        # again, skipping error bc it's found with bootstrap


        # plot the fitted gaussians:
        # TEST: hard code values
        # a2 = 0.2 / 1.5
        # a5 = 0.2
        # a2 = 0.4 / 1.5
        # a5 = 0.4
        """
        x_range = np.arange(range1[0], range1[1], .01)
        y1 = np.exp(-np.power(x_range - peak_l, 2)/(2*peak_width**2))
        y2 = np.exp(-np.power(x_range - peak_l_2, 2) / (2 * peak_width ** 2))
        alpha_pred = a1+a2*y1+1.5*a2*y2
        plt.plot(x_range, alpha_pred, '.')
        # plot the peaks separately
        plt.plot(x_range, a1+a2*y1, '.', label='size:' + str(a2))
        plt.plot(x_range, a1+1.5*a2*y2, '.', label='size:' + str(1.5*a2))
        plt.plot(rel_lambdas, rel_alphas, 'r.')
        plt.vlines([peak_l, peak_l_2], 0, 2)
        plt.legend()
        plt.show()
        """


        if not bootstrap:
            print("o2 widths:")
            print(width1, width2)


        return width1 + width2  # combined width

#return H-alpha and H-beta ratios, with errors
def get_ratios(width, err, halpha, a_err, hbeta, b_err):
    a_ratio = width/halpha
    b_ratio = width/hbeta
    a_ratio_err = a_ratio*(a_err/halpha + err/width)
    b_ratio_err = b_ratio*(b_err/hbeta + err/width)
    return a_ratio, a_ratio_err, b_ratio, b_ratio_err

# calculate 4000 A break
def break_4000(alphas, stds, wav):

    idx_4000 = np.argmin(np.abs(wav-4000))
    idx_3850 = np.argmin(np.abs(wav-3850))
    idx_4150 = np.argmin(np.abs(wav - 4150))
    left_break = np.sum(np.diff(wav)[idx_3850:idx_4000] * alphas[idx_3850:idx_4000] / wav[idx_3850:idx_4000])
    right_break = np.sum(np.diff(wav)[idx_4000:idx_4150] * alphas[idx_4000:idx_4150] / wav[idx_4000:idx_4150])

    # TEST
    # idx_3950 = np.argmin(np.abs(wav-3950))
    # idx_4100 = np.argmin(np.abs(wav-4100))
    # left_break = np.sum(np.diff(wav)[idx_3850:idx_3950] * alphas[idx_3850:idx_3950] / wav[idx_3850:idx_3950])
    # right_break = np.sum(np.diff(wav)[idx_4000:idx_4100] * alphas[idx_4000:idx_4100] / wav[idx_4000:idx_4100])

    delta = left_break / right_break
    beta_width = -2.2 * delta + .17
    alpha_width = -1.5 * delta - .19

    left_break_err = np.sqrt(np.sum(np.power(stds[idx_3850:idx_4000], 2)))
    right_break_err = np.sqrt(np.sum(np.power(stds[idx_4000:idx_4150], 2)))
    frac_delta_err = np.sqrt((left_break_err/left_break)**2 + (right_break_err/right_break)**2)
    delta_err = frac_delta_err*delta
    beta_err = 2.2 * delta_err
    alpha_err = 1.5 * delta_err
    
    return delta, alpha_width, beta_width, alpha_err, beta_err

# calculate 4000 A break:
deltas = np.zeros(len(alphas))
bws = np.zeros(len(alphas))
aws = np.zeros(len(alphas))
bw_errs = np.zeros(len(alphas))
aw_errs = np.zeros(len(alphas))
for i in range(len(alphas)):
    deltas[i], aws[i], bws[i], aw_errs[i], bw_errs[i] = break_4000(alphas[i], alpha_stds[i], wavelength)
    if not bootstrap:
        print("4000 A break:")
        print(deltas[i])

# calculate equivalent widths:
#width and error for each wavelength
#n ratio, err
#temp ratio and err for N and S

if bootstrap:
    bootstrap_widths = np.zeros((len(peaks), alphas.shape[0]))

for i in range(len(alphas)):

    if i % 100 == 0:
        print("PROGRESS:", i)

    widths = np.zeros((len(peaks), 1))
    width_errs = np.zeros((len(peaks), 1))
    ratios = np.zeros((len(peaks), 4))

    for j, peak_l in enumerate(peaks):
        widths[j], width_errs[j] = equiv_width(peaks[j], alphas[i], alpha_stds[i], left_ranges[j], right_ranges[j])
        if j == idx_6565:
            widths[j] = widths[j] - aws[i]
            width_errs[j] = np.sqrt(width_errs[j]**2 + aw_errs[i]**2)
        if j == idx_4863:
            widths[j] = widths[j] - bws[i]
            width_errs[j] = np.sqrt(width_errs[j]**2 + bw_errs[i]**2)
    for j, peak_l in enumerate(peaks):
        ratios[j] = np.array(get_ratios(widths[j], width_errs[j], widths[idx_6565], width_errs[idx_6565], widths[idx_4863], width_errs[idx_4863])).T[0]
    # find the o2 doublet and overwrite
    o2_combined = equiv_width(peaks[0], alphas[i], alpha_stds[i], left_ranges[0], right_ranges[0], two_peaks=True, peak_l_2=peaks[1])
    widths[0], widths[1] = o2_combined, o2_combined

    if bootstrap:
        bootstrap_widths[:, i] = widths.ravel()
    else:
        #N ratio: 6585/6550
        N_ratio = widths[idx_6585] / widths[idx_6550]
        N_err = N_ratio*(width_errs[idx_6585]/widths[idx_6585] + width_errs[idx_6550]/widths[idx_6550])
        #ISM:
        temp_ratio_N = (widths[idx_6585] + widths[idx_6550])/widths[idx_6565]
        temp_ratio_N_err = temp_ratio_N * (np.sqrt(pow(width_errs[idx_6585], 2) + pow(width_errs[idx_6550], 2))/(widths[idx_6585] + \
                            widths[idx_6550]) + (width_errs[idx_6565]/widths[idx_6565]))
        temp_ratio_S = (widths[idx_6718] + widths[idx_6733])/widths[idx_6565]
        temp_ratio_S_err = temp_ratio_S * (np.sqrt(pow(width_errs[idx_6718], 2) + pow(width_errs[idx_6733], 2))/(widths[idx_6718] + \
                            widths[idx_6733]) + (width_errs[idx_6565]/widths[idx_6565]))

        #tables:
        titles = ["4863 (H-beta)", "4960 (OIII)", "5008 (OIII)", "5877 (HeI)", "6550 (NII)", "6565 (H-alpha)", "6585 (NII)", "6718 (SII)", "6733 (SII)"]
        if boss:
            titles.insert(0, "3729 (OII)")
            titles.insert(0, "3726 (OII)")
        titles = np.reshape(titles, (len(titles), 1))
        headers=["Wavelength", "Width", "Err", "H-alpha ratio", "Err", "H-beta ratio", "Err"]

        cell_text = np.concatenate((titles, widths, width_errs, ratios), axis=1)
        print("\n" + tabulate(cell_text, headers))
        print("N ratio:", N_ratio, "+/-", N_err)
        print("N to alpha: from", temp_ratio_N - temp_ratio_N_err, "to", temp_ratio_N + temp_ratio_N_err)
        print("S to alpha: from", temp_ratio_S - temp_ratio_S_err, "to", temp_ratio_S + temp_ratio_S_err)


    '''
    #plot the widths with error bars:
    if boss:
        widths = widths[1:, :]
        width_errs = width_errs[1:, :]
    plt.errorbar(np.arange(widths.shape[0]), widths, fmt='.', yerr=width_errs)
    plt.show()
    '''

if bootstrap:
    # find percentiles for the widths
    print(bootstrap_widths)
    lower_bounds = np.percentile(bootstrap_widths, 16, axis=1)
    upper_bounds = np.percentile(bootstrap_widths, 84, axis=1)
    errors = (upper_bounds - lower_bounds) / 2
    print(errors)

    # and for 4000 A break:
    lower_delta = np.percentile(deltas, 16)
    upper_delta = np.percentile(deltas, 84)
    error_delta = (upper_delta - lower_delta) / 2
    print("delta error:", error_delta)

    all_percentiles = np.percentile(deltas, [0.15, 2.5, 16, 50, 84, 97.5, 99.85])
    print("delta percentiles [0.15, 2.5, 16, 50, 84, 97.5, 99.85]")
    print(all_percentiles)

# find 4000 A break for models:
cst_coeff = 0.3  # TEST, was 0.5
cst_coeff_wd = 0.6
load_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/convergence_test/'
wd_t5e9_filename = 't5e9_grid_150_wd.npy'
zd_t5e9_filename = 't5e9_grid_150_zd.npy'
wd_t9e9_filename = 't9e9_grid_150_wd.npy'
zd_t9e9_filename = 't9e9_grid_150_zd.npy'
wd_cst_filename = 'cst_grid_150_wd.npy'
zd_cst_filename = 'cst_grid_150_zd.npy'
wd01_t5e9_alphas = np.load(load_path + wd_t5e9_filename)
zd01_t5e9_alphas = np.load(load_path + zd_t5e9_filename)
wd01_t9e9_alphas = np.load(load_path + wd_t9e9_filename)
zd01_t9e9_alphas = np.load(load_path + zd_t9e9_filename)
wd01_cst_alphas = np.load(load_path + wd_cst_filename)
zd01_cst_alphas = np.load(load_path + zd_cst_filename)
wd01_combo = np.load(load_path + wd_t5e9_filename) + cst_coeff_wd * np.load(load_path + wd_cst_filename)
zd01_combo = np.load(load_path + zd_t5e9_filename) + cst_coeff * np.load(load_path + zd_cst_filename)
break_zd_t5e9, _, _, _, _ = break_4000(zd01_t5e9_alphas, np.ones(len(zd01_t5e9_alphas)), wavelength)
break_zd_combo, _, _, _, _ = break_4000(zd01_combo, np.ones(len(zd01_t5e9_alphas)), wavelength)
break_wd_t5e9, _, _, _, _ = break_4000(wd01_t5e9_alphas, np.ones(len(zd01_t5e9_alphas)), wavelength)
break_wd_combo, _, _, _, _ = break_4000(wd01_combo, np.ones(len(zd01_t5e9_alphas)), wavelength)
print("t5e9 ZD:", break_zd_t5e9)
print("combo ZD:", break_zd_combo)
print("t5e9 WD:", break_wd_t5e9)
print("combo WD:", break_wd_combo)



