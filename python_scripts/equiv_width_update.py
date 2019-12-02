#calculations to find equivalent widths, line ratios, temperatures

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
import sys #for command line args
from scipy.optimize import curve_fit #for fitting gaussian to equiv width
from numpy.linalg import lstsq #for fitting gaussian
from scipy.integrate import quad
from tabulate import tabulate #for tables

import warnings
warnings.simplefilter("ignore")

#command line options
if len(sys.argv) != 2:
    print("Usage: equiv_width.py [boss: 0, 1]")
    exit(0)  
boss = int(sys.argv[1])
    
# load wavelengths
if boss:
    wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')
else:
    hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
    wavelength = np.array(hdulist[1].data)

#load in alphas
#boss: iris, 2d
#sdss: replica of brandt paper (see end of file)
if boss:
    alphas = [np.load('../alphas_and_stds/alphas_boss_iris_91119_10.npy')]
    alpha_stds = [np.load('../alphas_and_stds/alpha_stds_boss_iris_91119_10.npy')]
else:
    alphas = [np.load('../alphas_and_stds/alphas_91019_10.npy'), \
              np.load('../alphas_and_stds/alphas_sdss_iris_2d_102019.npy')]
    alpha_stds = [np.load('../alphas_and_stds/alpha_stds_91019_10.npy'), \
                  np.load('../alphas_and_stds/alpha_stds_sdss_iris_2d_102019.npy')]

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
    peaks.insert(0, 3727)
    left_ranges.insert(0, (3700, 3772))
    right_ranges.insert(0, None)
    idx_6565 += 1
    idx_4863 += 1
    idx_6585 += 1
    idx_6550 += 1
    idx_6718 += 1
    idx_6733 += 1

#gaussian for fitting
def gaussian_func(x, a1, a2, a3, a4):
        return a1 + a2* np.exp(-np.power(x-a3, 2)/(2*a4**2))

#helper for fitting gaussian
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
    
#fitting a line + gaussian to peak, calculate equiv width
def equiv_width(peak_l, alphas, alpha_stds, range1, range2=None):

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
    
    peak_width = 1.85 #1.85 width from H-alpha
    a1, a2, a1_std, a2_std = linear_fit(rel_alphas, rel_lambdas, rel_sigmas, peak_l, peak_width) 

    #plot the fitted gaussian:
    #x_range = np.arange(range1[0], range1[1], .01)
    #y = np.exp(-np.power(x_range-peak_l, 2)/(2*peak_width**2))
    #alpha_pred = a1+a2*y
    #plt.plot(x_range, alpha_pred, '.')
    #plt.plot(rel_lambdas, rel_alphas, 'r.')
    #plt.show()
    
    #integrate (over 10 sigma)
    num_sigmas = 10
    width, _ = quad(width_helper, peak_l-num_sigmas*peak_width, peak_l + num_sigmas*peak_width, args=(a1, a2, peak_l, peak_width))
    #error:
    gauss_integral, _ = quad(width_helper, peak_l-num_sigmas*peak_width, peak_l + num_sigmas*peak_width, args=(1, 1, peak_l, peak_width))
    frac_err = np.sqrt(np.power(a1_std/a1, 2) + np.power(a2_std/a2, 2))
    err = frac_err*(a2/a1)*gauss_integral

    #left off here on 12/1/19
    
    #calculate errors:
    #continuum: remove indices of peak
    close_peak_idx = np.argmin(np.abs(wavelength-peak_l))
    for i in range(close_peak_idx-3, close_peak_idx+4):
        range_indices = np.delete(range_indices, np.argwhere(range_indices==i))
    rel_alphas = alphas[range_indices]
    rel_sigmas = alpha_stds[range_indices]
    rel_lambdas = wavelength[range_indices] #needed?
    
    cont = np.average(rel_alphas, weights = 1/np.power(rel_sigmas, 2))
    stdev = np.sqrt(np.sum(np.power(rel_sigmas, 2))) / len(rel_alphas)

    delta_lambda = wavelength[close_peak_idx+1] - wavelength[close_peak_idx]
    bars = (alphas-cont)/cont * delta_lambda

    #calculate widths the other way:
    #width = np.sum(bars[left_peak_idx+l_off:left_peak_idx+r_off])

    # each bar has error bar*sigma*(1/peak + 1/continuum)*delta_lambda
    # then add in quadrature (because we can simplify to alpha/continuum - 1)
    bar_errs = delta_lambda*(alphas/cont)*np.sqrt(np.power(np.divide(alpha_stds, alphas), 2) + (stdev/cont)**2)
    #err = np.sqrt(np.sum(np.power(bar_errs[close_peak_idx-3:close_peak_idx+4], 2)))

    return width, err

#return H-alpha and H-beta ratios, with errors
def get_ratios(width, err, halpha, a_err, hbeta, b_err):
    a_ratio = width/halpha
    b_ratio = width/hbeta
    a_ratio_err = a_ratio*(a_err/halpha + err/width)
    b_ratio_err = b_ratio*(b_err/hbeta + err/width)
    return a_ratio, a_ratio_err, b_ratio, b_ratio_err

# calculate 4000 A break
def break_4000(alphas):

    #temp:
    idx_4000 = np.argmin(np.abs(wavelength-4000))
    idx_3850 = np.argmin(np.abs(wavelength-3850))
    idx_4150 = np.argmin(np.abs(wavelength-4150))
    delta_lambda = wavelength[idx_4000+1] - wavelength[idx_4000]
    left_break = delta_lambda*np.sum(alphas[idx_3850:idx_4000])
    right_break = delta_lambda*np.sum(alphas[idx_4000:idx_4150])

    #alpha_func = interpolate.interp1d(wavelength, alphas)
    #left_break = integrate.quad(alpha_func, 3850, 4000)[0]
    #right_break = integrate.quad(alpha_func, 4000, 4150)[0]
    delta = left_break / right_break
    beta_width = -2.2*delta + .17
    alpha_width = -1.5*delta - .19
    return alpha_width, beta_width

# calculate 4000 A break:
bws = np.zeros(len(alphas))
aws = np.zeros(len(alphas))
for i in range(len(alphas)):
    aws[i], bws[i] = break_4000(alphas[i])

# calculate equivalent widths:
#width and error for each wavelength
#n ratio, err
#temp ratio and err for N and S

for i in range(len(alphas)):

    widths = np.zeros((len(peaks), 1))
    width_errs = np.zeros((len(peaks), 1))
    ratios = np.zeros((len(peaks), 4))
    
    for j, peak_l in enumerate(peaks):
        widths[j], width_errs[j] = equiv_width(peaks[j], alphas[i], alpha_stds[i], left_ranges[j], right_ranges[j])
        if j == idx_6565:
            widths[j] = widths[j] - aws[i]
        if j == idx_4863:
            widths[j] = widths[j] - bws[i]
    for j, peak_l in enumerate(peaks):
        ratios[j] = np.array(get_ratios(widths[j], width_errs[j], widths[idx_6565], width_errs[idx_6565], widths[idx_4863], width_errs[idx_4863])).T[0]

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
        titles.insert(0, "3727 (OII)")
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


