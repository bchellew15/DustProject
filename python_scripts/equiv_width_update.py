# calculations to find equivalent widths, line ratios, temperatures
# also plot spectra (overall and certain sections) for comparison

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

import warnings
warnings.simplefilter("ignore")

#command line options
if len(sys.argv) != 3:
    print("Usage: equiv_width.py [boss: 0, 1] [plots_on: 0, 1]")
    exit(0)
    
boss = int(sys.argv[1])
plots_on = int(sys.argv[2])
    
# load wavelengths
if boss:
    wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')
else:
    hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
    wavelength = np.array(hdulist[1].data)


# load in npy files
# original, tao, tao AND iris, iris
if boss:
    alphas = [np.load('../alphas_and_stds/alphas_1d_boss.npy'), np.load('../alphas_and_stds/alphas_2d_boss.npy'), \
              np.load('../alphas_and_stds/alphas_iris_boss.npy'), np.load('../alphas_and_stds/alphas_iris_1d_boss.npy')]
    alpha_stds = [np.load('../alphas_and_stds/alphas_1d_boss_stds.npy'), np.load('../alphas_and_stds/alphas_2d_boss_stds.npy'), \
                  np.load('../alphas_and_stds/alphas_iris_boss_stds.npy'), np.load('../alphas_and_stds/alphas_iris_1d_boss_stds.npy')]
else:
    alphas = [np.load('../alphas_and_stds/alphas_1d.npy'), np.load('../alphas_and_stds/alphas_2d.npy'), \
              np.load('../alphas_and_stds/alphas_iris.npy'), np.load('../alphas_and_stds/alphas_iris_1d.npy')]
    alpha_stds = [np.load('../alphas_and_stds/alphas_1d_stds.npy'), np.load('../alphas_and_stds/alphas_2d_stds.npy'), \
                  np.load('../alphas_and_stds/alphas_iris_stds.npy'), np.load('../alphas_and_stds/alphas_iris_stds_1d.npy')]
num_arrays = len(alphas)
num_regions = 8 #I think this is number of peaks to look at

# integrate using 3 wavelength elements on either side
def equiv_width2(peak_l, alphas, alpha_stds, cont, stdev):

    #find indices of 2 closest wavelengths
    peak_idx1 = np.argpartition(np.abs(wavelength-peak_l), 0)[0]
    peak_idx2 = np.argpartition(np.abs(wavelength-peak_l), 1)[1]
    left_peak_idx = min(peak_idx1, peak_idx2)
    
    bars = (alphas-cont)/cont
    width = np.sum(bars[left_peak_idx-2:left_peak_idx+4]) 

    # each bar has error bar*sigma*(1/peak + 1/continuum)
    # then add in quadrature
    # (because we can simplify to alpha/continuum - 1)
    bar_errs = np.array([bars[i]*(alpha_stds[i]/alphas[i] + stdev/cont) for i in range(left_peak_idx-2, left_peak_idx+4)])
    err = np.sqrt(np.sum(np.power(bar_errs, 2)))

    return width, err

#return H-alpha and H-beta ratios, with errors
def get_ratios(width, err, halpha, a_err, hbeta, b_err):
    a_ratio = width/halpha
    b_ratio = width/hbeta
    a_ratio_err = a_ratio*(a_err/halpha + err/width)
    b_ratio_err = b_ratio*(b_err/hbeta + err/width)
    return a_ratio, a_ratio_err, b_ratio, b_ratio_err

#calculate continuum and stds
continuum_6600s = np.zeros(num_arrays)
continuum_6600_stds = np.zeros(num_arrays)
continuum_4800s = np.zeros(num_arrays)
continuum_4800_stds = np.zeros(num_arrays)
continuum_o3s = np.zeros(num_arrays)
continuum_o3_stds = np.zeros(num_arrays)
for i in range(num_arrays): 
    # noise for 6600 range:
    alphas_in_range = alphas[i][(wavelength>6600) * (wavelength < 6700)]
    continuum_6600s[i] = np.mean(alphas_in_range)
    continuum_6600_stds[i] = np.mean(alpha_stds[i][(wavelength>6600) * (wavelength < 6700)])
    #np.std(alphas_in_range, ddof=1)
    # noise for H-beta
    alphas_in_range = alphas[i][(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
    continuum_4800s[i] = np.mean(alphas_in_range)
    continuum_4800_stds[i] = np.mean(alpha_stds[i][(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)])
    #np.std(alphas_in_range, ddof=1)
    # noise for OIII
    alphas_in_range = alphas[i][(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
    continuum_o3s[i] = np.mean(alphas_in_range)
    continuum_o3_stds[i] = np.mean(alpha_stds[i][(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)])
    #np.std(alphas_in_range, ddof=1)

# calculate 4000 A break
def break_4000(alphas):
    alpha_func = interpolate.interp1d(wavelength, alphas)
    left_break = integrate.quad(alpha_func, 3850, 4000)[0]
    right_break = integrate.quad(alpha_func, 4000, 4150)[0]
    delta = left_break / right_break
    beta_width = -2.2*delta + .17
    alpha_width = -1.5*delta - .19
    return alpha_width, beta_width

# calculate 4000 A break:
bws = np.zeros(num_arrays)
aws = np.zeros(num_arrays)
for i in range(num_arrays):
    bws[i], aws[i] = break_4000(alphas[i])

# calculate equivalent widths:

width_4863s = np.zeros(num_arrays)
err_4863s = np.zeros(num_arrays)
width_4960s = np.zeros(num_arrays)
err_4960s = np.zeros(num_arrays)
width_5008s = np.zeros(num_arrays)
err_5008s = np.zeros(num_arrays)
width_6550s = np.zeros(num_arrays)
err_6550s = np.zeros(num_arrays)
width_6565s = np.zeros(num_arrays)
err_6565s = np.zeros(num_arrays)
width_6585s = np.zeros(num_arrays)
err_6585s = np.zeros(num_arrays)
width_6718s = np.zeros(num_arrays)
err_6718s = np.zeros(num_arrays)
width_6733s = np.zeros(num_arrays)
err_6733s = np.zeros(num_arrays)
N_ratios = np.zeros(num_arrays)
N_errs = np.zeros(num_arrays)
temp_ratio_Ns = np.zeros(num_arrays)
temp_ratio_N_errs = np.zeros(num_arrays)
temp_ratio_Ss = np.zeros(num_arrays)
temp_ratio_S_errs = np.zeros(num_arrays)
widths = [] #array of width arrays
errs = [] #array of err arrays
ratios = np.zeros((num_arrays, num_regions, 4)) #array of ratio arrays

for i in range(num_arrays):
    # 4863 peak (H-beta)
    width_4863s[i], err_4863s[i] = equiv_width2(4863, alphas[i], alpha_stds[i], continuum_4800s[i], continuum_4800_stds[i])
    width_4863s[i] = width_4863s[i] - bws[i] #correction
    width_4960s[i], err_4960s[i] = equiv_width2(4960, alphas[i], alpha_stds[i], continuum_o3s[i], continuum_o3_stds[i])
    width_5008s[i], err_5008s[i] = equiv_width2(5008, alphas[i], alpha_stds[i], continuum_o3s[i], continuum_o3_stds[i])
    width_6550s[i], err_6550s[i] = equiv_width2(6550, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
    # 6565 peak (H-alpha)
    width_6565s[i], err_6565s[i] = equiv_width2(6565, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
    width_6565s[i] = width_6565s[i] - aws[i] #correction
    width_6585s[i], err_6585s[i] = equiv_width2(6585, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
    width_6718s[i], err_6718s[i] = equiv_width2(6718, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
    width_6733s[i], err_6733s[i] = equiv_width2(6733, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
    width_arr = np.array([width_4863s[i], width_4960s[i], width_5008s[i], width_6550s[i], \
                          width_6565s[i], width_6585s[i], width_6718s[i], width_6733s[i]])
    widths.append(width_arr.reshape(width_arr.shape[0],1))
    err_arr = np.array([err_4863s[i], err_4960s[i], err_5008s[i], err_6550s[i], \
                        err_6565s[i], err_6585s[i], err_6718s[i],  err_6733s[i]])
    errs.append(err_arr.reshape(err_arr.shape[0], 1))
    ratios[i] = np.array([get_ratios(widths[i][j], errs[i][j], width_6565s[i], err_6565s[i], width_4863s[i], err_4863s[i]) \
                          for j in range(num_regions)])[:,:,0]
    #N ratio: 6585/6550
    N_ratios[i] = width_6585s[i] / width_6550s[i]
    N_errs[i] = N_ratios[i]*(err_6585s[i]/width_6585s[i] + err_6550s[i]/width_6550s[i])
    #ISM:
    temp_ratio_Ns[i] = (width_6585s[i] + width_6550s[i])/width_6565s[i]
    temp_ratio_N_errs[i] = temp_ratio_Ns[i] * (np.sqrt(pow(err_6585s[i], 2) + \
                        pow(err_6550s[i], 2))/(width_6585s[i] + width_6550s[i]) + (err_6565s[i]/width_6565s[i]))
    temp_ratio_Ss[i] = (width_6718s[i] + width_6733s[i])/width_6565s[i]
    temp_ratio_S_errs[i] = temp_ratio_Ss[i] * (np.sqrt(pow(err_6718s[i], 2) + pow(err_6733s[i], 2))/(width_6718s[i] + \
                        width_6733s[i]) + (err_6565s[i]/width_6565s[i]))


#tables:
from tabulate import tabulate
titles = ["4863 (H-beta)", "4960 (OIII)", "5008 (OIII)", "6550 (NII)", "6565 (H-alpha)", "6585 (NII)", "6718 (SII)", "6733 (SII)"]
titles = np.reshape(titles, (len(titles), 1))
headers=["Wavelength", "Width", "Err", "H-alpha ratio", "Err", "H-beta ratio", "Err"]

for i in range(num_arrays):
    cell_text = np.concatenate((titles, widths[i], errs[i], ratios[i]), axis=1)
    print("\n" + tabulate(cell_text, headers))
    print("N ratio:", N_ratios[i], "+/-", N_errs[i])
    print("N to alpha: from", temp_ratio_Ns[i] - temp_ratio_N_errs[i], "to", temp_ratio_Ns[i] + temp_ratio_N_errs[i])
    print("S to alpha: from", temp_ratio_Ss[i] - temp_ratio_S_errs[i], "to", temp_ratio_Ss[i] + temp_ratio_S_errs[i])

#plot unbinned spectra
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

#plot unbinned spectra:
if plots_on:
    plot_emissions(alphas[0], alphas[1], alpha_stds[0], alpha_stds[1], "SFD", r"With $\tau$ Correction")
    plt.show()
    plot_emissions(alphas[0], alphas[3], alpha_stds[0], alpha_stds[3], "SFD", "With IRIS data")
    plt.show()
    plot_emissions(alphas[0], alphas[2], alpha_stds[0], alpha_stds[2], "SFD", r"With $\tau$ and IRIS" )
    plt.show()

#remove this?
'''
from scipy import stats
bin_means, bin_edges, binnumber = stats.binned_statistic(wavelength, alphas1, statistic='mean', bins=108)
bin_width = (bin_edges[1] - bin_edges[0])
print("bin width: ", bin_width)
print(bin_edges[0], bin_edges[-1])
plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g')
if plots_on:
    plt.show()
'''

#plot binned alpha vs wavelength (original) 
binned_lambdas = np.linspace(wavelength[0], wavelength[-1], 109) #108 bins
binned_lambdas = np.arange(3900, 9200, 50)
binned_alphas = []
binned_stds = []

for i in range(num_arrays):
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

plt.figure(figsize=(6, 14))
    
plt.subplot(3, 1, 1)    
#compare original to tao
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[1], c='r', drawstyle='steps', label=r'With $\tau$ Correction')
plt.plot(binned_lambdas, binned_stds[1], c='r', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(3850, 9200)
plt.ylim(0, 0.4)
plt.legend(frameon=False)

plt.subplot(3, 1, 2)
#compare original to iris
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[3], c='r', drawstyle='steps', label='With IRIS Data')
plt.plot(binned_lambdas, binned_stds[3], c='r', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(3850, 9200)
plt.ylim(0, 0.4)
plt.legend(frameon=False)

plt.subplot(3, 1, 3)
#compare original to tao AND iris
plt.plot(binned_lambdas, binned_alphas[0], c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_stds[0], c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas[2], c='r', drawstyle='steps', label=r'With IRIS and $\tau$')
plt.plot(binned_lambdas, binned_stds[2], c='r', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(3850, 9200)
plt.ylim(0, 0.4)
plt.legend(frameon=False)

plt.tight_layout()
#plt.savefig("vert_3_panels.png")
if plots_on:
    plt.show()

#all 4 together:
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
plt.xlim(3850, 9200)
plt.ylim(0, 0.4)
plt.legend(frameon=False)
if plots_on:
    plt.show()


'''
#range1 and range2 for noise
def equiv_width(peak_l, alphas, range1, range2=(0, 0)):
    #look at noise:
    alphas_in_range = alphas[(wavelength>range1[0])*(wavelength<range1[1]) \
                             +(wavelength>range2[0])*(wavelength<range2[1])]
    mean_val = np.mean(alphas_in_range)
    std = np.std(alphas_in_range)
    #integrate part above line:
    alpha_func = interpolate.interp1d(wavelength, (alphas-mean_val)/mean_val)
    #find bounds of peak:
    peak_idx = np.argmin(np.abs(wavelength-peak_l))
    u_idx = peak_idx
    while alphas[u_idx] > mean_val:
        u_idx += 1
    u_bound = wavelength[u_idx]
    l_idx = peak_idx
    while alphas[l_idx] > mean_val:
        l_idx -= 1
    l_bound = wavelength[l_idx]
    integral_result = integrate.quad(alpha_func, l_bound, u_bound)
    return integral_result, mean_val
'''

#analyze wavelength array:
'''
#min and max:
min_l = min(wavelength)
max_l =  max(wavelength)
print("interval:", (max_l-min_l)/4000) #interval: 1.35 angstroms

#intervals in wavelength array:
wavelength_diffs = [wavelength[i+1] - wavelength[i] for i in range(len(wavelength)-1)]
plt.hist(wavelength_diffs)
plt.show()
'''


