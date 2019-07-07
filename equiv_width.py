# load in alphas and make a scatterplot
# also do calculations to find equivalent widths, line ratios, temperatures, etc.

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

#option to turn plots off
plots_on = True
if len(sys.argv) == 2:
    if sys.argv[1] == "0":
        plots_on = False
    elif sys.argv[1] != 1:
        print("Usage: equiv_width.py [plots_on]")
        exit(0)


# load wavelengths
hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
wavelength = np.array(hdulist[1].data)  # Angstroms

# load in npy files
# original
alphas1 = np.load('alphas_1d.npy')
alpha_std1 = np.load('alphas_1d_stds.npy')
# tao
alphas2 = np.load('alphas_2d.npy')
alpha_std2 = np.load('alphas_2d_stds.npy')
# tao and iris
alphas3 = np.load('alphas_iris.npy')
alpha_std3 = np.load('alphas_iris_stds.npy')
# iris
alphas4 = np.load('alphas_iris_1d.npy')
alpha_std4 = np.load('alphas_iris_stds_1d.npy')

# integrate using 6-7 wavelength elements on either side
def equiv_width2(peak_l, alphas, alpha_stds, cont, stdev):
    #find indices of 2 closest wavelengths
    peak_idx1 = np.argpartition(np.abs(wavelength-peak_l), 0)[0]
    peak_idx2 = np.argpartition(np.abs(wavelength-peak_l), 1)[1]
    left_peak_idx = min(peak_idx1, peak_idx2)
    
    bars = (alphas-cont)/cont #this is width 7, not 6
    width = np.sum(bars[left_peak_idx-2:left_peak_idx+4]) 

    # each bar has error bar*sigma*(1/peak + 1/continuum)
    # then add in quadrature
    # (because we can simplify to alpha/continuum - 1)
    bar_errs = np.array([bars[i]*(alpha_stds[i]/alphas[i] + stdev/cont) for i in range(left_peak_idx-2, left_peak_idx+4)])
    err = np.sqrt(np.sum(np.power(bar_errs, 2)))

    return width, err

#return H-alpha and H-beta ratios, with errors
def ratios(width, err, halpha, a_err, hbeta, b_err):
    a_ratio = width/halpha
    b_ratio = width/hbeta
    a_ratio_err = a_ratio*(a_err/halpha + err/width)
    b_ratio_err = b_ratio*(b_err/hbeta + err/width)
    return a_ratio, a_ratio_err, b_ratio, b_ratio_err

# noise for alpha1:
alphas_in_range = alphas1[(wavelength>6600) * (wavelength < 6700)]
continuum_6600_1 = np.mean(alphas_in_range)
continuum_6600_std_1 = np.std(alphas_in_range, ddof=1)
# noise for H-beta
alphas_in_range = alphas1[(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
continuum_4800_1 = np.mean(alphas_in_range)
continuum_4800_std_1 = np.std(alphas_in_range, ddof=1)
# noise for OIII
alphas_in_range = alphas1[(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
continuum_o3_1 = np.mean(alphas_in_range)
continuum_o3_std_1 = np.std(alphas_in_range, ddof=1)

# noise for alpha2:
alphas_in_range = alphas2[(wavelength>6600) * (wavelength < 6700)]
continuum_6600_2 = np.mean(alphas_in_range)
continuum_6600_std_2 = np.std(alphas_in_range, ddof=1)
# noise for H-beta
alphas_in_range = alphas2[(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
continuum_4800_2 = np.mean(alphas_in_range)
continuum_4800_std_2 = np.std(alphas_in_range, ddof=1)
# noise for OIII
alphas_in_range = alphas2[(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
continuum_o3_2 = np.mean(alphas_in_range)
continuum_o3_std_2 = np.std(alphas_in_range, ddof=1)

# noise for alpha3:
alphas_in_range = alphas3[(wavelength>6600) * (wavelength < 6700)]
continuum_6600_3 = np.mean(alphas_in_range)
continuum_6600_std_3 = np.std(alphas_in_range, ddof=1)
# noise for H-beta
alphas_in_range = alphas3[(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
continuum_4800_3 = np.mean(alphas_in_range)
continuum_4800_std_3 = np.std(alphas_in_range, ddof=1)
# noise for OIII
alphas_in_range = alphas3[(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
continuum_o3_3 = np.mean(alphas_in_range)
continuum_o3_std_3 = np.std(alphas_in_range, ddof=1)

# noise for alpha4:
alphas_in_range = alphas4[(wavelength>6600) * (wavelength < 6700)]
continuum_6600_4 = np.mean(alphas_in_range)
continuum_6600_std_4 = np.std(alphas_in_range, ddof=1)
# noise for H-beta
alphas_in_range = alphas4[(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
continuum_4800_4 = np.mean(alphas_in_range)
continuum_4800_std_4 = np.std(alphas_in_range, ddof=1)
# noise for OIII
alphas_in_range = alphas4[(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
continuum_o3_4 = np.mean(alphas_in_range)
continuum_o3_std_4 = np.std(alphas_in_range, ddof=1)

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
bw1, aw1 = break_4000(alphas1)
bw2, aw2 = break_4000(alphas2)
bw3, aw3 = break_4000(alphas3)
bw4, aw4 = break_4000(alphas4)

# calculate equivalent widths:

# 4863 peak (H-beta)
width_4863_1, err_4863_1 = equiv_width2(4863, alphas1, alpha_std1, continuum_4800_1, continuum_4800_std_1)
width_4863_1 = width_4863_1 - bw1 #correction
width_4960_1, err_4960_1 = equiv_width2(4960, alphas1, alpha_std1, continuum_o3_1, continuum_o3_std_1)
width_5008_1, err_5008_1 = equiv_width2(5008, alphas1, alpha_std1, continuum_o3_1, continuum_o3_std_1)
width_6550_1, err_6550_1 = equiv_width2(6550, alphas1, alpha_std1, continuum_6600_1, continuum_6600_std_1)
# 6565 peak (H-alpha)
width_6565_1, err_6565_1 = equiv_width2(6565, alphas1, alpha_std1, continuum_6600_1, continuum_6600_std_1)
width_6565_1 = width_6565_1 - aw1 #correction
width_6585_1, err_6585_1 = equiv_width2(6585, alphas1, alpha_std1, continuum_6600_1, continuum_6600_std_1)
width_6718_1, err_6718_1 = equiv_width2(6718, alphas1, alpha_std1, continuum_6600_1, continuum_6600_std_1)
width_6733_1, err_6733_1 = equiv_width2(6733, alphas1, alpha_std1, continuum_6600_1, continuum_6600_std_1)
widths1 = [width_4863_1, width_4960_1, width_5008_1, width_6550_1, width_6565_1, width_6585_1, width_6718_1, width_6733_1]
errs1 = [err_4863_1, err_4960_1, err_5008_1, err_6550_1, err_6565_1, err_6585_1, err_6718_1,  err_6733_1]
widths1 = np.reshape(widths1, (len(widths1), 1))
errs1 = np.reshape(errs1, (len(errs1), 1))
ratios1 = np.array([ratios(widths1[i], errs1[i], width_6565_1, err_6565_1, width_4863_1, err_4863_1) for i in range(len(widths1))])
ratios1 = ratios1.reshape((ratios1.shape[0], ratios1.shape[1]))
#N ratio: 6585/6550
N_ratio_1 = width_6585_1 / width_6550_1
N_err_1 = N_ratio_1*(err_6585_1/width_6585_1 + err_6550_1/width_6550_1)
#ISM:
temp_ratio_N_1 = (width_6585_1 + width_6550_1)/width_6565_1
temp_ratio_N_err_1 = temp_ratio_N_1 * (np.sqrt(pow(err_6585_1, 2) + pow(err_6550_1, 2))/(width_6585_1 + width_6550_1) + (err_6565_1/width_6565_1))
temp_ratio_S_1 = (width_6718_1 + width_6733_1)/width_6565_1
temp_ratio_S_err_1 = temp_ratio_S_1 * (np.sqrt(pow(err_6718_1, 2) + pow(err_6733_1, 2))/(width_6718_1 + width_6733_1) + (err_6565_1/width_6565_1))


#4863 peak (H-beta)
width_4863_2, err_4863_2 = equiv_width2(4863, alphas2, alpha_std2, continuum_4800_2, continuum_4800_std_2)
width_4863_2 = width_4863_2 - bw2 #correction
width_4960_2, err_4960_2 = equiv_width2(4960, alphas2, alpha_std2, continuum_o3_2, continuum_o3_std_2)
width_5008_2, err_5008_2 = equiv_width2(5008, alphas2, alpha_std2, continuum_o3_2, continuum_o3_std_2)
width_6550_2, err_6550_2 = equiv_width2(6550, alphas2, alpha_std2, continuum_6600_2, continuum_6600_std_2)
#6565 peak (H-alpha)
width_6565_2, err_6565_2 = equiv_width2(6565, alphas2, alpha_std2, continuum_6600_2, continuum_6600_std_2)
width_6565_2 = width_6565_2 - aw2 #correction
width_6585_2, err_6585_2 = equiv_width2(6585, alphas2, alpha_std2, continuum_6600_2, continuum_6600_std_2)
width_6718_2, err_6718_2 = equiv_width2(6718, alphas2, alpha_std2, continuum_6600_2, continuum_6600_std_2)
width_6733_2, err_6733_2 = equiv_width2(6733, alphas2, alpha_std2, continuum_6600_2, continuum_6600_std_2)
widths2 = [width_4863_2, width_4960_2, width_5008_2, width_6550_2, width_6565_2, width_6585_2, width_6718_2, width_6733_2]
errs2 = [err_4863_2, err_4960_2, err_5008_2, err_6550_2, err_6565_2, err_6585_2, err_6718_2,  err_6733_2]
widths2 = np.reshape(widths2, (len(widths2), 1))
errs2 = np.reshape(errs2, (len(errs2), 1))
ratios2 = np.array([ratios(widths2[i], errs2[i], width_6565_2, err_6565_2, width_4863_2, err_4863_2) for i in range(len(widths2))])
ratios2 = ratios2.reshape((ratios2.shape[0], ratios2.shape[1]))
#N ratio: 6585/6550
N_ratio_2 = width_6585_2 / width_6550_2
N_err_2 = N_ratio_2*(err_6585_2/width_6585_2 + err_6550_2/width_6550_2)
#ISM:
temp_ratio_N_2 = (width_6585_2 + width_6550_2)/width_6565_2
temp_ratio_N_err_2 = temp_ratio_N_2 * (np.sqrt(pow(err_6585_2, 2) + pow(err_6550_2, 2))/(width_6585_2 + width_6550_2) + (err_6565_2/width_6565_2))
temp_ratio_S_2 = (width_6718_2 + width_6733_2)/width_6565_2
temp_ratio_S_err_2 = temp_ratio_S_2 * (np.sqrt(pow(err_6718_2, 2) + pow(err_6733_2, 2))/(width_6718_2 + width_6733_2) + (err_6565_2/width_6565_2))

#4863 peak (H-beta)
width_4863_3, err_4863_3 = equiv_width2(4863, alphas3, alpha_std3, continuum_4800_3, continuum_4800_std_3)
width_4863_3 = width_4863_3 - bw3 #correction
width_4960_3, err_4960_3 = equiv_width2(4960, alphas3, alpha_std3, continuum_o3_3, continuum_o3_std_3)
width_5008_3, err_5008_3 = equiv_width2(5008, alphas3, alpha_std3, continuum_o3_3, continuum_o3_std_3)
width_6550_3, err_6550_3 = equiv_width2(6550, alphas3, alpha_std3, continuum_6600_3, continuum_6600_std_3)
#6565 peak (H-alpha)
width_6565_3, err_6565_3 = equiv_width2(6565, alphas3, alpha_std3, continuum_6600_3, continuum_6600_std_3)
width_6565_3 = width_6565_3 - aw3 #correction
width_6585_3, err_6585_3 = equiv_width2(6585, alphas3, alpha_std3, continuum_6600_3, continuum_6600_std_3)
width_6718_3, err_6718_3 = equiv_width2(6718, alphas3, alpha_std3, continuum_6600_3, continuum_6600_std_3)
width_6733_3, err_6733_3 = equiv_width2(6733, alphas3, alpha_std3, continuum_6600_3, continuum_6600_std_3)
widths3 = [width_4863_3, width_4960_3, width_5008_3, width_6550_3, width_6565_3, width_6585_3, width_6718_3, width_6733_3]
errs3 = [err_4863_3, err_4960_3, err_5008_3, err_6550_3, err_6565_3, err_6585_3, err_6718_3,  err_6733_3]
widths3 = np.reshape(widths3, (len(widths3), 1))
errs3 = np.reshape(errs3, (len(errs3), 1))
ratios3 = np.array([ratios(widths3[i], errs3[i], width_6565_3, err_6565_3, width_4863_3, err_4863_3) for i in range(len(widths3))])
ratios3 = ratios3.reshape((ratios3.shape[0], ratios3.shape[1]))
#N ratio: 6585/6550
N_ratio_3 = width_6585_3 / width_6550_3
N_err_3 = N_ratio_3*(err_6585_3/width_6585_3 + err_6550_3/width_6550_3)
#ISM:
temp_ratio_N_3 = (width_6585_3 + width_6550_3)/width_6565_3
temp_ratio_N_err_3 = temp_ratio_N_3 * (np.sqrt(pow(err_6585_3, 2) + pow(err_6550_3, 2))/(width_6585_3 + width_6550_3) + (err_6565_3/width_6565_3))
temp_ratio_S_3 = (width_6718_3 + width_6733_3)/width_6565_3
temp_ratio_S_err_3 = temp_ratio_S_3 * (np.sqrt(pow(err_6718_3, 2) + pow(err_6733_3, 2))/(width_6718_3 + width_6733_3) + (err_6565_3/width_6565_3))

#4863 peak (H-beta)
width_4863_4, err_4863_4 = equiv_width2(4863, alphas4, alpha_std4, continuum_4800_4, continuum_4800_std_4)
width_4863_4 = width_4863_4 - bw4 #correction
width_4960_4, err_4960_4 = equiv_width2(4960, alphas4, alpha_std4, continuum_o3_4, continuum_o3_std_4)
width_5008_4, err_5008_4 = equiv_width2(5008, alphas4, alpha_std4, continuum_o3_4, continuum_o3_std_4)
width_6550_4, err_6550_4 = equiv_width2(6550, alphas4, alpha_std4, continuum_6600_4, continuum_6600_std_4)
#6565 peak (H-alpha)
width_6565_4, err_6565_4 = equiv_width2(6565, alphas4, alpha_std4, continuum_6600_4, continuum_6600_std_4)
width_6565_4 = width_6565_4 - aw4 #correction
width_6585_4, err_6585_4 = equiv_width2(6585, alphas4, alpha_std4, continuum_6600_4, continuum_6600_std_4)
width_6718_4, err_6718_4 = equiv_width2(6718, alphas4, alpha_std4, continuum_6600_4, continuum_6600_std_4)
width_6733_4, err_6733_4 = equiv_width2(6733, alphas4, alpha_std4, continuum_6600_4, continuum_6600_std_4)
widths4 = [width_4863_4, width_4960_4, width_5008_4, width_6550_4, width_6565_4, width_6585_4, width_6718_4, width_6733_4]
errs4 = [err_4863_4, err_4960_4, err_5008_4, err_6550_4, err_6565_4, err_6585_4, err_6718_4,  err_6733_4]
widths4 = np.reshape(widths4, (len(widths4), 1))
errs4 = np.reshape(errs4, (len(errs4), 1))
ratios4 = np.array([ratios(widths4[i], errs4[i], width_6565_4, err_6565_4, width_4863_4, err_4863_4) for i in range(len(widths4))])
ratios4 = ratios4.reshape((ratios4.shape[0], ratios4.shape[1]))
#N ratio: 6585/6550
N_ratio_4 = width_6585_4 / width_6550_4
N_err_4 = N_ratio_4*(err_6585_4/width_6585_4 + err_6550_4/width_6550_4)
#ISM:
temp_ratio_N_4 = (width_6585_4 + width_6550_4)/width_6565_4
temp_ratio_N_err_4 = temp_ratio_N_4 * (np.sqrt(pow(err_6585_4, 2) + pow(err_6550_4, 2))/(width_6585_4 + width_6550_4) + (err_6565_4/width_6565_4))
temp_ratio_S_4 = (width_6718_4 + width_6733_4)/width_6565_4
temp_ratio_S_err_4 = temp_ratio_S_4 * (np.sqrt(pow(err_6718_4, 2) + pow(err_6733_4, 2))/(width_6718_4 + width_6733_4) + (err_6565_4/width_6565_4))


#tables:
from tabulate import tabulate
titles = ["4863 (H-beta)", "4960 (OIII)", "5008 (OIII)", "6550 (NII)", "6565 (H-alpha)", "6585 (NII)", "6718 (SII)", "6733 (SII)"]
titles = np.reshape(titles, (len(titles), 1))
cell_text1 = np.concatenate((titles, widths1, errs1, ratios1), axis=1)
headers=["Wavelength", "Width", "Err", "H-alpha ratio", "Err", "H-beta ratio", "Err"]
print("\n" + tabulate(cell_text1, headers))
print("N ratio:", N_ratio_1, "+/-", N_err_1)
print("N to alpha: from", temp_ratio_N_1 - temp_ratio_N_err_1, "to", temp_ratio_N_1 + temp_ratio_N_err_1)
print("S to alpha: from", temp_ratio_S_1 - temp_ratio_S_err_1, "to", temp_ratio_S_1 + temp_ratio_S_err_1)

cell_text2 = np.concatenate((titles, widths2, errs2, ratios2), axis=1)
print("\n" + tabulate(cell_text2, headers))
print("N ratio:", N_ratio_2, "+/-", N_err_2)
print("N to alpha: from", temp_ratio_N_2 - temp_ratio_N_err_2, "to", temp_ratio_N_2 + temp_ratio_N_err_2)
print("S to alpha: from", temp_ratio_S_2 - temp_ratio_S_err_2, "to", temp_ratio_S_2 + temp_ratio_S_err_2)

cell_text3 = np.concatenate((titles, widths3, errs3, ratios3), axis=1)
print("\n" + tabulate(cell_text3, headers))
print("N ratio:", N_ratio_3, "+/-", N_err_3)
print("N to alpha: from", temp_ratio_N_3 - temp_ratio_N_err_3, "to", temp_ratio_N_3 + temp_ratio_N_err_3)
print("S to alpha: from", temp_ratio_S_3 - temp_ratio_S_err_3, "to", temp_ratio_S_3 + temp_ratio_S_err_3)

cell_text4 = np.concatenate((titles, widths4, errs4, ratios4), axis=1)
print("\n" + tabulate(cell_text4, headers))
print("N ratio:", N_ratio_4, "+/-", N_err_4)
print("N to alpha: from", temp_ratio_N_4 - temp_ratio_N_err_4, "to", temp_ratio_N_4 + temp_ratio_N_err_4)
print("S to alpha: from", temp_ratio_S_4 - temp_ratio_S_err_4, "to", temp_ratio_S_4 + temp_ratio_S_err_4)


def plot_emissions(alphas1, alphas2, label1, label2):
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
    plot_emissions(alphas1, alphas2, "SFD", r"With $\tau$ Correction")
    plt.show()
    plot_emissions(alphas1, alphas4, "SFD", "With IRIS data")
    plt.show()
    plot_emissions(alphas1, alphas2, "SFD", r"With $\tau$ and IRIS" )
    plt.show()


from scipy import stats

bin_means, bin_edges, binnumber = stats.binned_statistic(wavelength, alphas1, statistic='mean', bins=108)
bin_width = (bin_edges[1] - bin_edges[0])
print("bin width: ", bin_width)
print(bin_edges[0], bin_edges[-1])
plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g')
if plots_on:
    plt.show()



#plot binned alpha vs wavelength (original) 
binned_lambdas = np.linspace(wavelength[0], wavelength[-1], 109) #108 bins
binned_lambdas = np.arange(3900, 9200, 50) #TEST
binned_alphas1 = np.zeros(binned_lambdas.shape)
binned_std1 = np.zeros(binned_lambdas.shape)
for i, lmda in enumerate(binned_lambdas):
    indices = np.where((wavelength > lmda) & (wavelength < lmda+50))[0]
    relevant_alphas = alphas1[indices]
    relevant_stds = alpha_std1[indices]
    #weighted average:
    variance = np.power(relevant_stds, 2)
    numerator = np.sum(np.divide(relevant_alphas, variance))
    denominator = np.sum(np.divide(1, variance))
    avg1 = numerator / denominator
    avg2 = 1 / denominator
    binned_alphas1[i] = avg1
    binned_std1[i] = np.sqrt(avg2)

#plot binned alpha vs wavelength (tao)
binned_alphas2 = np.zeros(binned_lambdas.shape)
binned_std2 = np.zeros(binned_lambdas.shape)
for i, lmda in enumerate(binned_lambdas):
    indices = np.where((wavelength > lmda-25) & (wavelength < lmda+25))[0]
    relevant_alphas = alphas2[indices]
    relevant_stds = alpha_std2[indices]
    #weighted average:
    variance = np.power(relevant_stds, 2)
    numerator = np.sum(np.divide(relevant_alphas, variance))
    denominator = np.sum(np.divide(1, variance))
    avg1 = numerator / denominator
    avg2 = 1 / denominator
    binned_alphas2[i] = avg1
    binned_std2[i] = np.sqrt(avg2)

#plot binned alpha vs wavelength (tao AND iris)
binned_alphas3 = np.zeros(binned_lambdas.shape)
binned_std3 = np.zeros(binned_lambdas.shape)
for i, lmda in enumerate(binned_lambdas):
    indices = np.where((wavelength > lmda-25) & (wavelength < lmda+25))[0]
    relevant_alphas = alphas3[indices]
    relevant_stds = alpha_std3[indices]
    #weighted average:
    variance = np.power(relevant_stds, 2)
    numerator = np.sum(np.divide(relevant_alphas, variance))
    denominator = np.sum(np.divide(1, variance))
    avg1 = numerator / denominator
    avg2 = 1 / denominator
    binned_alphas3[i] = avg1
    binned_std3[i] = np.sqrt(avg2)

#plot binned alpha vs wavelength (iris)
binned_alphas4 = np.zeros(binned_lambdas.shape)
binned_std4 = np.zeros(binned_lambdas.shape)
for i, lmda in enumerate(binned_lambdas):
    indices = np.where((wavelength > lmda-25) & (wavelength < lmda+25))[0]
    relevant_alphas = alphas4[indices]
    relevant_stds = alpha_std4[indices]
    #weighted average:
    variance = np.power(relevant_stds, 2)
    numerator = np.sum(np.divide(relevant_alphas, variance))
    denominator = np.sum(np.divide(1, variance))
    avg1 = numerator / denominator
    avg2 = 1 / denominator
    binned_alphas4[i] = avg1
    binned_std4[i] = np.sqrt(avg2)

plt.figure(figsize=(6, 14))
    
plt.subplot(3, 1, 1)    
#compare original to tao
plt.plot(binned_lambdas, binned_alphas1, c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_std1, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas2, c='r', drawstyle='steps', label=r'With $\tau$ Correction')
plt.plot(binned_lambdas, binned_std2, c='r', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(3850, 9200)
plt.ylim(0, 0.4)
plt.legend(frameon=False)

plt.subplot(3, 1, 2)
#compare original to iris
plt.plot(binned_lambdas, binned_alphas1, c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_std1, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas4, c='r', drawstyle='steps', label='With IRIS Data')
plt.plot(binned_lambdas, binned_std4, c='r', drawstyle='steps')
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"$\alpha_\lambda$")
plt.xlim(3850, 9200)
plt.ylim(0, 0.4)
plt.legend(frameon=False)

plt.subplot(3, 1, 3)
#compare original to tao AND iris
plt.plot(binned_lambdas, binned_alphas1, c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_std1, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas3, c='r', drawstyle='steps', label=r'With IRIS and $\tau$')
plt.plot(binned_lambdas, binned_std3, c='r', drawstyle='steps')
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
plt.plot(binned_lambdas, binned_alphas1, c='k', drawstyle='steps', label='SFD')
plt.plot(binned_lambdas, binned_std1, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas3, c='r', drawstyle='steps', label=r'With IRIS and $\tau$')
plt.plot(binned_lambdas, binned_std3, c='r', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas2, c='b', drawstyle='steps', label='With tao corr.')
plt.plot(binned_lambdas, binned_std2, c='b', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas4, c='g', drawstyle='steps', label=r'with IRIS')
plt.plot(binned_lambdas, binned_std4, c='g', drawstyle='steps')
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


