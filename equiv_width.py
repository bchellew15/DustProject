#load in alphas and make a scatterplot

#std dev: remember to divide by N-1...
#sanity check on bounds: maybe plot vert lines in red

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy import interpolate
from scipy import integrate

import warnings #for integration warnings... maybe investigate
warnings.simplefilter("ignore")

#load wavelengths
hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
wavelength = np.array(hdulist[1].data)  # Angstroms

#load in npy files
alphas1 = np.load('alphas_1d.npy')
alpha_std1 = np.load('alphas_1d_stds.npy')
alphas2 = np.load('alphas_2d.npy')
alpha_std2 = np.load('alphas_2d_stds.npy')
alphas3 = np.load('alphas_iris.npy')
alpha_std3 = np.load('alphas_iris_stds.npy')

#integrate using 6-7 wavelength elements on either side
def equiv_width2(peak_l, alphas, cont, stdev):
    peak_idx = np.argmin(np.abs(wavelength-peak_l))
    width = np.sum((alphas[peak_idx-3:peak_idx+4]-cont)/cont) #this is width 7, not 6

    #each bar has error sigma*(1/peak + 1/continuum)
    #then add in quadrature
    bar_errs = np.array([stdev*(1/alphas[i] + 1/cont) for i in range(peak_idx-3, peak_idx+4)])
    err = np.sqrt(np.sum(np.power(bar_errs, 2)))

    return width, err

#return H-alpha and H-beta ratios, with errors
def ratios(width, err, halpha, a_err, hbeta, b_err):
    a_ratio = width/halpha
    b_ratio = width/hbeta
    a_ratio_err = a_ratio*(a_err/halpha + err/width)
    b_ratio_err = b_ratio*(b_err/hbeta + err/width)
    return a_ratio, a_ratio_err, b_ratio, b_ratio_err

#noise for alpha1:
alphas_in_range = alphas1[(wavelength>6600) * (wavelength < 6700)]
continuum_6600_1 = np.mean(alphas_in_range) #should be weighted mean?
continuum_6600_std_1 = np.std(alphas_in_range) #/N-1?
#noise for H-beta
alphas_in_range = alphas1[(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
continuum_4800_1 = np.mean(alphas_in_range)
continuum_4800_std_1 = np.std(alphas_in_range)
#noise for OIII
alphas_in_range = alphas1[(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
continuum_o3_1 = np.mean(alphas_in_range)
continuum_o3_std_1 = np.std(alphas_in_range)

#noise for alpha2:
alphas_in_range = alphas2[(wavelength>6600) * (wavelength < 6700)]
continuum_6600_2 = np.mean(alphas_in_range) #should be weighted mean?
continuum_6600_std_2 = np.std(alphas_in_range) #/N-1?
#noise for H-beta
alphas_in_range = alphas2[(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
continuum_4800_2 = np.mean(alphas_in_range)
continuum_4800_std_2 = np.std(alphas_in_range)
#noise for OIII
alphas_in_range = alphas2[(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
continuum_o3_2 = np.mean(alphas_in_range)
continuum_o3_std_2 = np.std(alphas_in_range)

#noise for alpha3:
alphas_in_range = alphas3[(wavelength>6600) * (wavelength < 6700)]
continuum_6600_3 = np.mean(alphas_in_range) #should be weighted mean?
continuum_6600_std_3 = np.std(alphas_in_range) #/N-1?
#noise for H-beta
alphas_in_range = alphas3[(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
continuum_4800_3 = np.mean(alphas_in_range)
continuum_4800_std_3 = np.std(alphas_in_range)
#noise for OIII
alphas_in_range = alphas3[(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
continuum_o3_3 = np.mean(alphas_in_range)
continuum_o3_std_3 = np.std(alphas_in_range)

#calculate 4000 A break
def break_4000(alphas):
    alpha_func = interpolate.interp1d(wavelength, alphas)
    left_break = integrate.quad(alpha_func, 3850, 4000)[0]
    right_break = integrate.quad(alpha_func, 4000, 4150)[0]
    delta = left_break / right_break
    beta_width = -2.2*delta + .17
    alpha_width = -1.5*delta - .19
    return alpha_width, beta_width

#calculate 4000 A break: [error...]
bw1, aw1 = break_4000(alphas1)
bw2, aw2 = break_4000(alphas2)
bw3, aw3 = break_4000(alphas3)

#calculate equivalent widths:
#5008 looks like double peak
#error prop from H corrections?

#4863 peak (H-beta)
width_4863_1, err_4863_1 = equiv_width2(4863, alphas1, continuum_4800_1, continuum_4800_std_1)
width_4863_1 = width_4863_1 - bw1 #correction
width_4960_1, err_4960_1 = equiv_width2(4960, alphas1, continuum_o3_1, continuum_o3_std_1)
width_5008_1, err_5008_1 = equiv_width2(5008, alphas1, continuum_o3_1, continuum_o3_std_1)
width_6550_1, err_6550_1 = equiv_width2(6550, alphas1, continuum_6600_1, continuum_6600_std_1)
#6565 peak (H-alpha)
width_6565_1, err_6565_1 = equiv_width2(6565, alphas1, continuum_6600_1, continuum_6600_std_1)
width_6565_1 = width_6565_1 - aw1 #correction
width_6585_1, err_6585_1 = equiv_width2(6585, alphas1, continuum_6600_1, continuum_6600_std_1)
width_6718_1, err_6718_1 = equiv_width2(6718, alphas1, continuum_6600_1, continuum_6600_std_1)
width_6733_1, err_6733_1 = equiv_width2(6733, alphas1, continuum_6600_1, continuum_6600_std_1)
widths1 = [width_4863_1, width_4960_1, width_5008_1, width_6550_1, width_6565_1, width_6585_1, width_6718_1, width_6733_1]
errs1 = [err_4863_1, err_4960_1, err_5008_1, err_6550_1, err_6565_1, err_6585_1, err_6718_1,  err_6733_1]
widths1 = np.reshape(widths1, (len(widths1), 1))
errs1 = np.reshape(errs1, (len(errs1), 1))
ratios1 = np.array([ratios(widths1[i], errs1[i], width_6565_1, err_6565_1, width_4863_1, err_4863_1) for i in range(len(widths1))])
ratios1 = ratios1.reshape((ratios1.shape[0], ratios1.shape[1]))

#4863 peak (H-beta)
width_4863_2, err_4863_2 = equiv_width2(4863, alphas2, continuum_4800_2, continuum_4800_std_2)
width_4863_2 = width_4863_2 - bw2 #correction
width_4960_2, err_4960_2 = equiv_width2(4960, alphas2, continuum_o3_2, continuum_o3_std_2)
width_5008_2, err_5008_2 = equiv_width2(5008, alphas2, continuum_o3_2, continuum_o3_std_2)
width_6550_2, err_6550_2 = equiv_width2(6550, alphas2, continuum_6600_2, continuum_6600_std_2)
#6565 peak (H-alpha)
width_6565_2, err_6565_2 = equiv_width2(6565, alphas2, continuum_6600_2, continuum_6600_std_2)
width_6565_2 = width_6565_2 - aw2 #correction
width_6585_2, err_6585_2 = equiv_width2(6585, alphas2, continuum_6600_2, continuum_6600_std_2)
width_6718_2, err_6718_2 = equiv_width2(6718, alphas2, continuum_6600_2, continuum_6600_std_2)
width_6733_2, err_6733_2 = equiv_width2(6733, alphas2, continuum_6600_2, continuum_6600_std_2)
widths2 = [width_4863_2, width_4960_2, width_5008_2, width_6550_2, width_6565_2, width_6585_2, width_6718_2, width_6733_2]
errs2 = [err_4863_2, err_4960_2, err_5008_2, err_6550_2, err_6565_2, err_6585_2, err_6718_2,  err_6733_2]
widths2 = np.reshape(widths2, (len(widths2), 1))
errs2 = np.reshape(errs2, (len(errs2), 1))
ratios2 = np.array([ratios(widths2[i], errs2[i], width_6565_2, err_6565_2, width_4863_2, err_4863_2) for i in range(len(widths2))])
ratios2 = ratios2.reshape((ratios2.shape[0], ratios2.shape[1]))

#4863 peak (H-beta)
width_4863_3, err_4863_3 = equiv_width2(4863, alphas3, continuum_4800_3, continuum_4800_std_3)
width_4863_3 = width_4863_3 - bw3 #correction
width_4960_3, err_4960_3 = equiv_width2(4960, alphas3, continuum_o3_3, continuum_o3_std_3)
width_5008_3, err_5008_3 = equiv_width2(5008, alphas3, continuum_o3_3, continuum_o3_std_3)
width_6550_3, err_6550_3 = equiv_width2(6550, alphas3, continuum_6600_3, continuum_6600_std_3)
#6565 peak (H-alpha)
width_6565_3, err_6565_3 = equiv_width2(6565, alphas3, continuum_6600_3, continuum_6600_std_3)
width_6565_3 = width_6565_3 - aw3 #correction
width_6585_3, err_6585_3 = equiv_width2(6585, alphas3, continuum_6600_3, continuum_6600_std_3)
width_6718_3, err_6718_3 = equiv_width2(6718, alphas3, continuum_6600_3, continuum_6600_std_3)
width_6733_3, err_6733_3 = equiv_width2(6733, alphas3, continuum_6600_3, continuum_6600_std_3)
widths3 = [width_4863_3, width_4960_3, width_5008_3, width_6550_3, width_6565_3, width_6585_3, width_6718_3, width_6733_3]
errs3 = [err_4863_3, err_4960_3, err_5008_3, err_6550_3, err_6565_3, err_6585_3, err_6718_3,  err_6733_3]
widths3 = np.reshape(widths3, (len(widths3), 1))
errs3 = np.reshape(errs3, (len(errs3), 1))
ratios3 = np.array([ratios(widths3[i], errs3[i], width_6565_3, err_6565_3, width_4863_3, err_4863_3) for i in range(len(widths3))])
ratios3 = ratios3.reshape((ratios3.shape[0], ratios3.shape[1]))

#tables:
from tabulate import tabulate
titles = ["4863 (H-beta)", "4960 (OIII)", "5008 (OIII)", "6550 (NII)", "6565 (H-alpha)", "6585 (NII)", "6718 (SII)", "6733 (SII)"]
titles = np.reshape(titles, (len(titles), 1))
cell_text1 = np.concatenate((titles, widths1, errs1, ratios1), axis=1)
headers=["Wavelength", "Width", "Err", "H-alpha ratio", "Err", "H-beta ratio", "Err"]
print("\n" + tabulate(cell_text1, headers))

#N ratio: 6585/6550
N_ratio_1 = width_6585_1 / width_6550_1
N_err_1 = N_ratio_1*(err_6585_1/width_6585_1 + err_6550_1/width_6550_1)
print("N ratio:", N_ratio_1, "+/-", N_err_1)

cell_text2 = np.concatenate((titles, widths2, errs2, ratios2), axis=1)
print("\n" + tabulate(cell_text2, headers))

#N ratio: 6585/6550
N_ratio_2 = width_6585_2 / width_6550_2
N_err_2 = N_ratio_2*(err_6585_2/width_6585_2 + err_6550_2/width_6550_2)
print("N ratio:", N_ratio_2, "+/-", N_err_2)

cell_text3 = np.concatenate((titles, widths3, errs3, ratios3), axis=1)
print("\n" + tabulate(cell_text3, headers))

#N ratio: 6585/6550
N_ratio_3 = width_6585_3 / width_6550_3
N_err_3 = N_ratio_3*(err_6585_3/width_6585_3 + err_6550_3/width_6550_3)
print("N ratio:", N_ratio_3, "+/-", N_err_3)


#plot 4830 - 5040 (original vs tao)
plt.plot(wavelength, alphas1, c='k', drawstyle='steps')
plt.plot(wavelength, alpha_std1, c='k', drawstyle='steps')
plt.plot(wavelength, alphas2, c='r', drawstyle='steps')
plt.plot(wavelength, alpha_std2, c='r', drawstyle='steps')

plt.xlim(4830, 5040)
plt.ylim(0, 0.6)
xcoords = [4863, 4960, 5008]
for xc in xcoords:
    plt.axvline(x=xc, color='k', linewidth=1, linestyle='--')
plt.title("4830 - 5040 (original vs tao)")
plt.show()

#plot 6530 - 6770 (original vs tao)
plt.plot(wavelength, alphas1, c='k', drawstyle='steps')
plt.plot(wavelength, alpha_std1, c='k', drawstyle='steps')
plt.plot(wavelength, alphas2, c='r', drawstyle='steps')
plt.plot(wavelength, alpha_std2, c='r', drawstyle='steps')

plt.xlim(6530, 6770)
plt.ylim(0, 1.1)
xcoords = [6550, 6565, 6585, 6718, 6733]
for xc in xcoords:
    plt.axvline(x=xc, color='k', linewidth=1, linestyle='--')
plt.title("6530 - 6770 (original vs tao)")
#vert lines:
peak_idx = np.argmin(np.abs(wavelength-6585))
plt.axvline(x=wavelength[peak_idx-3], color='r', linewidth=1, linestyle='--')
plt.axvline(x=wavelength[peak_idx+3], color='r', linewidth=1, linestyle='--')
plt.show()

#plot 4830 - 5040 (iris vs tao)
plt.plot(wavelength, alphas2, c='k', drawstyle='steps')
plt.plot(wavelength, alpha_std2, c='k', drawstyle='steps')
plt.plot(wavelength, alphas3, c='r', drawstyle='steps')
plt.plot(wavelength, alpha_std3, c='r', drawstyle='steps')

plt.xlim(4830, 5040)
plt.ylim(0, 0.6)
xcoords = [4863, 4960, 5008]
for xc in xcoords:
    plt.axvline(x=xc, color='k', linewidth=1, linestyle='--')
plt.title("4830 - 5040 (iris vs tao)")
plt.show()

#plot 6530 - 6770 (iris vs tao)
plt.plot(wavelength, alphas2, c='k', drawstyle='steps')
plt.plot(wavelength, alpha_std2, c='k', drawstyle='steps')
plt.plot(wavelength, alphas3, c='r', drawstyle='steps')
plt.plot(wavelength, alpha_std3, c='r', drawstyle='steps')

plt.xlim(6530, 6770)
plt.ylim(0, 1.1)
xcoords = [6550, 6565, 6585, 6718, 6733]
for xc in xcoords:
    plt.axvline(x=xc, color='k', linewidth=1, linestyle='--')
plt.title("6530 - 6770 (iris vs tao)")
#vert lines:
peak_idx = np.argmin(np.abs(wavelength-6585))
plt.axvline(x=wavelength[peak_idx-3], color='r', linewidth=1, linestyle='--')
plt.axvline(x=wavelength[peak_idx+3], color='r', linewidth=1, linestyle='--')
plt.show()



#plot binned alpha vs wavelength (original) 
binned_lambdas = np.arange(3900, 9200, 50)
binned_alphas1 = np.zeros(binned_lambdas.shape)
binned_std1 = np.zeros(binned_lambdas.shape)
for i, lmda in enumerate(binned_lambdas):
    indices = np.where((wavelength > lmda-25) & (wavelength < lmda+25))[0]
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

#plot binned alpha vs wavelength (optical correction)
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

#plot binned alpha vs wavelength (optical correction AND iris)
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

#compare original to tao
plt.plot(binned_lambdas, binned_alphas1, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_std1, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas2, c='r', drawstyle='steps')
plt.plot(binned_lambdas, binned_std2, c='r', drawstyle='steps')
plt.title("1 vs 2")
plt.show()

#compare original to IRIS
plt.plot(binned_lambdas, binned_alphas1, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_std1, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas3, c='r', drawstyle='steps')
plt.plot(binned_lambdas, binned_std3, c='r', drawstyle='steps')
plt.title("1 vs 3")
plt.show()

#compare original to IRIS
plt.plot(binned_lambdas, binned_alphas2, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_std2, c='k', drawstyle='steps')
plt.plot(binned_lambdas, binned_alphas3, c='r', drawstyle='steps')
plt.plot(binned_lambdas, binned_std3, c='r', drawstyle='steps')
plt.title("2 vs 3")
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


