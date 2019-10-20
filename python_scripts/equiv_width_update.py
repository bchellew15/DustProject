# calculations to find equivalent widths, line ratios, temperatures

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

    alphas[0] = np.load('../alphas_and_stds/alphas_sdss_83019.npy') / 1.38 #temp
    alpha_stds[0] = np.load('../alphas_and_stds/alphas_sdss_stds_83019.npy') / 1.38 #temp

    alphas[1] = [0.18972603, 0.18287671, 0.1630137, 0.14657534, 0.09863014, 0.1609589, 0.18150685, 0.1890411, 0.16986301, 0.17328767, \
                  0.20547945, 0.20205479, 0.20205479, 0.23835616, 0.26438356, 0.24383562, 0.20136986, 0.16986301, 0.1390411, 0.13630137, 0.17876712, \
                  0.1739726, 0.28150685, 0.53835616, 0.61438356, 0.49452055, 0.31780822, 0.17465753, 0.19041096, 0.16643836, 0.16849315, \
                  0.16164384, 0.19383562, 0.19520548, 0.1760274, 0.18356164, 0.24109589, 0.33630137, 0.44315068, 0.4, 0.25, 0.20616438, \
                  0.1760274, 0.16780822, 0.17123288, 0.16712329, 0.17671233, 0.17260274, 0.17328767, 0.20205479, 0.17671233, 0.17191781, \
                  0.17808219, 0.20479452, 0.19383562, 0.1739726, 0.19383562, 0.19863014, 0.18630137, 0.19109589, 0.19246575, 0.18013699, \
                  0.17945205, 0.17671233, 0.16917808, 0.17260274, 0.16369863, 0.17671233, 0.18013699, 0.19726027, 0.16643836, 0.18561644, \
                  0.1739726, 0.17123288, 0.16712329, 0.18219178, 0.18630137, 0.17945205, 0.18287671, 0.17808219, 0.1760274, 0.19383562, \
                  0.16506849, 0.14931507, 0.16506849, 0.16575342, 0.17671233, 0.18630137, 0.18356164, 0.17534247, 0.18561644, 0.1869863, \
                  0.1609589, 0.17808219, 0.19109589, 0.18561644, 0.13424658, 0.17191781, 0.17945205, 0.18561644, 0.17534247, 0.1869863, \
                  0.17945205, 0.17739726, 0.17808219, 0.18630137, 0.18356164, 0.19726027, 0.20342466, 0.17739726, 0.15890411, 0.17260274, \
                  0.17945205, 0.17534247, 0.18082192, 0.16849315, 0.18219178, 0.20205479, 0.19794521, 0.16712329, 0.16849315, 0.18493151, \
                  0.17808219, 0.19520548, 0.17876712, 0.18835616, 0.19246575, 0.26712329, 0.37671233, 0.40410959, 0.3239726, 0.21849315, \
                  0.17465753, 0.19041096, 0.19383562, 0.18287671, 0.20273973, 0.26643836, 0.33972603, 0.34383562, 0.26712329, 0.19863014, \
                  0.17945205, 0.1869863, 0.19657534, 0.17260274, 0.18493151, 0.16369863, 0.17876712, 0.20273973, 0.20890411, 0.1869863, \
                  0.18013699, 0.16849315, 0.16643836, 0.18493151, 0.19041096, 0.17260274, 0.17739726, 0.17328767, 0.15958904, 0.18287671, \
                  0.19657534, 0.18561644, 0.19178082] #from Brandt paper
    alphas[1] = np.pad(alphas[1], pad_width=((2449, 1386)), mode='constant', constant_values=(1))

num_arrays = len(alphas)
num_regions = 8 #I think this is number of peaks to look at

# integrate using 3 wavelength elements on either side
def equiv_width2(peak_l, alphas, alpha_stds, cont, stdev):

    #find indices of 2 closest wavelengths
    #peak_idx1 = np.argpartition(np.abs(wavelength-peak_l), 0)[0]
    #peak_idx2 = np.argpartition(np.abs(wavelength-peak_l), 1)[1]
    #left_peak_idx = min(peak_idx1, peak_idx2)
    left_peak_idx = np.argmin(np.abs(wavelength-peak_l)) #temp. this is just closest to peak

    delta_lambda = wavelength[left_peak_idx+1] - wavelength[left_peak_idx]
    
    bars = (alphas-cont)/cont * delta_lambda
    
    l_off = -2 # -2: this is shifter half-unit to the right
    r_off = 4 #l_off + 6

    
    #peak fitting
    
    def gaussian_func(x, a1, a2, a3, a4):
        return a1 + a2* np.exp(-np.power(x-a3, 2)/(2*a4**2))

    def width_helper(x, a1, a2, a3, a4):
        return (gaussian_func(x, a1, a2, a3, a4) - a1) / a1
    
    #better peak fitting: (move this outside of this function)
    def linear_fit(rel_alphas, rel_lambdas, rel_sigmas, a3, a4):

        rel_ivars = 1 / np.power(rel_sigmas, 2)
        y = np.exp(-np.power(rel_lambdas-a3, 2)/(2*a4**2))
        A = [[np.sum(rel_ivars), np.sum(np.multiply(y, rel_ivars))], \
             [np.sum(np.multiply(y, rel_ivars)), np.sum(np.multiply(np.power(y, 2), rel_ivars))]]
        b = [np.sum(np.multiply(rel_alphas, rel_ivars)), np.sum(np.multiply(np.multiply(rel_alphas, y), rel_ivars))]
        
        a1, a2 = lstsq(A, b)[0]
        return a1, a2, y

    def nonlinear_helper(alphas, sigmas):
        def nonlinear_fit(rel_lambdas, a3, a4):
            rel_sigmas = sigmas
            rel_alphas = alphas
            a1, a2, y = linear_fit(rel_alphas, rel_lambdas, rel_sigmas, a3, a4)
            return a1 + a2*y
        return nonlinear_fit

    #guesses?
    #try it just for a couple lines
    #6585 nII peak (tall one): 6574.6 to 6708
    n2_left_idx = np.argmin(np.abs(wavelength-6574.6))
    n2_right_idx = np.argmin(np.abs(wavelength-6708))
    rel_alphas = alphas[n2_left_idx:n2_right_idx+1]
    rel_sigmas = alpha_stds[n2_left_idx:n2_right_idx+1]
    rel_lambdas = wavelength[n2_left_idx:n2_right_idx+1]
    popt,pcov = curve_fit(nonlinear_helper(rel_alphas, rel_sigmas), rel_lambdas, rel_alphas, p0=[6580, 10], sigma=rel_sigmas)

    #recover a1 and a2:
    a1, a2, y = linear_fit(rel_alphas, rel_lambdas, rel_sigmas, popt[0], popt[1])

    print("params:", a1, a2, popt[0], popt[1])

    x_range = np.arange(6550, 6750, .01)
    y = np.exp(-np.power(x_range-popt[0], 2)/(2*popt[1]**2))
    alpha_pred = a1+a2*y
    plt.plot(x_range, alpha_pred, '.')
    plt.plot(rel_lambdas, rel_alphas, 'r.')
    plt.show()

    print("left:")
    print(popt[0] - 3*popt[1])
    print("right:")
    print(popt[0] + 3*popt[1])
    
    #now integrate (over about 3 sigma I think)
    width = quad(width_helper, popt[0]-10*popt[1], popt[0] + 10*popt[1], args=(a1, a2, popt[0], popt[1]))
    print("n2 width:")
    print(width)

    exit(0)
    
    width = np.sum(bars[left_peak_idx+l_off:left_peak_idx+r_off])

    #print("wavelengths and bars:")
    #print(wavelength[left_peak_idx+l_off-1:left_peak_idx+r_off+1]) #temp
    #print(bars[left_peak_idx+l_off-1:left_peak_idx+r_off+1]) #temp

    # each bar has error bar*sigma*(1/peak + 1/continuum)*delta_lambda
    # then add in quadrature
    # (because we can simplify to alpha/continuum - 1)
    #bar_errs = delta_lambda*np.array([bars[i]*(alpha_stds[i]/alphas[i] + stdev/cont) for i in range(left_peak_idx+l_off, left_peak_idx+r_off)])
    
    bar_errs = delta_lambda*(alphas/cont)*np.sqrt(np.power(np.divide(alpha_stds, alphas), 2) + (stdev/cont)**2)
    err = np.sqrt(np.sum(np.power(bar_errs[left_peak_idx+l_off:left_peak_idx+r_off], 2)))

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
    stds_in_range = alpha_stds[i][(wavelength>6600) * (wavelength < 6700)]
    continuum_6600s[i] = np.average(alphas_in_range, weights = 1/np.power(stds_in_range, 2))
    continuum_6600_stds[i] = np.sqrt(np.sum(np.power(stds_in_range, 2))) / len(alphas_in_range)   
    #np.std(alphas_in_range, ddof=1)
    # noise for H-beta
    alphas_in_range = alphas[i][(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
    stds_in_range = alpha_stds[i][(wavelength>4829)*(wavelength<4858) + (wavelength>4864)*(wavelength<4893)]
    continuum_4800s[i] = np.average(alphas_in_range, weights = 1/np.power(stds_in_range, 2))
    continuum_4800_stds[i] = np.sqrt(np.sum(np.power(stds_in_range, 2))) / len(alphas_in_range) 
    #np.std(alphas_in_range, ddof=1)
    # noise for OIII
    alphas_in_range = alphas[i][(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
    stds_in_range = alpha_stds[i][(wavelength>4900)*(wavelength<4958) + (wavelength>4963)*(wavelength<5000)]
    continuum_o3s[i] = np.average(alphas_in_range, weights = 1/np.power(stds_in_range, 2))
    continuum_o3_stds[i] = np.sqrt(np.sum(np.power(stds_in_range, 2))) / len(alphas_in_range)  
    #np.std(alphas_in_range, ddof=1)
    

    #Updated O3 range:
    alphas_in_range = alphas[i][(wavelength>4914)*(wavelength<4957) + (wavelength>4963)*(wavelength<5000)] #temp
    stds_in_range = alpha_stds[i][(wavelength>4914)*(wavelength<4957) + (wavelength>4963)*(wavelength<5000)]
    #in original o3 range it should be 4957, I think
    continuum_o3s[i] = np.average(alphas_in_range, weights = 1/np.power(stds_in_range, 2)) #temp
    continuum_o3_stds[i] = np.sqrt(np.sum(np.power(stds_in_range, 2))) / len(alphas_in_range)
    
    print("n2 continuum;")
    print(continuum_6600s[i])
    

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
bws = np.zeros(num_arrays)
aws = np.zeros(num_arrays)
for i in range(num_arrays):
    aws[i], bws[i] = break_4000(alphas[i])

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

l_4863 = 4863 #4862 is highest values
l_4960 = 4960
l_5008 = 5008
l_6550 = 6550
l_6565 = 6565
l_6585 = 6585
l_6718 = 6718
l_6733 = 6733


for i in range(num_arrays):
    # 4863 peak (H-beta)
    width_4863s[i], err_4863s[i] = equiv_width2(l_4863, alphas[i], alpha_stds[i], continuum_4800s[i], continuum_4800_stds[i])
    width_4863s[i] = width_4863s[i] - bws[i] #correction
    width_4960s[i], err_4960s[i] = equiv_width2(l_4960, alphas[i], alpha_stds[i], continuum_o3s[i], continuum_o3_stds[i])
    width_5008s[i], err_5008s[i] = equiv_width2(l_5008, alphas[i], alpha_stds[i], continuum_o3s[i], continuum_o3_stds[i])
    width_6550s[i], err_6550s[i] = equiv_width2(l_6550, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
    # 6565 peak (H-alpha)
    width_6565s[i], err_6565s[i] = equiv_width2(l_6565, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
    width_6565s[i] = width_6565s[i] - aws[i] #correction
    width_6585s[i], err_6585s[i] = equiv_width2(l_6585, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
    width_6718s[i], err_6718s[i] = equiv_width2(l_6718, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
    width_6733s[i], err_6733s[i] = equiv_width2(l_6733, alphas[i], alpha_stds[i], continuum_6600s[i], continuum_6600_stds[i])
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

    if i != 0:
        break #temp

    cell_text = np.concatenate((titles, widths[i], errs[i], ratios[i]), axis=1)
    print("\n" + tabulate(cell_text, headers))
    print("N ratio:", N_ratios[i], "+/-", N_errs[i])
    print("N to alpha: from", temp_ratio_Ns[i] - temp_ratio_N_errs[i], "to", temp_ratio_Ns[i] + temp_ratio_N_errs[i])
    print("S to alpha: from", temp_ratio_Ss[i] - temp_ratio_S_errs[i], "to", temp_ratio_Ss[i] + temp_ratio_S_errs[i])

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

