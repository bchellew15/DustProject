#code to deal with the star formation files
#calculate Lick indices
#find best-fit spectra

import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import interp1d
from scipy.optimize import minimize

paths = glob.glob('/Users/blakechellew/Documents/DustProject/BrandtFiles/bc03/*.spec')

#print out indices and filenames:
#for i in range(len(paths)):
#    print(i, paths[i])

#alphas and stds (boss):
alphas = [np.load('../alphas_and_stds/alphas_boss_102019.npy'), \
               np.load('../alphas_and_stds/alphas_boss_2d_102019.npy')]
alpha_stds = [np.load('../alphas_and_stds/alpha_stds_boss_102019.npy'), \
              np.load('../alphas_and_stds/alpha_stds_boss_2d_102019.npy')]
wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')

#mask emission lines
emission_line_mask = np.zeros(len(wavelength), dtype=int)
emission_lines = [3727, 4863, 4960, 5008, 5877, 6550, 6565, 6585, 6718, 6733]
for line in emission_lines:
    peak_idx = np.argmin(np.abs(wavelength-line))
    emission_line_mask[peak_idx-3:peak_idx+4] = 1
for i in range(len(alphas)):
    alphas[i] = alphas[i][np.logical_not(emission_line_mask)]
    alpha_stds[i] = alpha_stds[i][np.logical_not(emission_line_mask)]
wavelength = wavelength[np.logical_not(emission_line_mask)]


#find best-fit linear combination [plus polynomial] using max likelihood estimator
#and plot the best-fit spectrum against actual spectrum
'''
#for adding a polynomial
poly_degree = 0 #1 is linear, etc.
poly_order = poly_degree + 1

#load the spectra and interpolate to same wavelength range:
model_spectra = np.zeros((len(paths)+poly_order, len(wavelength)))
for i, p in enumerate(paths):
    a = np.loadtxt(p)
    wav = a[:,0]
    values = a[:,1]
    f = interp1d(wav, values, kind='cubic')
    model_spectra[i] = np.array([f(w) for w in wavelength])
#scale them based on total flux from 4000 to 9000
for i in range(len(paths)):
    model_spectra[i] = model_spectra[i] / np.sum(model_spectra[i]) * np.sum(model_spectra[18])
for i in range(poly_order):
    model_spectra[-i-1] = np.power(wavelength, i)
    
num_spectra = len(paths)+poly_order
A = np.zeros((num_spectra, num_spectra))
for i in range(num_spectra):
    for j in range(num_spectra):
        A[i, j] = np.nansum(model_spectra[i]*model_spectra[j]/np.power(alpha_stds[1], 2))
b = np.zeros(num_spectra)
for i in range(num_spectra):
    b[i] = np.nansum(alphas[1]*model_spectra[i]/np.power(alpha_stds[1], 2))

print(A)
print(b)
    
coeffs = np.linalg.lstsq(A, b)[0]
print("coefficients")
print(coeffs)

#plot the best fit spectrum:
#best_fit_model = 
weighted_model_spectra = np.copy(model_spectra)
for i in range(len(coeffs)):
    weighted_model_spectra[i] *= coeffs[i]
best_fit_model = np.sum(weighted_model_spectra, axis=0)

#find the chai squared:
print("chai squared")
print(np.nansum(np.power(best_fit_model - alphas[1], 2)))
    
plt.plot(wavelength, best_fit_model, 'r', drawstyle='steps')
plt.plot(wavelength, alphas[1], 'k', drawstyle='steps')
plt.xlim(3500, 10000)
plt.ylim(0, 1)
plt.show()
'''


#plot the model spectra 

for p in paths[:10]:
    if 'ssp_5Myr_z008' in p:
        a = np.loadtxt(p)
        
        wav = a[:,0]
        values = a[:,1]
        
        plt.plot(wav, values, drawstyle='steps')
        plt.xlim(3500, 10000)
        plt.show()
#exit(0)


#some important lick index tuples:
h_beta_bounds = (4847.875, 4876.625, 4827.875, 4847.875, 4876.625, 4891.625, 0)
mg1_bounds = (5069.125, 5134.125, 4895.125, 4957.625, 5301.125, 5366.125, 1)
mg2_bounds = (5154.125, 5196.625, 4895.125, 4957.625, 5301.125, 5366.125, 1)
h_gamma_a_bounds = (4319.750, 4363.500, 4283.500, 4319.750, 4367.250, 4419.750, 0)
h_delta_a_bounds = (4083.500, 4122.250, 4041.600, 4079.750, 4128.500, 4161.000, 0)
fe_4531_bounds = (4514.250, 4559.250, 4504.250, 4514.250, 4560.500, 4579.250, 0)
fe_5015_bounds = (4977.750, 5054.000, 4946.500, 4977.750, 5054.000, 5065.250, 0)
fe_5270_bounds = (5245.650, 5285.650, 5233.150, 5248.150, 5285.650, 5318.150, 0)
fe_5335_bounds = (5312.125, 5352.125, 5304.625, 5315.875, 5353.375, 5363.375, 0)
mg_b_bounds = (5160.125, 5192.625, 5142.625, 5161.375, 5191.375, 5206.375, 0)


#calculate lick indices:
#
def calc_lick_idx(alphas, alpha_stds, bounds, wavelength):

    (c1, c2, l1, l2, r1, r2, idx_type) = bounds

    idx_c1 = np.argmin(np.abs(wavelength-c1))
    idx_c2 = np.argmin(np.abs(wavelength-c2))
    idx_l1 = np.argmin(np.abs(wavelength-l1))
    idx_l2 = np.argmin(np.abs(wavelength-l2))
    idx_r1 = np.argmin(np.abs(wavelength-r1))
    idx_r2 = np.argmin(np.abs(wavelength-r2))

    center_width = c2 - c1
    variance = np.power(alpha_stds, 2)

    alphas_center = alphas[(wavelength > c1) * (wavelength < c2)]
    variance_center = variance[(wavelength > c1) * (wavelength < c2)]
    num_center = np.nansum(np.divide(alphas_center, variance_center))
    denom_center = np.nansum(np.divide(1, variance_center))
    avg_center = num_center / denom_center
    
    alphas_side = alphas[(wavelength > l1) * (wavelength < l2) + (wavelength > r1) * (wavelength < r2)]
    variance_side = variance[(wavelength > l1) * (wavelength < l2) + (wavelength > r1) * (wavelength < r2)]
    num_side = np.nansum(np.divide(alphas_side, variance_side))
    denom_side = np.nansum(np.divide(1, variance_side))
    avg_side = num_side / denom_side

    if idx_type == 0: #units: angstroms
        return (1 - avg_center / avg_side) * center_width
    else: #units: magnitudes
        return -2.5 * np.log(1 - avg_center / avg_side)

def calc_lick(alphas, alpha_stds, wavelength):
    
    h_beta = calc_lick_idx(alphas, alpha_stds, h_beta_bounds, wavelength)
    mg1 = calc_lick_idx(alphas, alpha_stds, mg1_bounds, wavelength)
    mg2 = calc_lick_idx(alphas, alpha_stds, mg2_bounds, wavelength)
    h_gamma_a = calc_lick_idx(alphas, alpha_stds, h_gamma_a_bounds, wavelength)
    h_delta_a = calc_lick_idx(alphas, alpha_stds, h_delta_a_bounds, wavelength)
    fe_4531 = calc_lick_idx(alphas, alpha_stds, fe_4531_bounds, wavelength)
    fe_5015 = calc_lick_idx(alphas, alpha_stds, fe_5015_bounds, wavelength)
    fe_5270 = calc_lick_idx(alphas, alpha_stds, fe_5270_bounds, wavelength)
    fe_5335 = calc_lick_idx(alphas, alpha_stds, fe_5335_bounds, wavelength)
    mg_b = calc_lick_idx(alphas, alpha_stds, mg_b_bounds, wavelength)
    
    mg1_fe = 0.6*mg1 + 0.4*np.log(fe_4531 + fe_5015)
    mg2_fe = 0.6*mg2 + 0.4*np.log(fe_4531 + fe_5015)
    mg_fe_p = np.sqrt(mg_b * (.72*fe_5270 + .28*fe_5335))

    #4000A discontinuity index
    #ratio of avg in 3850-3950 to 4000-4100
    idx_3850 = np.argmin(np.abs(wavelength-3850))
    idx_3950 = np.argmin(np.abs(wavelength-3950))
    idx_4000 = np.argmin(np.abs(wavelength-4000))
    idx_4100 = np.argmin(np.abs(wavelength-4100))
    left_break = np.sum(alphas[idx_3850:idx_3950])
    right_break = np.sum(alphas[idx_4000:idx_4100])
    delta = left_break / right_break

    #print("selected lines")
    #print(mg_b)
    #print(fe_5270)
    #print(fe_5335)
    #print(mg_fe_p)
    
    return np.array([h_beta, h_gamma_a, h_delta_a, mg_fe_p, mg1_fe, mg2_fe, delta])

#metric to compare index results from different spectra
def distance_metric(indices_1, indices_2):
    default_indices = [3.168, 0.972, 0.477, 2.066, 5.014, 4.02, 0.738]
    indices_1_scale = indices_1 / default_indices
    indices_2_scale = indices_2 / default_indices
    
    return np.sum(np.power(indices_1_scale - indices_2_scale, 2))


measured_indices = calc_lick(alphas[1], alpha_stds[1], wavelength)
print(measured_indices)


#calculate lick indices for the simulated spectra
#(need to interpolate?)
for p in paths:
    a = np.loadtxt(p)
    wav = a[:,0]
    values = a[:,1]
    stds = np.ones(len(wav)) #all same, no effect

    indices = calc_lick(values, stds, wav)

    if np.isnan(indices[4]):
        print("found NAN")
        print(p)

    #compare to alphas:
    measured_indices = calc_lick(alphas[1], alpha_stds[1], wavelength)
    dist = distance_metric(indices, measured_indices)

    #print out and plot the ones that are a good match
    #if dist < 10:
    #    plt.plot(wav, values, drawstyle='steps')
    #    plt.plot(wavelength, alphas[1], drawstyle='steps')
    #    plt.xlim(3500, 10000)
    #    plt.ylim(0, 1)
    #    plt.show()

    #print("\npath and distance")
    #print(p)
    #print(dist)


def linear_combo_2(coeffs, p1, p2):

    a1 = np.loadtxt(p1)
    a2 = np.loadtxt(p2)
    values1 = a1[:,1]
    values2 = a2[:,1]
    wav = a1[:,0]
    stds = np.ones(len(wav)) #all same, no effect
    
    values = coeffs[0]*values1 + coeffs[1]*values2
    indices = calc_lick(values, stds, wav)

    measured_indices = calc_lick(alphas[1], alpha_stds[1], wavelength)
    return distance_metric(indices, measured_indices)
    
    
#best combos of 2:

for p1 in paths:
    for p2 in paths:
        min_res = minimize(linear_combo_2, [1, 1], args=(p1, p2))
        if min_res.fun < 2:
            print("\nsuccess:", min_res.fun)
            print(min_res.x)
            print(p1)
            print(p2)

        


#max likelihood fitting (linear combo of spectra):
