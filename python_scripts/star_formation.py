##############################################
# SEE STAR_FORMATION_V2 FOR CODE WHERE EVERY MODEL SPECTRUM IS DIVIDED BY SAME POLYNOMIAL.
# THIS VERSION HAS LICK INDEX CODE.
##############################################


# code to compare alphas to star formation files
# in order to find best-fit spectra
# (Lick index stuff at the end)

import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import interp1d
from scipy import stats # TEST
from scipy import integrate
from generate_plots import generate_binned_alphas

import warnings
warnings.filterwarnings("ignore")  # TEST (remove eventually)

# parameters
# (poly_degree is for the polynomial added to the linear combo)
# (continuum_degree is for ???)
poly_degree = -1  # 1 is linear, etc.
poly_order = poly_degree + 1
continuum_deg = 5
continuum_order = continuum_deg + 1
bootstrap = False
num_bootstrap_samples = 200
num_spectra = 3
min_wav = 4000
max_wav = 10000  # 10000

# paths for bc03 spectra
paths_bc03 = glob.glob('/Users/blakechellew/Documents/DustProject/BrandtFiles/bc03/*.spec')

if bootstrap:
    alphas_bootstrap = np.load('../alphas_and_stds/bootstrap_alphas_boss_iris_2d_012720.npy')
    # alpha_stds_bootstrap = np.load('../alphas_and_stds/bootstrap_alpha_stds_boss_iris_2d_012720.npy')
# alphas and stds (boss):
alphas_boss = np.load('../alphas_and_stds/alphas_boss_iris_2d_012720_10.npy')
alpha_stds_boss = np.load('../alphas_and_stds/alpha_stds_boss_iris_2d_012720_10.npy')
wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')

if bootstrap:
    # calculate bootstrap errors (ignoring flux calibration):
    bootstrap_lower = np.nanpercentile(alphas_bootstrap, 16, axis=0)
    bootstrap_upper = np.nanpercentile(alphas_bootstrap, 84, axis=0)
    bootstrap_stds = (bootstrap_upper - bootstrap_lower) / 2

# TEST set everything = 1
# alphas_boss[:] = 1
# alpha_stds_boss[:] = 1

"""
# TEST: plot to see if bootstrap alphas match
for i in range(num_bootstrap_samples):
    plt.plot(wavelength, alphas_bootstrap[i], 'g')
plt.plot(wavelength, alphas_boss, 'k')
plt.ylim(0, .6)
plt.show()
"""

# TEST: 1d alphas
#alphas = np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_10.npy')
#alpha_stds = np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_10.npy')

def alphas_to_coeffs(alphas, alpha_stds, wavelength, paths, showPlots=True):

    # limit the wavelength range to 4000 to 10000 A:
    # (BOSS goes from 3550 to 10,400)
    alphas = alphas[(wavelength > min_wav) & (wavelength < max_wav)]
    alpha_stds = alpha_stds[(wavelength > min_wav) & (wavelength < max_wav)]
    wavelength = wavelength[(wavelength > min_wav) & (wavelength < max_wav)]

    # mask emission lines
    emission_line_mask = np.zeros(len(wavelength), dtype=int)
    emission_lines = [3727, 4863, 4960, 5008, 5877, 6550, 6565, 6585, 6718, 6733]
    for line in emission_lines:
        peak_idx = np.argmin(np.abs(wavelength-line))
        emission_line_mask[peak_idx-3:peak_idx+4] = 1
    alphas = alphas[np.logical_not(emission_line_mask)]
    alpha_stds = alpha_stds[np.logical_not(emission_line_mask)]
    wavelength = wavelength[np.logical_not(emission_line_mask)]

    # if only 3 model spectra:
    # cst_6gyr_z02 [idx 29]
    # t9e9_12gyr_z02 [idx 21]
    # t5e9_12gyr_z02 [idx 22]
    paths = [paths[29],
             paths[21],
             paths[22]]

    print("verify paths:")
    print(paths)

    # load the spectra and interpolate to same wavelength range
    # emission lines are effectively masked bc those wavelengths are removed
    model_spectra = np.zeros((len(paths)+poly_order, len(wavelength)))
    for i, p in enumerate(paths):
        a = np.loadtxt(p)
        wav = a[:, 0]
        values = a[:, 1]
        f = interp1d(wav, values, kind='cubic')
        model_spectra[i] = np.array([f(w) for w in wavelength])

    model_spectra_continuum = np.copy(model_spectra)

    # subtract continuum from the model spectra
    for i in range(len(paths)):
        coeffs_i = np.polyfit(wavelength, model_spectra[i], continuum_deg)
        continuum_fit = 0
        for j in range(continuum_order):
            continuum_fit += coeffs_i[-1-j] * np.power(wavelength, j)
        assert continuum_fit.shape == wavelength.shape
        model_spectra_continuum[i] = model_spectra[i] / continuum_fit

    # subtract continuum from the alphas
    coeffs_i = np.polyfit(wavelength, alphas, continuum_deg, w=1/alpha_stds)
    continuum_fit = 0
    for j in range(continuum_order):
        continuum_fit += coeffs_i[-1 - j] * np.power(wavelength, j)
    alphas_continuum = alphas / continuum_fit
    alpha_stds_continuum = alpha_stds / continuum_fit

    # TEST: plot each of the spectra alongside alphas
    colors = ['r', 'g', 'b']
    #plt.plot(wavelength, alphas_continuum, 'k', drawstyle='steps')
    for i in range(len(paths)):
        plt.plot(wavelength, model_spectra[i], colors[i], drawstyle='steps')
    plt.show()
    for i in range(len(paths)):
        plt.plot(wavelength, model_spectra_continuum[i], colors[i], drawstyle='steps')
    plt.show()

    # find best-fit linear combination using max likelihood estimator
    # and plot the best-fit spectrum against actual spectrum
    num_spectra = len(paths)+poly_order
    A = np.zeros((num_spectra, num_spectra))
    for i in range(num_spectra):
        for j in range(num_spectra):
            A[i, j] = np.nansum(model_spectra_continuum[i]*model_spectra_continuum[j]/np.power(alpha_stds_continuum, 2))
    b = np.zeros(num_spectra)
    for i in range(num_spectra):
        b[i] = np.nansum(alphas_continuum*model_spectra_continuum[i]/np.power(alpha_stds_continuum, 2))

    print("rank", np.linalg.lstsq(A, b)[2])
    coeffs = np.linalg.solve(A, b)
    # coeffs = np.linalg.lstsq(A, b)[0]
    print("coefficients")
    print(coeffs)

    # A happens to be the inverse covariance matrix.
    # (see paper notes)
    Cov = np.linalg.inv(A)
    print("covariance")
    print(Cov)

    # variance of the function combining the 3 best-fit coeffs
    # calculate separately at each wavelength
    # here the "coeffs" are value of model spectrum at given wavelength
    fitted_model_vars = np.zeros(wavelength.shape)
    for i in range(len(wavelength)):
        a = model_spectra_continuum[:, i]
        a_col = a.reshape(len(a), 1)
        a_row = a.reshape(1, len(a))
        var_i = a_row.dot(Cov).dot(a_col)
        fitted_model_vars[i] = var_i
    fitted_model_stds = np.sqrt(fitted_model_vars)

    # combine to get the best fit spectrum:
    best_fit_model = np.sum(model_spectra_continuum * coeffs[:, None], axis=0)
    best_fit_uncorrected = np.sum(model_spectra * coeffs[:, None], axis=0)

    # find the chai squared:
    print("chai squared")
    chai_squared = np.nansum(np.divide(np.power(best_fit_model - alphas_continuum, 2), np.power(alpha_stds_continuum, 2)))
    print(chai_squared)
    print("per degree of freedom (for 3 spectra)")
    chai_df = chai_squared / (len(wavelength) + poly_order)
    print(chai_df)

    # find correlation:
    corr = stats.pearsonr(alphas_continuum, best_fit_model)
    print("correlation coeff:", corr)

    if showPlots:  # if bootstrapping there could be too many plots
        plt.plot(wavelength, alphas_continuum, 'k', drawstyle='steps')
        plt.plot(wavelength, best_fit_model, 'r', drawstyle='steps')
        plt.fill_between(wavelength, best_fit_model + 3*fitted_model_stds, best_fit_model - 3*fitted_model_stds, color='r', alpha=0.5, step='pre')
        plt.plot(wavelength, alpha_stds_continuum, 'g', drawstyle='steps')
        #plt.plot(wavelength, best_fit_model + 3*fitted_model_stds, 'r', drawstyle='steps')
        #plt.plot(wavelength, best_fit_model - 3*fitted_model_stds, 'r', drawstyle='steps')
        plt.xlim(3500, 10000)
        plt.ylim(0, 2)
        plt.show()

        #try plotting without adjusting:
        plt.plot(wavelength, alphas, 'k', drawstyle='steps')
        plt.plot(wavelength, best_fit_uncorrected, 'r', drawstyle='steps')
        plt.xlim(3500, 10000)
        plt.ylim(0, 2)
        plt.show()

    print("shapes")
    return coeffs, wavelength, best_fit_uncorrected

# do some fitting
"""
if bootstrap:
    coeffs_array = np.zeros((num_bootstrap_samples, num_spectra))
    # best to use a variable here. len(wavelength) changes due due masking / cutoffs.
    best_fit_array = np.zeros((num_bootstrap_samples, 3916))
    best_fit_wavelengths = None

    for i in range(num_bootstrap_samples):
        print("sample number", i)
        # alpha_stds_boss is not the same as the bootstrap stds, but it's basically equal
        coeffs_array[i], best_fit_wavelengths, best_fit_array[i] = alphas_to_coeffs(alphas_bootstrap[i], alpha_stds_boss, wavelength, paths_bc03, showPlots=False)  # TEST: was using alphas_bootstrap[i]
    print(coeffs_array)

    # combine the coefficients:
    avg_coeffs = np.mean(coeffs_array, axis=0)
    coeff_stds = np.std(coeffs_array, axis=0)
    print("coeffs with stds:")
    print(avg_coeffs)
    print(coeff_stds)

    # bin the bootstrap stuff
    # I think the std input here can be ones, but maybe there is something better?
    print("best fit array:", best_fit_array.shape)
    lambdas_bin, best_fit_bin, _ = generate_binned_alphas([best_fit_array], [np.ones(best_fit_array.shape)], best_fit_wavelengths, bin_offset=0)
    best_fit_bin = best_fit_bin[0]
    print("best fit bin:", best_fit_bin.shape)
    # bin the boss stuff
    lambdas_bin_boss, boss_bin, _ = generate_binned_alphas([alphas_boss], [alpha_stds_boss], wavelength, bin_offset=0)
    boss_bin = boss_bin[0]
    # percentiles
    bootstrap_binned_lower = np.nanpercentile(best_fit_bin, 16, axis=0)  # 68 percent confidence interval
    bootstrap_binned_upper = np.nanpercentile(best_fit_bin, 84, axis=0)
    # TODO: need to run with the original alphas. for now take the mean.
    plt.plot(lambdas_bin, np.mean(best_fit_bin, axis=0), 'k', drawstyle='steps')
    plt.fill_between(lambdas_bin, bootstrap_binned_lower, bootstrap_binned_upper, linewidth=0.0, color='k',
                     alpha=0.2, step='pre')
    plt.plot(lambdas_bin_boss, boss_bin, 'r', drawstyle='steps')
    plt.show()

else:
    coeffs, _, _ = alphas_to_coeffs(alphas_boss, alpha_stds_boss, wavelength, paths_bc03)
    print("returned coeffs:")
    print(coeffs)
"""

################

"""
# try fitting with curve_fit instead:
def func(x, a, b, c):
    spec_1 = interp1d(wavelength, model_spectra[0], kind='cubic')
    spec_2 = interp1d(wavelength, model_spectra[1], kind='cubic')
    spec_3 = interp1d(wavelength, model_spectra[2], kind='cubic')
    return a*spec_1(x) + b*spec_2(x) + c*spec_3(x)
from scipy.optimize import curve_fit
coeffs, pcov = curve_fit(func, wavelength, alphas, bounds=(-30, 30), sigma=alpha_stds)
print("new coeffs")
print(coeffs)
"""

# LICK INDEX STUFF BELOW THIS POINT #############################################

# some important lick index tuples:
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


# calculate lick indices:
# following procedure in Worthey 1994
def calc_lick_idx(alphas, alpha_stds, bounds, wavelength, show_plots=False):

    (c1, c2, l1, l2, r1, r2, idx_type) = bounds

    center_width = c2 - c1
    variance = np.power(alpha_stds, 2)

    alphas_center = alphas[(wavelength > c1) * (wavelength < c2)]
    variance_center = variance[(wavelength > c1) * (wavelength < c2)]
    num_center = np.nansum(np.divide(alphas_center, variance_center))
    denom_center = np.nansum(np.divide(1, variance_center))
    avg_center = num_center / denom_center
    
    alphas_left = alphas[(wavelength > l1) * (wavelength < l2)]
    variance_left = variance[(wavelength > l1) * (wavelength < l2)]
    num_left = np.nansum(np.divide(alphas_left, variance_left))
    denom_left = np.nansum(np.divide(1, variance_left))
    avg_left = num_left / denom_left

    alphas_right = alphas[(wavelength > r1) * (wavelength < r2)]
    variance_right = variance[(wavelength > r1) * (wavelength < r2)]
    num_right = np.nansum(np.divide(alphas_right, variance_right))
    denom_right = np.nansum(np.divide(1, variance_right))
    avg_right = num_right / denom_right

    # define continuum line between centers of left and right continua
    left_center = (l1 + l2) / 2
    right_center = (r1 + r2) / 2
    continuum_slope = (avg_right - avg_left) / (right_center - left_center)
    continuum_intercept = avg_right - continuum_slope * right_center
    def continuum(w):
        return continuum_slope * w + continuum_intercept
    def flux_ratio(w):
        alpha_idx = np.argmin(np.abs(w - wavelength))
        return alphas[alpha_idx] / continuum(w)

    # now integrate (continuum - spectrum) / continuum
    ratio_integral, err = integrate.quad(flux_ratio, c1, c2)

    if c1 == 4319.75:
        # plot to help with interpretation:
        plt.plot(wavelength, alphas)
        plt.xlim(l1 - 20, r2 + 20)
        plt.ylim(0, 1.5)
        plt.vlines([c1, c2, l1, l2, r1, r2], 0, 1.5)
        plt.plot(wavelength, continuum(wavelength), color='red')
        plt.title(bounds)
        plt.show()

    if idx_type == 0: #units: angstroms
        # return (avg_side - avg_center) / avg_side * center_width
        equiv_width = (c2 - c1) - ratio_integral
        return equiv_width
    else: #units: magnitudes
        # return -2.5 * np.log(1 - avg_center / avg_side)
        mag = -2.5 * np.log( (1 / (c2 - c1)) * ratio_integral)
        return mag

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


measured_indices = calc_lick(alphas_boss, alpha_stds_boss, wavelength)
print("measured indices")
print(measured_indices)
plt.plot(wavelength, alphas_boss)
plt.title("BOSS alphas")
plt.ylim(0, .6)
plt.xlim(4000, 10000)
plt.show()


#calculate lick indices for the simulated spectra
#(need to interpolate?)
for p in [paths_bc03[29], paths_bc03[21]]:
    a = np.loadtxt(p)
    wav = a[:, 0]
    values = a[:, 1]
    stds = np.ones(len(wav)) #all same, no effect

    indices = calc_lick(values, stds, wav)
    print("indices for path: ", p)
    print(indices)

    if np.isnan(indices[4]):
        print("found NAN")
        print(p)

    #compare to alphas:
    measured_indices = calc_lick(alphas_boss, alpha_stds_boss, wavelength)
    dist = distance_metric(indices, measured_indices)

    #print out and plot the ones that are a good match
    #if dist < 10:
    #    plt.plot(wav, values, drawstyle='steps')
    #    plt.plot(wavelength, alphas[0], drawstyle='steps')
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

    measured_indices = calc_lick(alphas[0], alpha_stds[0], wavelength)
    return distance_metric(indices, measured_indices)
    

"""   
#best combos of 2:
for p1 in paths_bc03:
    for p2 in paths_bc03:
        min_res = minimize(linear_combo_2, [1, 1], args=(p1, p2))
        if min_res.fun < 2:
            print("\nsuccess:", min_res.fun)
            print(min_res.x)
            print(p1)
            print(p2)
"""