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
from scipy.optimize import minimize
from scipy import stats # TEST
from scipy import integrate
from scipy import ndimage
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

# TEST: look at south only
alphas_boss = np.load('../alphas_and_stds/alphas_south011720.npy')
alpha_stds_boss = np.load('../alphas_and_stds/alpha_stds_south011720.npy')

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

def truncate_wav(w, a, s):
    # limit the wavelength range to 4000 to 10000 A:
    # (BOSS goes from 3550 to 10,400)
    alphas = a[(w > min_wav) & (w < max_wav)]
    alpha_stds = s[(w > min_wav) & (w < max_wav)]
    wavelength = w[(w > min_wav) & (w < max_wav)]
    return wavelength, alphas, alpha_stds

def mask_emission(w, a, s):
    # mask emission lines
    emission_line_mask = np.zeros(len(w), dtype=int)
    emission_lines = [3727, 4863, 4960, 5008, 5877, 6550, 6565, 6585, 6718, 6733]
    for line in emission_lines:
        peak_idx = np.argmin(np.abs(w - line))
        emission_line_mask[peak_idx - 3:peak_idx + 4] = 1
    alphas = a[np.logical_not(emission_line_mask)]
    alpha_stds = s[np.logical_not(emission_line_mask)]
    wavelength = w[np.logical_not(emission_line_mask)]
    return wavelength, alphas, alpha_stds

# mask emission for the boss alphas
wav_masked, alphas_boss_masked, alpha_stds_boss_masked = mask_emission(wavelength, alphas_boss, alpha_stds_boss)

def alphas_to_coeffs(alphas, alpha_stds, wavelength, paths, showPlots=True):
    wavelength, alphas, alpha_stds = truncate_wav(wavelength, alphas, alpha_stds)
    wavelength, alphas, alpha_stds = mask_emission(wavelength, alphas, alpha_stds)

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

# LICK INDEX STUFF #############################################

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
# following procedure is from Worthey 1994
def calc_lick_idx(alphas, alpha_stds, bounds, wavelength):

    (c1, c2, l1, l2, r1, r2, idx_type) = bounds

    center_width = c2 - c1
    variance = np.power(alpha_stds, 2)

    alphas_center = alphas[(wavelength > c1) * (wavelength < c2)]
    sigmas_center = alpha_stds[(wavelength > c1) * (wavelength < c2)]
    variance_center = variance[(wavelength > c1) * (wavelength < c2)]
    wavelength_center = wavelength[(wavelength > c1) * (wavelength < c2)]
    delta_wav_center = np.diff(wavelength_center)
    delta_wav_center = np.append(delta_wav_center, delta_wav_center[-1])
    num_center = np.sum(np.divide(alphas_center, variance_center))
    denom_center = np.sum(np.divide(1, variance_center))
    avg_center = num_center / denom_center
    
    alphas_left = alphas[(wavelength > l1) * (wavelength < l2)]
    sigmas_left = alpha_stds[(wavelength > l1) * (wavelength < l2)]
    variance_left = variance[(wavelength > l1) * (wavelength < l2)]
    num_left = np.sum(np.divide(alphas_left, variance_left))
    denom_left = np.sum(np.divide(1, variance_left))
    avg_left = num_left / denom_left

    alphas_right = alphas[(wavelength > r1) * (wavelength < r2)]
    sigmas_right = alpha_stds[(wavelength > r1) * (wavelength < r2)]
    variance_right = variance[(wavelength > r1) * (wavelength < r2)]
    num_right = np.sum(np.divide(alphas_right, variance_right))
    denom_right = np.sum(np.divide(1, variance_right))
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

    # error bars (see iPad notes)
    center_term = np.sum(np.power(delta_wav_center * sigmas_center / continuum(wavelength_center), 2))
    dm_dalph_left = (1 / variance_left) / (right_center - left_center) / np.sum(1 / variance_left)
    dm_dalph_right = (1 / variance_right) / (right_center - left_center) / np.sum(1 / variance_right)
    dR_dalph_left = dm_dalph_left * np.sum(delta_wav_center * alphas_center / continuum(wavelength_center)**2 * np.abs(wavelength_center - right_center))
    dR_dalph_right = dm_dalph_right * np.sum(delta_wav_center * alphas_center / continuum(wavelength_center)**2 * np.abs(wavelength_center - left_center))
    left_term = np.sum(np.power(sigmas_left * dR_dalph_left, 2))
    right_term = np.sum(np.power(sigmas_right * dR_dalph_right, 2))
    ratio_err = np.sqrt(center_term + left_term + right_term)

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
        equiv_width = (c2 - c1) - ratio_integral
        return equiv_width, ratio_err
    else: #units: magnitudes
        # return -2.5 * np.log(1 - avg_center / avg_side)
        mag = -2.5 * np.log( (1 / (c2 - c1)) * ratio_integral)
        mag_err = 2.5 / ratio_integral * ratio_err
        return mag, mag_err

# we are ignoring the balmer series because there are emission lines there in the observed data.
def calc_lick(alphas, alpha_stds, wavelength):
    
    # h_beta = calc_lick_idx(alphas, alpha_stds, h_beta_bounds, wavelength)
    mg1, mg1_err = calc_lick_idx(alphas, alpha_stds, mg1_bounds, wavelength)
    mg2, mg2_err = calc_lick_idx(alphas, alpha_stds, mg2_bounds, wavelength)
    # h_gamma_a = (alphas, alpha_stds, h_gamma_a_bounds, wavelength)
    # h_delta_a = calc_lick_idx(alphas, alpha_stds, h_delta_a_bounds, wavelength)
    fe_4531, fe_4531_err = calc_lick_idx(alphas, alpha_stds, fe_4531_bounds, wavelength)
    fe_5015, fe_5015_err = calc_lick_idx(alphas, alpha_stds, fe_5015_bounds, wavelength)
    fe_5270, fe_5270_err = calc_lick_idx(alphas, alpha_stds, fe_5270_bounds, wavelength)
    fe_5335, fe_5335_err = calc_lick_idx(alphas, alpha_stds, fe_5335_bounds, wavelength)
    mg_b, mg_b_err = calc_lick_idx(alphas, alpha_stds, mg_b_bounds, wavelength)

    mg1_fe = 0.6*mg1 + 0.4*np.log(fe_4531 + fe_5015)
    mg1_fe_err = np.sqrt(.36 * mg1_err**2 + (.16 / (fe_4531 + fe_5015))**2 * (fe_4531_err**2 + fe_5015_err**2))
    mg2_fe = 0.6*mg2 + 0.4*np.log(fe_4531 + fe_5015)
    mg2_fe_err = np.sqrt(.36 * mg2_err ** 2 + (.16 / (fe_4531 + fe_5015)) ** 2 * (fe_4531_err ** 2 + fe_5015_err ** 2))
    mg_fe_p = np.sqrt(mg_b * (.72*fe_5270 + .28*fe_5335))
    mg_fe_p_err = 1 / 2 / np.sqrt(mg_fe_p) * np.sqrt((.72*fe_5270 + .28*fe_5335)**2 * mg_b_err**2 + (.72**2*fe_5270_err**2 + .28**2*fe_5335_err**2) * mg_b**2)

    #4000A discontinuity index
    #ratio of avg in 3850-3950 to 4000-4100
    idx_3850 = np.argmin(np.abs(wavelength-3850))
    idx_3950 = np.argmin(np.abs(wavelength-3950))
    idx_4000 = np.argmin(np.abs(wavelength-4000))
    idx_4100 = np.argmin(np.abs(wavelength-4100))
    left_break = np.sum(alphas[idx_3850:idx_3950])
    left_errs = alpha_stds[idx_3850:idx_3950]
    right_break = np.sum(alphas[idx_4000:idx_4100])
    right_errs = alpha_stds[idx_4000:idx_4100]
    delta = left_break / right_break
    delta_err = np.sqrt(1 / right_break**2 * np.sum(left_errs**2) + (left_break**2 / right_break**4) * np.sum(right_errs**2))

    return_indices = np.array([mg_fe_p, mg1_fe, mg2_fe, delta])
    return_errors = np.array([mg_fe_p_err, mg1_fe_err, mg2_fe_err, delta_err])

    return (return_indices, return_errors)

def chi_squared(measured_indices, measured_errors, model_indices):
    # TEST: remove the first one
    measured_indices = measured_indices[1:]
    measured_errors = measured_errors[1:]
    model_indices = model_indices[1:]

    return np.sum((measured_indices - model_indices)**2 / measured_errors**2)

# get indices from measured alphas
measured_indices, measured_errs = calc_lick(alphas_boss, alpha_stds_boss, wavelength)
print("measured indices")
print(measured_indices)
print(measured_errs)

# calculate lick indices for simulated spectra
# (no radiative transfer)
for p in [paths_bc03[21], paths_bc03[22], paths_bc03[29]]:
# for p in paths_bc03:
    print("path:")

    a = np.loadtxt(p)
    wav = a[:, 0]
    values = a[:, 1]
    stds = np.ones(len(wav))  # all same, no effect

    indices, _ = calc_lick(values, stds, wav)
    print("measured vs. indices for path: ", p)
    print(measured_indices)
    print(indices)
    print("chi squared and per df")
    chi_sqr = chi_squared(measured_indices, measured_errs, indices)
    print(chi_sqr)
    print(chi_sqr / 4)

print("\n")

# now with radiative transfer
# (the spectrum here was already interpolated to the BOSS wavelengths)
# order: t9e9, t5e9, CST
bc03_rads =[np.load('/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/t9e9_12gyr_z02zd082221.npy'),
           np.load('/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/t5e9_12gyr_z02zd072121.npy'),
           np.load('/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/cst_6gyr_z02zd072121.npy')]
bc03_rads_masked = [mask_emission(wavelength, a, np.ones(len(wavelength)))[1] for a in bc03_rads]
stds = np.ones(len(wavelength))
for bc03_rad in bc03_rads:
    rad_indices, _ = calc_lick(bc03_rad, stds, wavelength)
    chi_sqr = chi_squared(measured_indices, measured_errs, rad_indices)
    print("chi squared and per df (radiative)")
    print(chi_sqr)
    print(chi_sqr / 4)
    print("measured indices and model:")
    print(measured_indices)
    print(rad_indices)

print("\n")

# test smoothing
boss_wavelength_7000 = wavelength[2955]
boss_diff_7000 = np.diff(wavelength)[2955]
boss_dlambda_over_lambda = boss_diff_7000 / boss_wavelength_7000  # should be .0002306

#interlude: plot the smoothed spectra
def top_hat_smooth(model_spectra, width):

    width_A = int(width / boss_diff_7000)
    print("width", width_A)
    boss_smoothed = ndimage.uniform_filter1d(alphas_boss_masked, width_A)
    boss_norm = alphas_boss_masked / boss_smoothed

    for i, model_spectrum in enumerate(model_spectra):
        model_smoothed = ndimage.uniform_filter1d(model_spectrum, width_A)
        model_norm = model_spectrum / model_smoothed
        # plt.plot(wav_masked, model_spectrum, label='bc03')
        # plt.plot(wav_masked, model_smoothed, label='bc03 smooth')
        plt.plot(wav_masked, model_norm, label='bc03 norm' + str(i))

    plt.plot(wav_masked, alphas_boss_masked, label='BOSS')
    plt.plot(wav_masked, boss_smoothed, label='BOSS smooth')
    plt.plot(wav_masked, boss_norm, label='BOSS norm')
    plt.ylim(0, 2)
    # plt.xlim(3850, 4100)  # TEST
    # plt.xlim(4945, 5065)  # TEST
    plt.legend()
    plt.show()


top_hat_smooth(bc03_rads_masked, 400)

exit(0)

# convolve with a gaussian
# v_dispersion should be in km/s
def test_smoothing(model_spectrum, v_dispersion=100, plot=False):
    sigma_over_lambda = v_dispersion / (3 * 10**5)
    sigma_pixels = sigma_over_lambda / boss_dlambda_over_lambda
    model_smoothed = ndimage.gaussian_filter(model_spectrum, sigma_pixels, mode='nearest')

    if plot:
        # plot to make sure it was smoothed
        plt.plot(wavelength, bc03_rads[1], label='original')
        plt.plot(wavelength, model_smoothed, label='smooth')
        plt.title("Dispersion: " + str(v_dispersion) + " km/s")
        plt.legend()
        plt.show()

    stds = np.ones(len(wavelength))
    rad_indices_smoothed, _ = calc_lick(model_smoothed, stds, wavelength)
    chi_sqr = chi_squared(measured_indices, measured_errs, rad_indices_smoothed)
    return chi_sqr, chi_sqr/4, rad_indices_smoothed

# disps = [70, 80, 90, 100, 110, 120, 130, 150, 170, 190, 220, 250, 280]
disps = [100]
chi_sqrds = []
for disp in disps:
    chi_sqr, per_df, indices = test_smoothing(bc03_rads[1], disp, plot=True)
    print("dispersion:", disp)
    print("chi sqr", chi_sqr)
    print("per df", per_df)
    print("indices")
    print(indices)
    chi_sqrds = np.append(chi_sqrds, per_df)

# plot chi^2 as fn of degree of smoothing
plt.plot(disps, chi_sqrds)
plt.title("Chi Squared Per Df")
plt.show()

print("\n")

# now compare models with radiative AND smoothing:
disp_fiducial = 150
for bc03_rad in bc03_rads:
    chi_sqr, per_df, indices = test_smoothing(bc03_rad, disp_fiducial)

    print("chi squared and per df (radiative and smoothed)")
    print(chi_sqr)
    print(chi_sqr / 4)
    print("measured indices and model:")
    print(measured_indices)
    print(indices)

# find best-fit coefficients for a linear combo of 2 model spectra
def linear_combo_2(coeffs, p1, p2):

    a1 = np.loadtxt(p1)
    a2 = np.loadtxt(p2)
    values1 = a1[:,1]
    values2 = a2[:,1]
    wav = a1[:,0]
    stds = np.ones(len(wav)) #all same, no effect
    
    values = coeffs[0]*values1 + coeffs[1]*values2
    indices, index_errs = calc_lick(values, stds, wav)

    return chi_squared(measured_indices, measured_errs, indices)


# try some linear combos
"""
print("\nTry some linear combos")
p1 = paths_bc03[29]
p2 = paths_bc03[22]
chi2 = linear_combo_2([0, 1], p1, p2)
print(str(chi2) + " for coeffs 0")
chi2 = linear_combo_2([.1, .9], p1, p2)
print(str(chi2) + " for coeffs 0.1")
chi2 = linear_combo_2([.2, .8], p1, p2)
print(str(chi2) + " for coeffs 0.2")
chi2 = linear_combo_2([.3, .7], p1, p2)
print(str(chi2) + " for coeffs 0.3")
chi2 = linear_combo_2([.4, .6], p1, p2)
print(str(chi2) + " for coeffs 0.4")
chi2 = linear_combo_2([.5, .5], p1, p2)
print(str(chi2) + " for coeffs 0.5")
chi2 = linear_combo_2([.6, .4], p1, p2)
print(str(chi2) + " for coeffs 0.6")
chi2 = linear_combo_2([.7, .3], p1, p2)
print(str(chi2) + " for coeffs 0.7")
chi2 = linear_combo_2([.8, .2], p1, p2)
print(str(chi2) + " for coeffs 0.8")
chi2 = linear_combo_2([.9, .1], p1, p2)
print(str(chi2) + " for coeffs 0.9")
chi2 = linear_combo_2([1, 0], p1, p2)
print(str(chi2) + " for coeffs 1")
"""

# find best-fit coeffs for const vs. exponential
# (takes a long time)
"""
min_res = minimize(linear_combo_2, [1, 1], args=(p1, p2))
print(min_res)
"""


#best combos of 2:
"""   
for p1 in paths_bc03:
    for p2 in paths_bc03:
        min_res = minimize(linear_combo_2, [1, 1], args=(p1, p2))
        if min_res.fun < 2:
            print("\nsuccess:", min_res.fun)
            print(min_res.x)
            print(p1)
            print(p2)
"""

# print out and plot the ones that are a good match
# if dist < 10:
#    plt.plot(wav, values, drawstyle='steps')
#    plt.plot(wavelength, alphas[0], drawstyle='steps')
#    plt.xlim(3500, 10000)
#    plt.ylim(0, 1)
#    plt.show()

# print("\npath and distance")
# print(p)
# print(dist)