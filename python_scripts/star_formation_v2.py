# code to compare alphas to star formation files
# in order to find best-fit spectra.

# in this version we'll divide each spectrum by the same polynomial.

import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import interp1d
from scipy import stats
from generate_plots import generate_binned_alphas
from scipy.optimize import curve_fit

# parameters
# (poly_degree is for the polynomial added to the linear combo)
# (continuum_degree is for ???)
poly_degree = -1  # 1 is linear, etc.
poly_order = poly_degree + 1
continuum_deg = 10
alpha_continuum_deg = 10
continuum_order = continuum_deg + 1
bootstrap = False
num_bootstrap_samples = 200
num_spectra = 3
min_wav = 4000
max_wav = 9000  # bc03 switches to lower resolution around 9300

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

# TEST: 1d alphas
# alphas = np.load('../alphas_and_stds/alphas_boss_iris_1d_91119_10.npy')
# alpha_stds = np.load('../alphas_and_stds/alpha_stds_boss_iris_1d_91119_10.npy')

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
    model_spectra = np.zeros((len(paths) + poly_order, len(wavelength)))
    for i, p in enumerate(paths):
        a = np.loadtxt(p)
        wav = a[:, 0]
        values = a[:, 1]
        f = interp1d(wav, values, kind='cubic')
        model_spectra[i] = np.array([f(w) for w in wavelength])

    # construct a linear combo and divide by the continuum
    def func(x, a, b, c):
        combined_spectrum = a*model_spectra[0] + b*model_spectra[1] + c*model_spectra[2]

        # divide by the continuum
        coeffs_i = np.polyfit(wavelength, combined_spectrum, continuum_deg)
        continuum_fit = 0
        for j in range(continuum_order):
            continuum_fit += coeffs_i[-1 - j] * np.power(wavelength, j)
        assert continuum_fit.shape == wavelength.shape
        combined_spectrum_continuum = combined_spectrum / continuum_fit

        return combined_spectrum_continuum

    # subtract continuum from the alphas
    coeffs_i = np.polyfit(wavelength, alphas, alpha_continuum_deg, w=1 / alpha_stds)
    continuum_fit = 0
    for j in range(alpha_continuum_deg+1):
        continuum_fit += coeffs_i[-1 - j] * np.power(wavelength, j)
    alphas_continuum = alphas / continuum_fit
    alpha_stds_continuum = alpha_stds / continuum_fit

    coeffs, pcov = curve_fit(func, wavelength, alphas_continuum, bounds=(-30, 30), p0=(1, 1, 1), sigma=alpha_stds_continuum)
    print("coefficients")
    print(coeffs)

    # reconstruct best-fit spectrum
    best_fit_model = func(wavelength, coeffs[0], coeffs[1], coeffs[2])
    plt.plot(wavelength, alphas_continuum, 'k', drawstyle='steps')
    plt.plot(wavelength, best_fit_model, 'r', drawstyle='steps')
    plt.title("Continuum subtracted; compare to best-fit spectrum")
    plt.show()

    # plot with just one spectrum
    single_model_corrected = func(wavelength, 0, 1, 0)
    plt.plot(wavelength, alphas_continuum, 'k', drawstyle='steps', label='BOSS')
    plt.plot(wavelength, best_fit_model, 'r', drawstyle='steps', label='BC03')
    plt.title("Continnum subtracted, no fitting, compare to t9e9")
    # plt.xlim(6000, 8000)  # big wiggles
    # plt.ylim(.55, 1.45)
    plt.xlim(5025, 5450)  # small wiggles
    plt.ylim(.7, 1.25)
    plt.legend()
    plt.show()

    # TODO: multiply by the continuum to reconstruct the spectrum...
    # TODO: (maybe) plot the corrected models, but separately
    best_fit_uncorrected = best_fit_model  # TEMP

    # find the chai squared:
    print("chai squared")
    chai_squared = np.nansum(
        np.divide(np.power(best_fit_model - alphas_continuum, 2), np.power(alpha_stds_continuum, 2)))
    print(chai_squared)
    print("per degree of freedom (for 3 spectra)")
    chai_df = chai_squared / (len(wavelength) + poly_order)
    print(chai_df)

    # find correlation:
    corr = stats.pearsonr(alphas_continuum, best_fit_model)
    print("correlation coeff:", corr)

    return coeffs, wavelength, best_fit_uncorrected


if bootstrap:
    coeffs_array = np.zeros((num_bootstrap_samples, num_spectra))
    # best to use a variable here. len(wavelength) changes due due masking / cutoffs.
    best_fit_array = np.zeros((num_bootstrap_samples, 3916))
    best_fit_wavelengths = None

    for i in range(num_bootstrap_samples):
        print("sample number", i)
        # alpha_stds_boss is not the same as the bootstrap stds, but it's basically equal
        coeffs_array[i], best_fit_wavelengths, best_fit_array[i] = alphas_to_coeffs(alphas_bootstrap[i],
                                                                                    alpha_stds_boss, wavelength,
                                                                                    paths_bc03,
                                                                                    showPlots=False)  # TEST: was using alphas_bootstrap[i]
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
    lambdas_bin, best_fit_bin, _ = generate_binned_alphas([best_fit_array], [np.ones(best_fit_array.shape)],
                                                          best_fit_wavelengths, bin_offset=0)
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