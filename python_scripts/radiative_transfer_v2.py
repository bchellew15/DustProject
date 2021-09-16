# take the BC03 model spectra and apply radiative transfer calculations.
# based on the appendix of BD12.

import glob
import numpy as np
from scipy.interpolate import interp1d
import scipy.special as special
import matplotlib.pyplot as plt
from generate_plots import generate_binned_alphas
from astropy.io import fits
import pickle

save = True  # save ERE figure
boss = True  # whether to interpolate to BOSS wavelengths or SDSS
wd01_model = False  # otherwise zda04 model
b = 40 * np.pi / 180  # latitude (should be 40 degrees -> radians)
sig_dust = 250  # 250 parsecs (from density eqn)   # TEST: change to 1 pc
tau_def = 0.15  # fiducial value (see BD12 paper) (z = 0)
V_band_wav = 5510  # A
# stellar scale heights 300 pc and 1350 pc
a_300 = 1  # 0.9
a_1350 = 0  # 0.1
#which bc03 model to use:
# 22 is t5e9, z=.02  (2 is z=.008)
# 21 is t9e9, 29 is const
p_num = 29
sig_star_1 = 300  # 300
sig_star_2 = 1350  # 1350
mistakes = False  # TEST
rho_scale = 1000  # 5 or 1000  #.1 matches brandt code for sig = 1

# number of grid points (for TIR and scattering):
n_lamb = 200  # not used currently
# n_beta needs to be even
n_u = 20  # 100
n_tau = 20  # 25
n_beta = 20  # should be even to avoid beta=pi, but hasn't been an issue since
n_theta = 20  # 25

paths = glob.glob('/Users/blakechellew/Documents/DustProject/BrandtFiles/bc03/*.spec')  # bc03 model spectra

# load wavelengths
wavelength_boss = np.load('../alphas_and_stds/wavelength_boss.npy')  # angstroms
hdulist_direc = '/Users/blakechellew/Documents/DustProject/BrandtFiles/'
hdulist_sdsswav = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
wavelength_sdss = np.array(hdulist_sdsswav[1].data)
if boss:
    wavelength = wavelength_boss
else:
    wavelength = wavelength_sdss

# load dust models
if wd01_model:
    dust_model_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/kext_albedo_WD_MW_3.1_60_D03.all'
    skip_rows = 80
    use_cols = (0, 1, 2, 3)
    # column 0: wavelength (microns, but converted to angstroms below)
    # column 1: albedo
    # column 2: <cos>
    # column 3: cross section (cm^2)
else:  # zda04 dust model
    dust_model_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/brandt_radiative/integrals/extcurv.out_ZDA04'
    skip_rows = 15
    use_cols = (0, 4, 3, 1)
dust_model = np.loadtxt(dust_model_path, skiprows=skip_rows, usecols=use_cols)

# flip because wavelength goes largest to smallest
dust_wav = np.flip(dust_model[:, 0].flatten()) * 10 ** 4  # convert to angstroms from um
dust_albedo = np.flip(dust_model[:, 1].flatten())
dust_cos = np.flip(dust_model[:, 2].flatten())
dust_cross = np.flip(dust_model[:, 3].flatten())  # cm^2
# sort by wavelength; some out of order:
sorted_idx = np.argsort(dust_wav)
dust_wav = dust_wav[sorted_idx]
dust_albedo = dust_albedo[sorted_idx]
dust_cos = dust_cos[sorted_idx]
dust_cross = dust_cross[sorted_idx]

if wd01_model:
    # remove duplicate wavelength
    bad_idx = 125
    dust_wav = np.delete(dust_wav, bad_idx)
    dust_albedo = np.delete(dust_albedo, bad_idx)
    dust_cos = np.delete(dust_cos, bad_idx)
    dust_cross = np.delete(dust_cross, bad_idx)

# interpolate dust model
dust_albedo_f = interp1d(dust_wav, dust_albedo, kind='cubic')
dust_cos_f = interp1d(dust_wav, dust_cos, kind='cubic')
dust_cross_f = interp1d(dust_wav, dust_cross, kind='cubic')  # cm^2

cross_sec_V = dust_cross_f(V_band_wav)  # cm^2

# set wavelength range
if wd01_model:  # 80 header lines
    wav_min, wav_max = len(dust_wav) - (631-80), len(dust_wav) - (429-80),
else:  # 15 header lines
    wav_min, wav_max = len(dust_wav) - (1023-15), len(dust_wav) - (690-15)
dust_wav = dust_wav[wav_min:wav_max]  # test. avoids interpolation error for I_tir currently.

# phase function
def henyey(cos_xi, g):
    if mistakes:
        g *= -1  # TEST
    result = (1 - g ** 2) / (1 + g ** 2 - 2 * g * cos_xi) ** (3 / 2) / (4 * np.pi)
    assert result.shape == cos_xi.shape
    return result

# eqn A1
# z should be in pc
# density: rho ~ exp(-z^2 / 2 / sig_dust^2).
# I took sigma out of the integral and integrated analytically.
# units of cross section should cancel out b/c it's also in the density prefactor.
def tau_f(lamb, z):
    cross_section = dust_cross_f(lamb)  # cm^2
    tau_0 = tau_def * cross_section / cross_sec_V
    if mistakes:
        tau_lambda = tau_0 * (1 + special.erf(z / sig_dust / np.sqrt(2)))
    else:
        tau_lambda = tau_0 * (1 - special.erf(z / sig_dust / np.sqrt(2)))
    return tau_lambda  # unitless


# reverse: find z given tau
# (see paper calculations)
# WARNING: it will multiply cross section across the last axis
def z_of_tau(lamb, tau):
    cross_section = dust_cross_f(lamb)  # cm^2
    tau_0 = tau_def * cross_section / cross_sec_V
    result = sig_dust * np.sqrt(2) * special.erfinv(1 - tau / tau_0)
    return result


# A2
def A_f(lamb, z, rho, beta):
    term1 = tau_f(lamb, z) - tau_f(lamb, z - rho * np.cos(beta))
    term2 = np.cos(beta)
    result = np.abs(term1 / term2)
    return result


# eqn A5
def cos_xi(z, rho, theta, beta):
    numer = rho * np.cos(beta)**2 - z * np.sin(beta) * np.cos(theta) / np.tan(b)
    denom_sqr = z**2 / np.tan(b)**2 + rho**2 * np.cos(beta)**2
    result = numer / np.sqrt(denom_sqr)

    beta_mask = np.abs(beta - np.pi/2) < 0.1
    if (beta_mask).any() and (z[beta_mask] == 0).any():
        print("BETA is not allowed")
        print(beta[beta_mask][z[beta_mask] == 0])
        print(result[beta_mask][z[beta_mask] == 0])
    if mistakes:
        return np.abs(result)  # TEST
    return -result



# surface power density (arbitrary units; any normalization factor will cancel out later)
# it's a function of height z_s
# the bc03 flux will be multiplied later
def surface_power_fn(z, rho, beta):  # z_s in pc
    return a_300 * np.exp(-np.abs(z - rho * np.cos(beta)) / sig_star_1) \
        + a_1350 * np.exp(-np.abs(z - rho * np.cos(beta)) / sig_star_2)


# eqn A4: integrand
# units: angstroms * units of sigma after integration
# broadcasting: lamb, z, rho, beta should all be 4d arrays
def i_tir_integrand(lamb, z, rho, beta):
    prefactor = 1 / (8 * np.pi) / np.sin(np.abs(b))
    term1 = 1 - dust_albedo_f(lamb)
    term2 = surface_power_fn(z, rho, beta)
    term3 = np.exp(-A_f(lamb, z, rho, beta))
    term4 = np.sin(beta)
    result = prefactor * term1 * term2 * term3 * term4
    result = result * rho_scale * np.exp(rho/rho_scale)  # transform rho to u
    return result

# eqn A4: integrate over everything except wavelength
def i_tir(bc03_f):

    # lamb_min = 100
    # lamb_max = 15*10**4   # 10 ** 6
    # h_lamb = (lamb_max - lamb_min) / 2 / n_lamb
    # lamb_grid = np.linspace(lamb_min + h_lamb, lamb_max - h_lamb, n_lamb)  # 91 to 10**4 A  # min, max, n
    lamb_grid = dust_wav  # TEST, for using lambda grid from dust model
    n_lamb = len(lamb_grid)

    tau_min = 0
    tau_max = tau_f(lamb_grid, 0)
    h_tau = (tau_max - tau_min) / 2 / n_tau
    # broadcasting: dim. of tau grid are (lambda, tau)
    tau_grid = np.linspace(tau_min + h_tau, tau_max - h_tau, n_tau)  # same as z, 0 to inf
    beta_min = 0
    beta_max = np.pi
    h_beta = (beta_max - beta_min) / 2 / n_beta
    beta_grid = np.linspace(beta_min + h_beta, beta_max - h_beta, n_beta)

    u_min = 0
    u_max = 1
    h_u = (u_max - u_min) / 2 / n_u
    u_grid = np.linspace(u_min + h_u, u_max - h_u, n_u)

    # meshgrid from the 3 easy dimensions (lambda, beta, rho)
    lambs, betas, us = np.meshgrid(lamb_grid, beta_grid, u_grid, indexing='ij')  # so input dims match output

    # broadcasting test
    assert lambs.shape == (n_lamb, n_beta, n_u)
    assert betas.shape == (n_lamb, n_beta, n_u)
    assert us.shape == (n_lamb, n_beta, n_u)
    assert tau_grid.shape == (n_tau, n_lamb)

    # add a new dimension for tau:
    new_shape = tuple([n_tau] + list(lambs.shape))  # lambs, betas, rhos have same shape
    zeros = np.zeros(new_shape)
    lambs = lambs[None, ...] + zeros
    betas = betas[None, ...] + zeros
    us = us[None, ...] + zeros
    taus = tau_grid[..., None, None] + zeros
    rhos = -rho_scale*np.log(us)
    assert lambs.shape == (n_tau, n_lamb, n_beta, n_u)

    ww = i_tir_integrand(lambs, z_of_tau(lambs, taus), rhos, betas)

    # integrate on simple grid
    # lamb_div = (lamb_max - lamb_min) / n_lamb
    beta_div = (beta_max - beta_min) / n_beta
    u_div = (u_max - u_min) / n_u
    tau_div = (tau_max - tau_min) / n_tau
    tau_div = tau_div[..., None, None] + zeros

    # temp, for using lambda grid from dust model
    lamb_div = np.diff(lamb_grid)
    lamb_div = np.append(lamb_div, lamb_div[-1])
    lamb_div = lamb_div[None, :, None, None]

    # multiply by bc03 and integrate over lambda
    bc03s = bc03_f(lamb_grid)[None, :, None, None]
    assert bc03s.shape == (1, ww.shape[1], 1, 1)
    result = np.sum(ww * bc03s * tau_div * beta_div * u_div * lamb_div)  # units of angstroms * sigma
    return result


# eqn A7, part 1: just the integrand.
# note that we integrate over z_s also.
# units: parsecs * units of sigma
def i_sca_integrand(theta, tau, rho, beta, lamb):
    prefactor = (1 / np.sin(np.abs(b))) * dust_albedo_f(lamb) / (4 * np.pi)
    z = z_of_tau(lamb, tau)

    term1 = np.exp((-1 / np.sin(np.abs(b))) * (tau_f(lamb, 0) - tau))
    term2 = henyey(cos_xi(z, rho, theta, beta), dust_cos_f(lamb))
    term3 = surface_power_fn(z, rho, beta)
    term4 = np.exp(-A_f(lamb, z, rho, beta))
    term5 = np.sin(beta)
    result = prefactor * term1 * term2 * term3 * term4 * term5
    result = result * rho_scale * np.exp(rho / rho_scale)  # transform rho to u

    return result


# eqn A7, part 2: do the integral for given wavelength
# still needs to be multiplied by bc03 flux
def i_sca(lamb):

    print("scattering integral:", lamb)  # to keep track of how long the code has left to run

    theta_min = 0
    theta_max = 2 * np.pi
    h_theta = (theta_max - theta_min) / 2 / n_theta
    theta_grid = np.linspace(theta_min + h_theta, theta_max - h_theta, n_theta)
    tau_min = 0
    tau_max = tau_f(lamb, 0)
    h_tau = (tau_max - tau_min) / 2 / n_tau
    tau_grid = np.linspace(tau_min + h_tau, tau_max - h_tau, n_tau)
    beta_min = 0
    beta_max = np.pi
    h_beta = (beta_max - beta_min) / 2 / n_beta
    beta_grid = np.linspace(beta_min + h_beta, beta_max - h_beta, n_beta)
    u_min = 0
    u_max = 1
    h_u = (u_max - u_min) / 2 / n_u
    u_grid = np.linspace(u_min + h_u, u_max - h_u, n_u)

    """
    # plot integrand for I_sca
    # start with placeholder values
    n_plot_pts = 50
    theta_temp = np.pi
    tau_temp = tau_f(lamb, 0) / 3
    u_temp = 0.5
    beta_temp = np.pi / 4
    # plot vs. u
    h_u_temp = (u_max - u_min) / 2 / n_plot_pts
    us_temp = np.linspace(u_min + h_u_temp, u_max - h_u_temp, n_plot_pts)
    integrand_vals = i_sca_integrand(theta_temp, tau_temp, -rho_scale * np.log(us_temp), beta_temp, lamb)
    print("us", integrand_vals)
    assert np.count_nonzero(np.isnan(integrand_vals)) == 0
    plt.plot(us_temp, integrand_vals, '.')
    plt.title('I_sca integrand vs. e^(-rho / ' + str(rho_scale) + ' pc)')
    plt.show()

    exit(0)
    """

    """
    h_theta_temp = (theta_max - theta_min) / 2 / n_plot_pts
    h_tau_temp = (tau_max - tau_min) / 2 / n_plot_pts
    h_beta_temp = (beta_max - beta_min) / 2 / n_plot_pts
    # plot vs. theta
    thetas_temp = np.linspace(theta_min + h_theta_temp, theta_max - h_theta_temp, n_plot_pts)
    integrand_vals = i_sca_integrand(thetas_temp, tau_temp, -1000 * np.log(u_temp), beta_temp, lamb, bc03_f)
    print("thetas", integrand_vals)
    assert np.count_nonzero(np.isnan(integrand_vals)) == 0
    plt.plot(thetas_temp, integrand_vals, '.')
    plt.title('I_sca integrand vs. theta')
    plt.show()
    # plot vs. tau
    taus_temp = np.linspace(tau_min + h_tau_temp, tau_max - h_tau_temp, n_plot_pts)
    integrand_vals = i_sca_integrand(theta_temp, taus_temp, -1000 * np.log(u_temp), beta_temp, lamb, bc03_f)
    print("taus", integrand_vals)
    assert np.count_nonzero(np.isnan(integrand_vals)) == 0
    plt.plot(taus_temp, integrand_vals, '.')
    plt.title('I_sca integrand vs. tau')
    plt.show()
    # plot vs. beta
    betas_temp = np.linspace(beta_min + h_beta_temp, beta_max + h_beta_temp, n_plot_pts)
    integrand_vals = i_sca_integrand(theta_temp, tau_temp, -1000 * np.log(u_temp), betas_temp, lamb, bc03_f)
    print("betas", integrand_vals)
    assert np.count_nonzero(np.isnan(integrand_vals)) == 0
    plt.plot(betas_temp, integrand_vals, '.')
    plt.title('I_sca integrand vs. beta')
    plt.show()
    # plot vs. lambda (same every time)
    lambs_temp = np.linspace(200, 12000, 150)
    integrand_vals = i_sca_integrand(theta_temp, tau_temp, -1000 * np.log(u_temp), beta_temp, lambs_temp, bc03_f)
    print("lambs", integrand_vals)
    assert np.count_nonzero(np.isnan(integrand_vals)) == 0
    plt.plot(lambs_temp, integrand_vals, '.')  # / (bc03_f(lambs_temp)
    plt.title('I_sca integrand vs. lambda')
    plt.show()
    print("shapes:")
    assert type(theta_temp) != np.ndarray
    assert type(tau_temp) != np.ndarray
    assert type(u_temp) != np.ndarray
    assert type(beta_temp) != np.ndarray
    print(lambs_temp.shape)
    """

    taus, betas, us, thetas = np.meshgrid(tau_grid, beta_grid, u_grid, theta_grid, indexing='ij')
    rhos = -rho_scale*np.log(us)
    ww = i_sca_integrand(thetas, taus, rhos, betas, lamb)

    # sum over grid
    theta_div = (theta_max - theta_min) / n_theta
    u_div = (u_max - u_min) / n_u
    tau_div = (tau_max - tau_min) / n_tau
    beta_div = (beta_max - beta_min) / n_beta
    result = np.sum(ww) * theta_div * u_div * tau_div * beta_div
    return result


# idx 21: t=9e9, 12 Gyr, z=.02 (.02 is closest to "solar" metallicity)
# idx 22: t=5e9, 12 Gyr, z=.02
# idx 29: cst, 6 Gyr, z=.02

# index 2: t=5e9, 12 Gyr, Z=.008
# index 34: t=9e9, 12 Gyr, Z=.008
# 5th to last: similar but t9e9

# for i, p in enumerate(paths):
#     print(i)
#     print(p)

# loop through the bc03 spectra
for p in paths[p_num:p_num+1]:

    ############################
    # Now the radiative transfer calculations
    ##########################

    print("verify path name:", p)

    # load and interpolate the bc03 spectra
    a = np.loadtxt(p)
    wav = a[:, 0]  # angstroms (91 A to 1.6 * 10^6)
    bc03 = a[:, 1]
    bc03_f = interp1d(wav, bc03, kind='cubic')

    """
    # plots of integrand for I_TIR: (lambda, z, rho, beta)
    num_div = 50
    tau_grid = np.linspace(0, tau_f(6000, 0), num_div)  # z_s grid, -inf to inf
    beta_grid = np.linspace(0, np.pi, num_div)
    u_grid = np.linspace(0, 1, num_div)
    lamb_grid = np.linspace(100, 15*10**4, num_div)

    # order: lambda, z, rho, beta
    lamb_plot = i_tir_integrand(lamb_grid, 1, 1, np.pi/4, bc03_f)
    print(lamb_plot)
    plt.plot(lamb_grid, lamb_plot)
    plt.title('I_TIR integrand vs. lambda')
    plt.show()
    tau_plot = i_tir_integrand(6000, z_of_tau(6000, tau_grid), 1, np.pi/4, bc03_f)
    print(tau_plot)
    plt.plot(tau_grid, tau_plot)
    plt.title('I_TIR integrand vs. tau')
    plt.show()
    beta_plot = i_tir_integrand(6000, 1, 1, beta_grid, bc03_f)
    print(beta_plot)
    plt.plot(beta_grid, beta_plot)
    plt.title('I_TIR integrand vs. beta')
    plt.show()
    u_plot = i_tir_integrand(6000, 1, -np.log(u_grid), np.pi / 4, bc03_f)
    print(u_plot)
    plt.plot(u_grid, u_plot)
    plt.title('I_TIR integrand vs. e^(-rho / 1000 pc)')
    plt.show()
    """

    # compute total infrared radiation
    # (range of lambda for dust model is 10^-4 to 10^4 microns, or 1 to 10^8 angstroms)
    # (range of lambda for bc03 is 91 to 1.6e6)

    I_tir = i_tir(bc03_f)
    print("I_tir result")
    print(I_tir)

    # convert to 100 micron (there is an associated uncertainty)
    nu_I_nu_100 = .52 * I_tir  # units of angstroms * sigma
    assert nu_I_nu_100.ndim == 0

    print("starting integral")

    wavelength_partial = dust_wav  # use wavelength grid from dust model file
    # BOSS is 3600 to 10000
    wavelength_partial = wavelength_partial[(wavelength_partial > 3400) & (wavelength_partial < 10500)]
    # wavelength_partial = wavelength[::40]
    # wavelength_partial = np.append(wavelength_partial, wavelength[-1])  # to avoid interpolation issues for SDSS array

    # wavelength_partial = wavelength

    # try to use broadcasting for this?
    i_sca_array = np.array([i_sca(lamb) for lamb in wavelength_partial])  # units of sigma * parsecs

    # interpolate the result of the scattering integral
    i_sca_f = interp1d(wavelength_partial, i_sca_array, kind='cubic')

    # then calculate alphas
    i_lam_array = bc03_f(wavelength) * i_sca_f(wavelength) * wavelength  # units of sigma * parsecs * angstroms
    alphas = i_lam_array / nu_I_nu_100
    assert alphas.shape == wavelength.shape

    print("wavelength and alphas")
    print(wavelength_partial)
    print(alphas)

    # scaling factors
    bd12_factor = 0.49 if wd01_model else 0.52
    boss_fluxfactor = 1.38
    sdss_fluxfactor = 1.38

    avg_bc03 = np.mean(bc03_f(wavelength_partial))
    avg_bc03_wav = np.mean(bc03_f(wavelength_partial) * wavelength_partial)

    # load bd12 plot for comparison
    if wd01_model:
        bd12_alphas = np.loadtxt('/Users/blakechellew/Documents/DustProject/alphas_and_stds/bd12_fig3_green_052621.csv',
                                 delimiter=",")
    else:
        bd12_alphas = np.loadtxt('/Users/blakechellew/Documents/DustProject/alphas_and_stds/bd12_fig3_blue_052621.csv',
                                 delimiter=",")

    # load plot from brandt code
    brandt_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/brandt_radiative/integrals/'
    if wd01_model:
        brandt_alphas = np.load(brandt_path + 'alphas_wd.npy')
    else:
        brandt_alphas = np.load(brandt_path + 'alphas_zd.npy')

    # bin some alphas
    lambdas_bin, alphas_bin, _ = generate_binned_alphas([alphas], [np.ones(alphas.shape)], wavelength, bin_offset=0)
    alphas_bin = alphas_bin[0]


    alphas_norad = [np.load('../alphas_and_stds/alphas_boss_iris_2d_012720.npy'),
                    np.load('../alphas_and_stds/alphas_north011720.npy'),
                    np.load('../alphas_and_stds/alphas_south011720.npy')]
    # apply correction factor
    correction_factors = [np.load('../alphas_and_stds/correction_factor_boss_iris_smooth.npy'),
                          np.load('../alphas_and_stds/correction_factor_boss_iris_north_smooth.npy'),
                          np.load('../alphas_and_stds/correction_factor_boss_iris_south_smooth.npy')]
    correction_factor_sdss = np.load('../alphas_and_stds/correction_factor_sdss_iris_smooth.npy')
    _, binned_corrections, _ = generate_binned_alphas(correction_factors, 3 * [np.ones(len(correction_factors[0]))], wavelength, boss=boss)
    _, binned_sdss_correction, _ = generate_binned_alphas([correction_factor_sdss], [np.ones(len(correction_factor_sdss))], wavelength_sdss, boss=False)
    binned_sdss_correction = binned_sdss_correction[0]

    alphas_norad = [a / boss_fluxfactor * corr for a, corr in zip(alphas_norad, correction_factors)]
    # bin
    alphas_norad_bin_wav, alphas_norad_bin, _ = generate_binned_alphas(alphas_norad, [np.ones(alphas_norad[0].shape) for i in range(len(alphas_norad))],
                                                                       wavelength, bin_offset=0)
    alphas_boss_bin = alphas_norad_bin[0]
    alphas_north_bin = alphas_norad_bin[1]
    alphas_south_bin = alphas_norad_bin[2]
    # bin sdss also
    alphas_sdss = [np.load('../alphas_and_stds/alphas_sdss_iris_2d_102019.npy') * correction_factor_sdss / sdss_fluxfactor]  # just for the shape
    sdss_bin_wav, sdss_bin, _ = generate_binned_alphas(alphas_sdss, [np.ones(alphas_sdss[0].shape)], wavelength_sdss, bin_offset=0)
    sdss_bin = sdss_bin[0]
    if boss:
        plt.plot(alphas_norad_bin_wav, alphas_boss_bin, 'orange', label='boss alphas (observed)', drawstyle='steps-mid')

    if not boss:
        _, brandt_alphas_bin, _ = generate_binned_alphas([brandt_alphas], [np.ones(alphas.shape)], wavelength, bin_offset=0)
        brandt_alphas_bin = brandt_alphas_bin[0]
        plt.plot(lambdas_bin, brandt_alphas_bin, 'cyan', label='brandt code', drawstyle='steps-mid')

        ratio_brandt_mine = brandt_alphas / (alphas * bd12_factor)
        ratio_bd12_mine = bd12_alphas[:, 1][:-1] / (alphas_bin * bd12_factor)

        plt.plot(wavelength, ratio_brandt_mine, 'purple', label='ratbinio of brandt alphas to mine')
        plt.plot(wavelength, np.ones(len(ratio_brandt_mine)), 'purple', linestyle='--')
        plt.plot(lambdas_bin, ratio_bd12_mine, 'pink', label='ratio of bd12 alphas to mine')

    plt.plot(bd12_alphas[:, 0], bd12_alphas[:, 1], 'green', label='bd12 alphas', drawstyle='steps-post')
    # plt.plot(wavelength, wavelength * bc03_f(wavelength) / avg_bc03_wav / 2, 'b', label='~ wav*bc03')
    plt.plot(lambdas_bin, alphas_bin * bd12_factor, 'k', label='~ alphas (radiative)', drawstyle='steps-mid')
    plt.plot(wavelength, alphas * bd12_factor / wavelength / bc03_f(wavelength) * avg_bc03 * 5000, 'r', label='~ alphas / bc03 / wav')

    plt.xlim(3600, 10400)  # later: expand to BOSS wavelength range
    plt.ylim(0, .34)
    if wd01_model:
        plt.title("WD01 Model")
    else:
        plt.title("ZD04 Model")
    plt.legend(frameon=False)
    plt.show()

    # save the result
    date = ''
    save_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
    save_path += p.split('/')[-1].rsplit('.spec')[0]
    if wd01_model:
         save_path += 'wd' + date + '.npy'
    else:
         save_path += 'zd' + date + '.npy'
    print(save_path)
    np.save(save_path, alphas)

    ##########################
    # PLOT BOSS WITH MULTIPLE DUST MODELS.
    # SEE WHAT SCALING LOOKS GOOD.
    ##########################
    if boss:

        show_full = False
        cst_coeff = 0.45

        load_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
        wd_t5e9_filename = 't5e9_12gyr_z02wd_070921.npy'
        zd_t5e9_filename = 't5e9_12gyr_z02zd_070921.npy'
        wd_t9e9_filename = 't9e9_12gyr_z02wd082221.npy'
        zd_t9e9_filename = 't9e9_12gyr_z02zd082221.npy'
        wd_cst_filename = 'cst_6gyr_z02wd082221.npy'
        zd_cst_filename = 'cst_6gyr_z02zd082221.npy'
        wd01_t5e9_alphas = np.load(load_path + wd_t5e9_filename)
        zd01_t5e9_alphas = np.load(load_path + zd_t5e9_filename)
        wd01_t9e9_alphas = np.load(load_path + wd_t9e9_filename)
        zd01_t9e9_alphas = np.load(load_path + zd_t9e9_filename)
        wd01_cst_alphas = np.load(load_path + wd_cst_filename)
        zd01_cst_alphas = np.load(load_path + zd_cst_filename)
        wd01_combo = np.load(load_path + wd_t5e9_filename) + cst_coeff * np.load(load_path + wd_cst_filename)
        zd01_combo = np.load(load_path + zd_t5e9_filename) + cst_coeff * np.load(load_path + zd_cst_filename)
        # bin the alphas
        _, wd01s_bin, _ = generate_binned_alphas([wd01_t5e9_alphas, wd01_t9e9_alphas, wd01_cst_alphas, wd01_combo],
                                                 [np.ones(wd01_t5e9_alphas.shape) for i in range(4)], wavelength,
                                                 bin_offset=0)
        _, zd01s_bin, _ = generate_binned_alphas([zd01_t5e9_alphas, zd01_t9e9_alphas, zd01_cst_alphas, zd01_combo],
                                                 [np.ones(wd01_t5e9_alphas.shape) for i in range(4)], wavelength,
                                                 bin_offset=0)

        # load bootstrap errors and bin them
        loadkey = '012720'
        with open('../alphas_and_stds/bootstrap_binned_stds_boss' + loadkey + '.p', 'rb') as pf:
            bootstrap_binned_stds = pickle.load(pf)
        with open('../alphas_and_stds/bootstrap_binned_upper_boss' + loadkey + '.p', 'rb') as pf:
            bootstrap_binned_upper = pickle.load(pf)
        with open('../alphas_and_stds/bootstrap_binned_lower_boss' + loadkey + '.p', 'rb') as pf:
            bootstrap_binned_lower = pickle.load(pf)
        with open('../alphas_and_stds/bootstrap_binned_stds_sdss'+ loadkey + '.p', 'rb') as pf:
            bootstrap_binned_stds_sdss = pickle.load(pf)
        bootstrap_binned_stds_north = bootstrap_binned_stds[0] / boss_fluxfactor * binned_corrections[1]
        bootstrap_binned_upper_north = bootstrap_binned_upper[0] / boss_fluxfactor * binned_corrections[1]
        bootstrap_binned_lower_north = bootstrap_binned_lower[0] / boss_fluxfactor * binned_corrections[1]
        bootstrap_binned_stds_south = bootstrap_binned_stds[1] / boss_fluxfactor * binned_corrections[2]
        bootstrap_binned_upper_south = bootstrap_binned_upper[1] / boss_fluxfactor * binned_corrections[2]
        bootstrap_binned_lower_south = bootstrap_binned_lower[1] / boss_fluxfactor * binned_corrections[2]
        bootstrap_binned_stds_boss = bootstrap_binned_stds[2] / boss_fluxfactor * binned_corrections[0]
        bootstrap_binned_upper_boss = bootstrap_binned_upper[2] / boss_fluxfactor * binned_corrections[0]
        bootstrap_binned_lower_boss = bootstrap_binned_lower[2] / boss_fluxfactor * binned_corrections[0]
        bootstrap_binned_stds_sdss = bootstrap_binned_stds_sdss[2] / boss_fluxfactor * binned_sdss_correction

        # calculate scaling factor: (4200 to 5000 for now)
        # and just avg value for now
        # and use the north spectrum
        mean_wds = [np.mean(wd01_bin[(lambdas_bin > 4200) & (lambdas_bin < 5000)]) for wd01_bin in wd01s_bin]
        mean_zds = [np.mean(zd01_bin[(lambdas_bin > 4200) & (lambdas_bin < 5000)]) for zd01_bin in zd01s_bin]
        mean_boss = np.mean(
            alphas_boss_bin[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)])
        mean_north = np.mean(
            alphas_north_bin[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)])
        mean_south = np.mean(
            alphas_south_bin[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)])
        scale_factors_wd_boss = [round(mean_boss / mean_wd, 2) for mean_wd in mean_wds]
        scale_factors_zd_boss = [round(mean_boss / mean_zd, 2) for mean_zd in mean_zds]
        scale_factors_wd_north = [round(mean_north / mean_wd, 2) for mean_wd in mean_wds]
        scale_factors_zd_north = [round(mean_north / mean_zd, 2) for mean_zd in mean_zds]
        scale_factors_wd_south = [round(mean_south / mean_wd, 2) for mean_wd in mean_wds]
        scale_factors_zd_south = [round(mean_south / mean_zd, 2) for mean_zd in mean_zds]

        # 2-panel plot
        plt.figure(figsize=(12, 5))
        colors = ['#4477AA', '#CCBB44', '#66CCEE', '#EE6677', '#228833']

        ax1 = plt.subplot(1, 2, 1)
        ax1.set_title('BOSS alphas and BC03 models (WD dust model)')
        # plt.plot(alphas_norad_bin_wav, alphas_boss_bin, 'orange', label='BOSS alphas (observed)', drawstyle='steps-mid')
        ax1.plot(alphas_norad_bin_wav, alphas_north_bin, colors[0], label='BOSS north',
                 drawstyle='steps-mid')
        ax1.plot(alphas_norad_bin_wav, alphas_south_bin, colors[3], label='BOSS south',
                 drawstyle='steps-mid')
        if show_full:
            ax1.plot(alphas_norad_bin_wav, alphas_boss_bin, 'k', label='BOSS overall',
                     drawstyle='steps-mid')

        """
        scale_factors_wd = scale_factors_wd_south
        scale_factors_zd = scale_factors_zd_south
        ax1.plot(lambdas_bin, wd01s_bin[0] * scale_factors_wd[0], 'k',
                 label='WD01 t5e9 (x ' + str(scale_factors_wd[0]) + ')', drawstyle='steps-mid')
        ax1.plot(lambdas_bin, zd01s_bin[0] * scale_factors_zd[0], 'blue',
                 label='ZDA01 t5e9 (x ' + str(scale_factors_zd[0]) + ')', drawstyle='steps-mid')
        ax1.plot(lambdas_bin, wd01s_bin[1] * scale_factors_wd[1], 'orange',
                 label='WD01 t9e9 (x ' + str(scale_factors_wd[1]) + ')', drawstyle='steps-mid')
        ax1.plot(lambdas_bin, zd01s_bin[1] * scale_factors_zd[1], 'yellow',
                 label='ZDA01 t9e9 (x ' + str(scale_factors_zd[1]) + ')', drawstyle='steps-mid')
        ax1.plot(lambdas_bin, wd01s_bin[2] * scale_factors_wd[2], 'pink',
                 label='WD01 cst (x ' + str(scale_factors_wd[2]) + ')', drawstyle='steps-mid')
        ax1.plot(lambdas_bin, zd01s_bin[2] * scale_factors_zd[2], 'cyan',
                 label='ZDA01 cst (x ' + str(scale_factors_zd[2]) + ')', drawstyle='steps-mid')
        ax1.plot(lambdas_bin, wd01s_bin[3] * scale_factors_wd[3], 'grey',
                 label='WD01 combo (x ' + str(scale_factors_wd[3]) + ')', drawstyle='steps-mid')
        ax1.plot(lambdas_bin, zd01s_bin[3] * scale_factors_zd[3], 'brown',
                 label='ZDA01 combo (x ' + str(scale_factors_zd[3]) + ')', drawstyle='steps-mid')
        """

        # plot some models, scaled to north and south:
        ax1.plot(lambdas_bin, wd01s_bin[0] * scale_factors_wd_north[0], colors[2],
                 label='WD01 t5e9 (x ' + str(scale_factors_wd_north[0]) + ')', drawstyle='steps-mid')
        ax1.plot(lambdas_bin, wd01s_bin[0] * scale_factors_wd_south[0], colors[1],
                 label='WD01 t5e9 (x ' + str(scale_factors_wd_south[0]) + ')', drawstyle='steps-mid')
        ax1.plot(lambdas_bin, wd01s_bin[3] * scale_factors_wd_south[3], colors[4],
                 label='WD01 combo (x ' + str(scale_factors_wd_south[3]) + ')', drawstyle='steps-mid')

        ax1.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_north,
                         bootstrap_binned_upper_north, linewidth=0.0, color=colors[0], alpha=0.2,
                         step='mid')
        ax1.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_south,
                         bootstrap_binned_upper_south, linewidth=0.0, color=colors[3], alpha=0.2,
                         step='mid')
        if show_full:
            ax1.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_boss,
                             bootstrap_binned_upper_boss, linewidth=0.0, color='k', alpha=0.2,
                             step='mid')

        ax1.set_ylim(0, 0.3)
        ax1.set_xlim(3650, 10200)
        ax1.legend(frameon=False)
        ax1.set_ylabel(r"$\alpha_\lambda \beta_\lambda$ = $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mu$m)")
        ax1.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")

        ax2 = plt.subplot(1, 2, 2)
        ax2.set_title('ERE Peak')

        # define functions to get everything to line up
        zd_fns = [interp1d(lambdas_bin, zd01_bin) for zd01_bin in zd01s_bin]
        wd_fns = [interp1d(lambdas_bin, wd01_bin) for wd01_bin in wd01s_bin]
        north_fn = interp1d(alphas_norad_bin_wav, alphas_north_bin)
        south_fn = interp1d(alphas_norad_bin_wav, alphas_south_bin)
        boss_fn = interp1d(alphas_norad_bin_wav, alphas_boss_bin)
        sdss_fn = interp1d(sdss_bin_wav, sdss_bin)
        bootstrap_binned_upper_north_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_upper_north)
        bootstrap_binned_lower_north_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_lower_north)
        bootstrap_binned_upper_south_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_upper_south)
        bootstrap_binned_lower_south_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_lower_south)
        wav_clipped = lambdas_bin[(alphas_norad_bin_wav > 4000) & (alphas_norad_bin_wav < 9000)]

        # also plot south - north
        # combine errors:
        south_errors_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_stds_south)
        north_errors_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_stds_north)
        boss_errors_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_stds_boss)
        sdss_errors_fn = interp1d(sdss_bin_wav, bootstrap_binned_stds_sdss)
        bootstrap_binned_stds_sdss_fn = interp1d(sdss_bin_wav, bootstrap_binned_stds_sdss)
        def north_south_errors_fn(args):
            return np.sqrt(south_errors_fn(args) ** 2 + north_errors_fn(args) ** 2)
        def north_south_diff_fn(args):
            return south_fn(args) - north_fn(args)
        # ax2.plot(wav_clipped, north_south_diff_fn(wav_clipped), 'grey', label='South - North', drawstyle='steps-mid')
        # plt.plot(wav_clipped, north_south_errors_fn(wav_clipped), label='Bootstrap errors', drawstyle='steps-mid')

        # find median signal to noise for boss vs. sdss
        boss_signal_to_noise =  boss_fn(wav_clipped) / boss_errors_fn(wav_clipped)
        sdss_signal_to_noise = sdss_fn(wav_clipped) / sdss_errors_fn(wav_clipped)
        print("median signal to noise:")
        print("BOSS:", np.median(boss_signal_to_noise))
        print("SDSS:", np.median(sdss_signal_to_noise))

        # TEST: was using ZD dust model for the last one

        ax2.plot(wav_clipped, north_fn(wav_clipped) - wd_fns[0](wav_clipped) * scale_factors_wd_north[0], colors[2],
                 label='north - wd (t5e9)', drawstyle='steps-mid')
        ax2.plot(wav_clipped, south_fn(wav_clipped) - wd_fns[3](wav_clipped) * scale_factors_wd_south[3], colors[4], label='south - wd (combo)',
                 drawstyle='steps-mid')
        ax2.plot(wav_clipped, south_fn(wav_clipped) - wd_fns[0](wav_clipped) * scale_factors_wd_south[0], colors[1],
                 label='south - wd (t5e9)', drawstyle='steps-mid')

        ax2.fill_between(wav_clipped,
                         bootstrap_binned_lower_north_fn(wav_clipped) - wd_fns[0](wav_clipped) * scale_factors_wd_north[0],
                         bootstrap_binned_upper_north_fn(wav_clipped) - wd_fns[0](wav_clipped) * scale_factors_wd_north[0],
                         linewidth=0.0, color=colors[2], alpha=0.2, step='mid')
        ax2.fill_between(wav_clipped,
                         bootstrap_binned_lower_south_fn(wav_clipped) - wd_fns[3](wav_clipped) * scale_factors_wd_south[3],
                         bootstrap_binned_upper_south_fn(wav_clipped) - wd_fns[3](wav_clipped) * scale_factors_wd_south[3],
                         linewidth=0.0, color=colors[4], alpha=0.2, step='mid')
        ax2.fill_between(wav_clipped,
                         bootstrap_binned_lower_south_fn(wav_clipped) - wd_fns[0](wav_clipped) * scale_factors_wd_south[0],
                         bootstrap_binned_upper_south_fn(wav_clipped) - wd_fns[0](wav_clipped) * scale_factors_wd_south[0],
                         linewidth=0.0, color=colors[1], alpha=0.2, step='mid')

        ax2.legend(frameon=False)
        ax2.hlines(0, 3650, 10200, color='grey')
        ax2.set_ylim(-0.05, 0.1)
        ax2.set_xlim(4000, 9000)
        ax2.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")

        plt.tight_layout()
        if save:
            plt.savefig('/Users/blakechellew/Documents/DustProject/paper_figures/ere_2plot_091621.pdf')

        plt.show()

        # some ERE calculations:

        # weighted i100 should be an average:
        i100_weighted = np.load('/Users/blakechellew/Documents/DustProject/alphas_and_stds/avg_i100_boss_iris_smooth.npy')[0]
        i100_weighted_north = np.load('/Users/blakechellew/Documents/DustProject/alphas_and_stds/avg_i100_boss_iris_north_smooth.npy')[0]
        i100_weighted_south = np.load('/Users/blakechellew/Documents/DustProject/alphas_and_stds/avg_i100_boss_iris_south_smooth.npy')[0]

        # integrate south - north:
        integrand_south_minus_north = north_south_diff_fn(wav_clipped)
        integrand_south_minus_north /= wav_clipped
        integral_south_minus_north = 50 * np.sum(integrand_south_minus_north)
        integral_south_minus_north_err = np.sqrt(np.sum(north_south_errors_fn(wav_clipped) ** 2))
        print("Integral of south - north:", integral_south_minus_north, "+/-", integral_south_minus_north_err)

        # integrate north - model and south - model (assuming no model errors)
        integrand_north_minus_model = north_fn(wav_clipped) - wd_fns[0](wav_clipped) * scale_factors_wd_north[0]
        integrand_north_minus_model /= wav_clipped  # to get a unitless integral
        integral_north_minus_model = 50 * np.sum(integrand_north_minus_model * i100_weighted_north * (3 * 10**12))  # multiply by bin width 50 A
        integral_north_minus_model *= 10**-17  # unit conversion (to erg / s / cm^2 / sr)
        integrand_south_minus_combo = south_fn(wav_clipped) - wd_fns[3](wav_clipped) * scale_factors_wd_south[3]
        integrand_south_minus_combo /= wav_clipped
        integral_south_minus_combo = 50 * np.sum(integrand_south_minus_combo * i100_weighted_south * (3 * 10**12))
        integral_south_minus_combo *= 10**-17  # unit conversion (to erg / s / cm^2 / sr)
        integrand_south_minus_model = south_fn(wav_clipped) - wd_fns[0](wav_clipped) * scale_factors_wd_south[0]
        integrand_south_minus_model /= wav_clipped
        integral_south_minus_model = 50 * np.sum(integrand_south_minus_model * i100_weighted_south * (3 * 10**12))
        integral_south_minus_model *= 10**-17  # unit conversion (to erg / s / cm^2 / sr)
        # multiplying by (nu*I_nu)_100 should give and actual value
        # units here are MJy / sr (bc i100 is MJy / sr, and the dlambda units cancel with division by wav)
        # convert:
        # start: (MJy / sr) * Hz
        # conversion factor: x10^-20 to get erg / s / cm^2 / sr

        integral_south_minus_model_err = np.sqrt(np.sum(south_errors_fn(wav_clipped) ** 2))
        integral_north_minus_model_err = np.sqrt(np.sum(north_errors_fn(wav_clipped) ** 2))
        print("Integral of south - model:", integral_south_minus_model, "+/-", integral_south_minus_model_err)
        print("Integral of north - model:", integral_north_minus_model, "+/-", integral_north_minus_model_err)
        print("Integral of south - combo:", integral_south_minus_combo)

        # total flux:
        total_integrand_south = south_fn(wav_clipped)
        total_integrand_south /= wav_clipped
        total_flux_south = np.sum(50 * total_integrand_south * i100_weighted_south * (3 * 10**12))
        total_flux_south *= 10**-17  # unit conversion
        total_integrand_north = north_fn(wav_clipped)
        total_integrand_north /= wav_clipped
        total_flux_north = np.sum(50 * total_integrand_north * i100_weighted_north * (3 * 10**12))
        total_flux_north *= 10 ** -17  # unit conversion
        print("Total flux south:", total_flux_south)
        print("Total flux north:", total_flux_north)

        # flux ratio:
        south_ERE_ratio = integral_south_minus_model / (total_flux_south - integral_south_minus_model)
        north_ERE_ratio = integral_north_minus_model / (total_flux_north - integral_north_minus_model)
        south_ERE_ratio_combo = integral_south_minus_combo / (total_flux_south - integral_south_minus_combo)
        print("ratio ERE to scattered, south:", south_ERE_ratio)
        print("ratio ERE to scattered, north:", north_ERE_ratio)
        print("ratio ERE to scattered, south combo:", south_ERE_ratio_combo)

        # comparing spectra (z-values):
        avg_north = np.mean(north_fn(wav_clipped))
        avg_south = np.mean(south_fn(wav_clipped))
        avg_north_err = np.median(north_errors_fn(wav_clipped))
        avg_south_err = np.median(south_errors_fn(wav_clipped))
        north_vs_south_err = np.sqrt(avg_north_err**2 + avg_south_err**2)
        z_north_vs_south = (avg_south - avg_north) / north_vs_south_err
        print("z north vs south:", z_north_vs_south)
        avg_sdss = np.mean(sdss_fn(wav_clipped))
        avg_boss = np.mean(boss_fn(wav_clipped))
        avg_sdss_err = np.median(bootstrap_binned_stds_sdss_fn(wav_clipped))
        sdss_vs_north_err = np.sqrt(avg_sdss_err**2 + avg_north_err**2)
        sdss_vs_south_err = np.sqrt(avg_sdss_err**2 + avg_south_err**2)
        z_sdss_vs_north = (avg_north - avg_sdss) / sdss_vs_north_err
        z_sdss_vs_south = (avg_south - avg_sdss) / sdss_vs_south_err
        print("z SDSS vs north:", z_sdss_vs_north)
        print("z SDSS vs south:", z_sdss_vs_south)
        print("avgs:")
        print("north:", avg_north)
        print("south:", avg_south)
        print("sdss:", avg_sdss)
        print("boss:", avg_boss)

        # find the peaks and quartiles (don't care about the scaling)
        # first the edges
        def find_quartiles(integrand):
            wavs_subzero = wav_clipped[np.nonzero(integrand < 0)]
            lower_edge = 6500 - np.min(np.abs(wavs_subzero[wavs_subzero < 6500] - 6500))
            upper_edge = 6500 + np.min(np.abs(wavs_subzero[wavs_subzero > 6500] - 6500))
            print("lower edge:", lower_edge)
            print("upper edge:", upper_edge)
            peak_integrand = integrand_south_minus_model[(wav_clipped > lower_edge) & (wav_clipped < upper_edge)]
            peak_wavs = wav_clipped[(wav_clipped > lower_edge) & (wav_clipped < upper_edge)]
            total_integrand = np.sum(peak_integrand)
            partial_sum = 0
            quartile_1 = None
            quartile_2 = None
            quartile_3 = None
            for w, e in zip(peak_wavs, peak_integrand):
                partial_sum += e
                if quartile_1 is None and partial_sum > total_integrand / 4:
                    quartile_1 = w
                if quartile_2 is None and partial_sum > 2 * total_integrand / 4:
                    quartile_2 = w
                if quartile_3 is None and partial_sum > 3 * total_integrand / 4:
                    quartile_3 = w
                    break  # no point continuing
            print("quartiles 1, 2, 3:", quartile_1, quartile_2, quartile_3)
            print("width:", quartile_3 - quartile_1)
        find_quartiles(integrand_south_minus_model)
        find_quartiles(integrand_south_minus_combo)

#############################

"""
# test: try a bunch of different binnings
# currently getting an error because shifting the binning changes the total number of elements
# I would expect the ratio to be smooth...
for b in range(60):
    lambdas_bin, alphas_bin, _ = generate_binned_alphas([alphas], [np.ones(alphas.shape)], wavelength,
                                                        bin_offset=b)
    alphas_bin = alphas_bin[0]
    plt.plot(lambdas_bin, alphas_bin * bd12_factor, 'k', label='~ alphas (radiative)', drawstyle='steps')
    plt.plot(bd12_alphas[:, 0], bd12_alphas[:, 1], 'green', label='bd12 wd01', drawstyle='steps-mid')

    ratio1 = bd12_alphas[:, 1][:100] / alphas_bin[:100]
    ratio2 = bd12_alphas[:, 1][1:101] / alphas_bin[:100]
    ratio3 = bd12_alphas[:, 1][2:102] / alphas_bin[:100]
    plt.plot(lambdas_bin[:100], ratio1, 'pink', label='ratio1')
    plt.plot(lambdas_bin[:100], ratio2, 'purple', label='ratio2')
    plt.plot(lambdas_bin[:100], ratio3, 'blue', label='ratio3')

    plt.title('Binning Offset: ' + str(b))
    plt.show()
"""

"""
# run for different longitude / latitude ranges
elif location == 3:
	ivar[b < 50] = 0 # above 50 degrees
elif location == 4:
	ivar[b < 50] = 0
	ivar[l > 120] = 0
elif location == 5:
	ivar[b < 50] = 0
	ivar[l < 120] = 0
	ivar[l > 240] = 0
elif location == 6:
	ivar[b < 50] = 0
	ivar[l < 240] = 0
"""