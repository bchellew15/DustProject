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
import sys

#command line args
#see above for explanation
if len(sys.argv) != 4:
    print("Usage: radiative_transfer_v2.py [p_num: 34, etc.] [wd01_model: 0, 1] [save: savekey]")
    exit(0)
p_num = int(sys.argv[1])
wd01_model = int(sys.argv[2])
savekey = sys.argv[3]  # t5e9_grid_200_wd.npy


save = False  # save ERE figure
boss = True  # whether to interpolate to BOSS wavelengths or SDSS
# wd01_model = True  # otherwise zda04 model
b = 40 * np.pi / 180  # latitude (should be 40 degrees -> radians)
sig_dust = 250  # 250  # 250 parsecs (from density eqn)   # (alt: 1)
tau_def = 0.15  # fiducial value 0.15 (see BD12 paper) (z = 0)
V_band_wav = 5510  # A
# stellar scale heights 300 pc and 1350 pc
a_300 = .9  # 0.9
a_1350 = .1  # 0.1
#which bc03 model to use:
sig_star_1 = 300  # 300 (alt 1.2)
sig_star_2 = 1350  # 1350 (alt 5.4)
mistakes = False
rho_scale = 500  # 500  # previous: 5 or 1000  #.1 matches brandt code for sig = 1
ahab = False
uv = True
uv_cutoff = 2 # 1 is 6 eV, 2 is 8 eV

# number of grid points (for TIR and scattering):
# n_beta needs to be even
# default: all 20s
n_u = 140  #150  # 100
n_tau = 50  #40  # 25
n_beta = 140  #100  # should be even to avoid beta=pi, but hasn't been an issue since
n_theta = 50  #30  # 25
if ahab:
    paths = sorted(glob.glob('/home/blakechellew/data/bc03/*.spec'))
else:
    paths = sorted(glob.glob('/Users/blakechellew/Documents/DustProject/BrandtFiles/bc03/*.spec'))  # bc03 model spectra

# load wavelengths
wavelength_boss = np.load('../alphas_and_stds/wavelength_boss.npy')  # angstroms
if ahab:
    hdulist_sdsswav = fits.open('/home/blakechellew/data/SDSS_allskyspec.fits')
else:
    hdulist_sdsswav = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
wavelength_sdss = np.array(hdulist_sdsswav[1].data)
if boss:
    wavelength = wavelength_boss
else:
    wavelength = wavelength_sdss

# load dust models
if wd01_model:
    if ahab:
        dust_model_path = '/home/blakechellew/data/wd01_dustmodel.all'
    else:
        dust_model_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/kext_albedo_WD_MW_3.1_60_D03.all'
    skip_rows = 80
    use_cols = (0, 1, 2, 3)
    # column 0: wavelength (microns, but converted to angstroms below)
    # column 1: albedo
    # column 2: <cos>
    # column 3: cross section (cm^2)
else:  # zda04 dust model
    if ahab:
        dust_model_path = '/home/blakechellew/data/zda04_dustmodel.txt'
    else:
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

# set wavelength range for I_tir
# 100 A to 42658 A
if wd01_model:
    wav_min = 405
    wav_max = 714
else:
    wav_min = 600
    wav_max = 1126
if uv:
    # bruce suggested cutoffs of 6 eV or 8 eV.
    # corresponding to (E = hf)
    wav_min = 600
    if uv_cutoff == 1:
        wav_max = 864  # index corresponding to 2066.4 A  (actual 2065.4)
    elif uv_cutoff == 2:
        wav_max = 839  # index corr. to 1549.8 A  (actual 1548.8)
dust_wav = dust_wav[wav_min:wav_max]

# print(dust_wav[-1])
# exit(0)

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
def i_tir_integrand(lamb, z, rho, beta, transform=1):
    prefactor = 1 / (8 * np.pi) / np.sin(np.abs(b))
    term1 = 1 - dust_albedo_f(lamb)
    term2 = surface_power_fn(z, rho, beta)
    term3 = np.exp(-A_f(lamb, z, rho, beta))
    term4 = np.sin(beta)
    result = prefactor * term1 * term2 * term3 * term4
    if transform == 1:
        result = result * rho_scale  # transformation 1
    elif transform == 2:
        result = result * rho**2 / rho_scale  # transformation 2
    # result = result * rho_scale * np.exp(rho/rho_scale)  # old: transform rho to u
    return result

# eqn A4: integrate over everything except wavelength
def i_tir(bc03_f):

    print("starting I_tir", flush=True)

    lamb_grid = dust_wav  # use wavelengths from dust model.
    # not worried about intervals because cutting off the edges anyway.
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
    # rhos = -rho_scale*np.log(us)  # TEST: remove
    rhos_1 = rho_scale * us  # transformation 1
    rhos_2 = rho_scale / us  # transformation 2
    assert lambs.shape == (n_tau, n_lamb, n_beta, n_u)

    print("starting first integral", flush=True)
    ww1 = i_tir_integrand(lambs, z_of_tau(lambs, taus), rhos_1, betas, transform=1)
    print("starting second integral", flush=True)
    ww2 = i_tir_integrand(lambs, z_of_tau(lambs, taus), rhos_2, betas, transform=2)
    print("finished second integral", flush=True)

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
    assert bc03s.shape == (1, ww1.shape[1], 1, 1)
    result = np.sum(ww1 * bc03s * tau_div * beta_div * u_div * lamb_div)  \
        + np.sum(ww2 * bc03s * tau_div * beta_div * u_div * lamb_div)  # units of angstroms * sigma
    return result


# eqn A7, part 1: just the integrand.
# note that we integrate over z_s also.
# units: parsecs * units of sigma
def i_sca_integrand(theta, tau, rho, beta, lamb, transform=1):
    prefactor = (1 / np.sin(np.abs(b))) * dust_albedo_f(lamb) / (4 * np.pi)
    z = z_of_tau(lamb, tau)

    term1 = np.exp((-1 / np.sin(np.abs(b))) * (tau_f(lamb, 0) - tau))
    term2 = henyey(cos_xi(z, rho, theta, beta), dust_cos_f(lamb))
    term3 = surface_power_fn(z, rho, beta)
    term4 = np.exp(-A_f(lamb, z, rho, beta))
    term5 = np.sin(beta)
    result = prefactor * term1 * term2 * term3 * term4 * term5
    if transform == 1:
        result = result * rho_scale  # transformation 1
    elif transform == 2:
        result = result * rho**2 / rho_scale  # transformation 2
    # result = result * rho_scale * np.exp(rho / rho_scale)  # old: transform rho to u

    return result


# eqn A7, part 2: do the integral for given wavelength
# still needs to be multiplied by bc03 flux
def i_sca(lamb):

    print("scattering integral:", lamb, flush=True)  # to keep track of how long the code has left to run

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
    # rhos = -rho_scale*np.log(us)  # old
    rhos1 = rho_scale * us  # transformation 1
    rhos2 = rho_scale / us  # transformation 2
    ww1 = i_sca_integrand(thetas, taus, rhos1, betas, lamb, transform=1)
    ww2 = i_sca_integrand(thetas, taus, rhos2, betas, lamb, transform=2)

    # sum over grid
    theta_div = (theta_max - theta_min) / n_theta
    u_div = (u_max - u_min) / n_u
    tau_div = (tau_max - tau_min) / n_tau
    beta_div = (beta_max - beta_min) / n_beta
    result = (np.sum(ww1) + np.sum(ww2)) * theta_div * u_div * tau_div * beta_div
    return result

"""
for i, p in enumerate(paths):
    print(i, p)
exit(0)
"""

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

    # TEST plotting bc03
    plt.plot(wav, bc03_f(wav))
    plt.show()
    print("lyman continuum level")
    print(np.mean(bc03_f(wav)[(wav>450) & (wav<900)]))
    exit(0)

    """
    # plots of integrand for I_TIR: (lambda, z, rho, beta)
    num_div = 50
    tau_grid = np.linspace(0, tau_f(6000, 0), num_div)  # z_s grid, -inf to inf
    beta_grid = np.linspace(0, np.pi, num_div)
    u_grid = np.linspace(0, 1, 400)  # TEST
    lamb_grid = dust_wav

    # order: lambda, z, rho, beta
    lamb_plot = i_tir_integrand(lamb_grid, 100, 100, np.pi/4) * bc03_f(lamb_grid)
    plt.vlines([np.min(dust_wav), 930, np.max(dust_wav)], 0, 100)
    plt.plot(lamb_grid, lamb_plot)
    plt.title('I_TIR integrand vs. lambda')
    plt.show()
    tau_plot = i_tir_integrand(6000, z_of_tau(6000, tau_grid), 100, np.pi/4)
    print(tau_plot)
    plt.plot(tau_grid, tau_plot)
    plt.title('I_TIR integrand vs. tau')
    plt.show()
    beta_plot = i_tir_integrand(6000, 100, 100, beta_grid)
    print(beta_plot)
    plt.plot(beta_grid, beta_plot)
    plt.title('I_TIR integrand vs. beta')
    plt.show()
    # TEST: different var transformations
    # rho_grid = -rho_scale*np.log(u_grid)
    rho_grid = np.linspace(0, 500, 500)
    rho_grid_1 = 100 * u_grid  # first 100
    rho_grid_2 = 100 / u_grid

    u_plot_1 = i_tir_integrand(6000, 100, rho_grid_1, np.pi / 4)
    u_plot_2 = i_tir_integrand(6000, 100, rho_grid_2, np.pi / 4)
    # print(u_plot)
    plt.plot(u_grid, u_plot_1)  # TEST: was u_grid
    plt.plot(u_grid, u_plot_2)  # TEST
    plt.title(r'I_TIR integrand vs. u (u = $\rho / \rho_0$ or u = $\rho_0 / \rho$)')
    plt.show()
    """

    # compute total infrared radiation
    # (range of lambda for dust model is 10^-4 to 10^4 microns, or 1 to 10^8 angstroms)
    # (range of lambda for bc03 is 91 to 1.6e6)
    I_tir = i_tir(bc03_f)
    print("I_tir result")
    print(I_tir)

    exit(0)

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

    """
    # check for convergence
    i_4500 = i_sca(4500) * bc03_f(4500) * 4500
    i_7000 = i_sca(7000) * bc03_f(7000) * 7000
    i_9500 = i_sca(9500) * bc03_f(9500) * 9500
    print(str(round(i_4500, 1)) + ', ' + str(round(i_7000, 1)) + ', ' + str(round(i_9500, 1)))
    exit(0)
    """

    # try to use broadcasting for this?
    i_sca_array = np.array([i_sca(lamb) for lamb in wavelength_partial])  # units of sigma * parsecs

    # interpolate the result of the scattering integral
    i_sca_f = interp1d(wavelength_partial, i_sca_array, kind='cubic')

    # then calculate alphas
    i_lam_array = bc03_f(wavelength) * i_sca_f(wavelength) * wavelength  # units of sigma * parsecs * angstroms
    alphas = i_lam_array / nu_I_nu_100
    assert alphas.shape == wavelength.shape

    # save: test
    if ahab:
        save_path = '/home/blakechellew/alphas_and_stds/'
    else:
        save_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/convergence_test/'
    # save_path += 't5e9_grid_200_wd.npy'
    save_path += savekey
    print(save_path)
    np.save(save_path, alphas)

    if not ahab:

        # scaling factors
        bd12_factor = 0.49 if wd01_model else 0.52
        boss_fluxfactor = 1.38
        sdss_fluxfactor = 1.38

        avg_bc03 = np.mean(bc03_f(wavelength_partial))

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
        _, binned_corrections, _ = generate_binned_alphas(correction_factors, 3 * [np.ones(len(correction_factors[0]))], wavelength, boss=boss)
        alphas_norad = [a / boss_fluxfactor * corr for a, corr in zip(alphas_norad, correction_factors)]

        # bin
        alphas_norad_bin_wav, alphas_norad_bin, _ = generate_binned_alphas(alphas_norad, [np.ones(alphas_norad[0].shape) for i in range(len(alphas_norad))],
                                                                           wavelength, bin_offset=0)
        alphas_boss_bin = alphas_norad_bin[0]
        alphas_north_bin = alphas_norad_bin[1]
        alphas_south_bin = alphas_norad_bin[2]

        # plots:
        if boss:
            plt.plot(alphas_norad_bin_wav, alphas_boss_bin, 'orange', label='boss alphas (observed)', drawstyle='steps-mid')

        if not boss:
            _, brandt_alphas_bin, _ = generate_binned_alphas([brandt_alphas], [np.ones(alphas.shape)], wavelength, bin_offset=0)
            brandt_alphas_bin = brandt_alphas_bin[0]
            plt.plot(lambdas_bin, brandt_alphas_bin, 'cyan', label='brandt code', drawstyle='steps-mid')

            ratio_brandt_mine = brandt_alphas / (alphas * bd12_factor)
            ratio_bd12_mine = bd12_alphas[:, 1][:-1] / (alphas_bin * bd12_factor)

            plt.plot(wavelength, ratio_brandt_mine, 'purple', label='ratio of brandt alphas to mine')
            plt.plot(wavelength, np.ones(len(ratio_brandt_mine)), 'purple', linestyle='--')
            plt.plot(lambdas_bin, ratio_bd12_mine, 'pink', label='ratio of bd12 alphas to mine')

        plt.plot(bd12_alphas[:, 0], bd12_alphas[:, 1], 'green', label='bd12 alphas', drawstyle='steps-post')
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

        """
        # save the result
        date = ''
        save_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
        save_path += p.split('/')[-1].rsplit('.spec')[0]
        save_path += 'default_'
        if wd01_model:
             save_path += 'wd' + date + '.npy'
        else:
             save_path += 'zd' + date + '.npy'
        """

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