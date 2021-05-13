# TODO
# check units again, for ex. cross section is cm^2
# still need to verify if the masking for numerically bad functions is working

# more things to verify:
# double check that henyey fn is used correctly

# verify density prefactor calculation

# verifying the integral:
# it should match the results from BC03
# make sure it converges; should stay same if I increase the bounds, and number of divs
# try to guestimate the value based on integrand plots?

# take the BC03 model spectra and apply radiative transfer calculations.
# based on the appendix of BD12.

import glob
import numpy as np
from scipy.interpolate import interp1d
import scipy.special as special
import matplotlib.pyplot as plt
from generate_plots import generate_binned_alphas

# load some files
wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')  # angstroms
paths = glob.glob('/Users/blakechellew/Documents/DustProject/BrandtFiles/bc03/*.spec')

# load dust models
dust_model_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/kext_albedo_WD_MW_3.1_60_D03.all'
# column 0: wavelength (microns, but converted to angstroms below)
# column 1: albedo
# column 2: <cos>
# column 3: cross section (cm^2)
# column 5: <cos^2>

dust_model = np.loadtxt(dust_model_path, skiprows=80, usecols=(0, 1, 2, 3, 4, 5))

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

# choose a latitude
b = 40 * np.pi / 180  # 40 degrees -> radians


# phase function
def henyey(cos_xi, g):
    g = g * (-1)  # TEST: not sure why * -1
    result = (1 - g ** 2) / (1 + g ** 2 - 2 * g * cos_xi) ** (3 / 2) / (4 * np.pi)
    assert result.shape == cos_xi.shape
    return result


sig_dust = 1  # parsecs (from density eqn)   # TEST: change to 1 pc (should be 250)
tau_def = 0.15  # fiducial value (see BD12 paper) (z = 0)
V_band_wav = 5510  # A
cross_sec_V = dust_cross_f(V_band_wav)  # cm^2

# eqn A1
# z should be in pc
# density: rho ~ exp(-z^2 / 2 / sig_dust^2).
# I took sigma out of the integral and integrated analytically.
# units of cross section should cancel out b/c it's also in the density prefactor.
def tau_f(lamb, z):
    cross_section = dust_cross_f(lamb)  # cm^2
    tau_0 = tau_def * cross_section / cross_sec_V
    tau_lambda = tau_0 * (1 + special.erf(z / sig_dust / np.sqrt(2))) # TEST: was 1 - erf
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
    return np.abs(result)  #TEST: apply abs value

# scale heights 300 pc and 1350 pc
a_300 = 0.9
a_1350 = 0.1


# surface power density (arbitrary units; any normalization factor will cancel out later)
# it's a function of height z_s
# the bc03 flux will be multiplied later
def surface_power_fn(z, rho, beta):  # z_s in pc
    return a_300 * np.exp(-np.abs(z - rho * np.cos(beta)) / 300) \
        + a_1350 * np.exp(-np.abs(z - rho * np.cos(beta)) / 1350)


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
    result = result * 1000 * np.exp(rho/1000)  # transform rho to u
    return result

# eqn A4: integrate over everything except wavelength
def i_tir(bc03_f):
    n_lamb = 200
    n_u = n_tau = n_beta = 30
    # n_u = 30
    # n_tau = 30
    # n_beta = 30  # should be even to avoid beta=pi

    lamb_min = 100
    lamb_max = 15*10**4   # 10 ** 6
    h_lamb = (lamb_max - lamb_min) / 2 / n_lamb
    lamb_grid = np.linspace(lamb_min + h_lamb, lamb_max - h_lamb, n_lamb)  # 91 to 10**4 A  # min, max, n
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
    rhos = -1000*np.log(us)
    assert lambs.shape == (n_tau, n_lamb, n_beta, n_u)

    ww = i_tir_integrand(lambs, z_of_tau(lambs, taus), rhos, betas)

    # integrate on simple grid
    lamb_div = (lamb_max - lamb_min) / n_lamb
    beta_div = (beta_max - beta_min) / n_beta
    u_div = (u_max - u_min) / n_u
    tau_div = (tau_max - tau_min) / n_tau
    tau_div = tau_div[..., None, None] + zeros

    # multiply by bc03
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
    result = result * 1000 * np.exp(rho / 1000)  # transform rho to u

    return result


# eqn A7, part 2: do the integral for given wavelength
# still needs to be multiplied by bc03 flux
def i_sca(lamb):

    # n_theta = 20
    # n_u = 20
    # n_tau = 20
    # n_beta = 20  # should be even to avoid beta=pi
    n_theta = n_u = n_tau = n_beta = 30
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
    h_theta_temp = (theta_max - theta_min) / 2 / n_plot_pts
    h_tau_temp = (tau_max - tau_min) / 2 / n_plot_pts
    h_u_temp = (u_max - u_min) / 2 / n_plot_pts
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
    # plot vs. u
    us_temp = np.linspace(u_min + h_u_temp, u_max - h_u_temp, n_plot_pts)
    integrand_vals = i_sca_integrand(theta_temp, tau_temp, -1000 * np.log(us_temp), beta_temp, lamb, bc03_f)
    print("us", integrand_vals)
    assert np.count_nonzero(np.isnan(integrand_vals)) == 0
    plt.plot(us_temp, integrand_vals, '.')
    plt.title('I_sca integrand vs. e^(-rho / 1000 pc)')
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
    rhos = -1000*np.log(us)
    ww = i_sca_integrand(thetas, taus, rhos, betas, lamb)

    # sum over grid
    theta_div = (theta_max - theta_min) / n_theta
    u_div = (u_max - u_min) / n_u
    tau_div = (tau_max - tau_min) / n_tau
    beta_div = (beta_max - beta_min) / n_beta
    result = np.sum(ww) * theta_div * u_div * tau_div * beta_div
    return result


# idx 21 and 22: exponential, 12 Gyr, z=.02 (.02 is correct)
# index 2: t=5e9, 12 Gyr, Z=.008
# 5th to last: similar but t9e9

# loop through the bc03 spectra
for p in paths[22:23]:
    print("verify path name:", p)

    # load and interpolate the bc03 spectra
    a = np.loadtxt(p)
    wav = a[:, 0]  # angstroms
    bc03 = a[:, 1]
    bc03_f = interp1d(wav, bc03, kind='cubic')

    """
    plt.plot(wav, bc03, label='bc03')
    plt.plot(wav, bc03 * wav, label='* wav')
    plt.plot(wav, bc03 / wav, label='/ wav')
    plt.semilogx()
    plt.legend()
    plt.xlim(800, 13000)
    plt.show()
    """

    # plot previously saved stuff:
    """
    load_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
    load_path += p.split('/')[-1].rsplit('.spec')[0] + '_021121.npy'  # 021121 or 022321
    alphas_radiative = np.load(load_path)
    bc03_factor = 3  # test
    bd12_factor = 0.49
    plt.plot(wavelength, alphas_radiative * bd12_factor, label='Radiative Transfer', drawstyle='steps')
    plt.plot(wavelength, bc03_f(wavelength) / bc03_factor * bd12_factor, label='BC03 Model', drawstyle='steps')
    plt.plot(wavelength, alphas_radiative / bc03_f(wavelength), drawstyle='steps')
    plt.xlim(3000, 10000)
    plt.ylim(0, 0.4)
    plt.xlabel("Wavelength")
    plt.ylabel("Alpha")
    plt.legend()
    plt.show()
    """

    # testing A4:
    # i_tir = i_tir_integrand(6000, 2, bc03_f)
    # print("I tir test")
    # print(i_tir)
    # exit(0)

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
    # (0 to inf throws an error)
    # (range of lambda for dust model is 10^-4 to 10^4 microns, or 1 to 10^8 angstroms)
    # (note: range of lambda for bc03 is 91 to 1.6e6)

    I_tir = i_tir(bc03_f)
    print("I_tir result")
    print(I_tir)

    # convert to 100 micron (there is an associated uncertainty)
    nu_I_nu_100 = .52 * I_tir  # units of angstroms * sigma
    assert nu_I_nu_100.ndim == 0

    print("starting integral")

    # wavelength_partial = wavelength[(wavelength > 4600) & (wavelength < 5600)]  #5380 to 5480
    wavelength_partial = wavelength[::40]
    # wavelength_partial = np.array([4996])
    # wavelength_partial = [3842, 3920, 4450, 4996, 6480, 6660, 8875, 9553]
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

    # bin the alphas
    lambdas_boss_bin, alphas_bin, _ = generate_binned_alphas([alphas], [np.ones(alphas.shape)], wavelength)
    alphas_bin = alphas_bin[0]

    # scaling factors
    bd12_factor = 0.49
    avg_bc03 = np.mean(bc03_f(wavelength_partial))
    avg_bc03_wav = np.mean(bc03_f(wavelength_partial) * wavelength_partial)

    alphas_norad = np.load('../alphas_and_stds/alphas_boss_iris_2d_012720.npy')

    plt.plot(wavelength, alphas_norad * 3 / 2, 'g', label='alphas (no radiative)')
    plt.plot(wavelength, wavelength * bc03_f(wavelength) / avg_bc03_wav / 2, 'b', label='~ wav*bc03')
    plt.plot(lambdas_boss_bin, alphas_bin * bd12_factor, 'k', label='~ alphas (radiative)', drawstyle='steps')
    plt.plot(wavelength, alphas * bd12_factor / wavelength / bc03_f(wavelength) * avg_bc03 * 5000, 'r', label='~ alphas / bc03 / wav')

    plt.xlim(3800, 9200)
    plt.ylim(0, 1)
    plt.legend()
    plt.show()

    bc03_test = wav * bc03_f(wav)
    print("bc03:", bc03_test[:100])

    # save the result
    save_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
    save_path += p.split('/')[-1].rsplit('.spec')[0] + '.npy'
    print(save_path)
    np.save(save_path, alphas)
