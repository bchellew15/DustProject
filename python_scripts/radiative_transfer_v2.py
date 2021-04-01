# TODO
# check units again, for ex. cross section is cm^2
# still need to verify if the masking for numerically bad functions is working

# more things to verify:
# double check that henyey fn is used correctly

# decide if trapezoid method is accurate enough
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
from scipy import integrate
import matplotlib.pyplot as plt

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
    return (1 - g ** 2) / (1 + g ** 2 - 2 * g * cos_xi) ** (3 / 2) / (4 * np.pi)


sig_dust = 250  # parsecs (from density eqn)
# prefactor: calculate based on tau(z=0) = 0.15. Use V band.
V_band_wav = 5510  # A
density_prefactor = 0.3 / dust_cross_f(V_band_wav) / sig_dust / np.sqrt(2 * np.pi)  # cm^-2 * pc^-1


# eqn A1
# z should be in pc
# density: rho ~ exp(-z^2 / 2 / sig_dust^2).
# I took sigma out of the integral and integrated analytically.
# units of cross section should cancel out b/c it's also in the density prefactor.
def tau_f(lamb, z):
    cross_section = dust_cross_f(lamb)  # cm^2
    tau_lambda = density_prefactor * cross_section * sig_dust * np.sqrt(2 * np.pi) / 2 \
        * (1 - special.erf(z / sig_dust / np.sqrt(2)))
    return tau_lambda  # unitless


# reverse: find z given tau
# (see paper calculations)
# WARNING: it will multiply cross section across the last axis
def z_of_tau(lamb, tau):
    cross_section = dust_cross_f(lamb)  # cm^2

    print("z of tau broadcasting")
    print(tau.shape)
    print(cross_section.shape)

    result = sig_dust * np.sqrt(2) * special.erfinv(
        1 - 2 * tau / (density_prefactor * cross_section * sig_dust * np.sqrt(2 * np.pi)))
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
    return result

# scale heights 300 pc and 1350 pc
a_300 = 0.9
a_1350 = 0.1


# surface power density (arbitrary units; any normalization factor will cancel out later)
# it's a function of height z_s
# def surface_power_fn(bc03, z_s, lamb):  # z_s in pc
#     return bc03(lamb) * (a_300 * np.exp(-np.abs(z_s) / 300) + a_1350 * np.exp(-np.abs(z_s) / 1350))


# I found the derivative analytically
def surface_power_deriv(bc03, z, rho, beta, lamb):
    return bc03(lamb) * \
           (-a_300 * np.exp(-np.abs(z - rho * np.cos(beta)) / 300) / 300
            - a_1350 * np.exp(-np.abs(z - rho * np.cos(beta)) / 1350) / 1350)


# eqn A4: integrand
# units: angstroms * units of sigma after integration
# broadcasting: lamb, z, rho, beta should all be 4d arrays
def i_tir_integrand(lamb, z, rho, beta, bc03):
    prefactor = 1 / (8 * np.pi) / np.sin(np.abs(b))
    term1 = 1 - dust_albedo_f(lamb)
    term2 = surface_power_deriv(bc03, z, rho, beta, lamb)
    term3 = np.exp(-A_f(lamb, z, rho, beta))
    term4 = np.sin(beta)
    result = prefactor * term1 * term2 * term3 * term4
    return result

# eqn A4: integrate over everything except wavelength
def i_tir():
    n_lamb = 30
    n_rho = 27
    n_tau = 28
    n_beta = 29

    lamb_min = 100
    lamb_max = 6*10**4   # 10 ** 6
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
    rho_min = 0
    rho_max = 500
    h_rho = (rho_max - rho_min) / 2 / n_rho
    rho_grid = np.linspace(rho_min + h_rho, rho_max - h_rho, n_rho)

    # meshgrid from the 3 easy dimensions (lambda, beta, rho)
    lambs, betas, rhos = np.meshgrid(lamb_grid, beta_grid, rho_grid, indexing='ij')  # so input dims match output

    # broadcasting test
    print("lambs shape", lambs.shape)
    print("betas shape", betas.shape)
    print("rhos shape", rhos.shape)
    print("tau grid", tau_grid.shape)

    # add a new dimension for tau:
    new_shape = tuple([n_tau] + list(lambs.shape))  # lambs, betas, rhos have same shape
    zeros = np.zeros(new_shape)
    lambs = lambs[None, ...] + zeros
    betas = betas[None, ...] + zeros
    rhos = rhos[None, ...] + zeros
    taus = tau_grid[..., None, None] + zeros

    ww = i_tir_integrand(lambs, z_of_tau(lambs, taus), rhos, betas, bc03_f)
    print("sum:", np.sum(ww))

    # broadcast test
    print("ww shape", ww.shape)

    # integrate on simple grid
    lamb_div = (lamb_max - lamb_min) / (n_lamb - 1)
    beta_div = (beta_max - beta_min) / (n_beta - 1)
    rho_div = (rho_max - rho_min) / (n_rho - 1)
    tau_div = (tau_max - tau_min) / (n_tau - 1)
    tau_div = tau_div[..., None, None] + zeros
    print("tau div shape:", tau_div.shape)

    result = np.sum(ww * tau_div * beta_div * rho_div * lamb_div)  # units of angstroms * sigma
    return result


# eqn A7, part 1: just the integrand.
# note that we integrate over z_s also.
# units: parsecs * units of sigma
def i_sca_integrand(theta, tau, rho, beta, lamb, bc03):
    prefactor = (1 / np.sin(np.abs(b))) * dust_albedo_f(lamb) / (4 * np.pi)
    z = z_of_tau(lamb, tau)

    # check for tau = 0
    if type(tau) != np.ndarray and tau == 0:
        return 0  # this happens when z = infinity
    z[tau == 0] = 1  # just so we don't have z = inf in cos_xi

    term1 = np.exp((-1 / np.sin(np.abs(b))) * (tau_f(lamb, 0) - tau))
    term2 = henyey(cos_xi(z, rho, theta, beta), dust_cos_f(lamb))
    term3 = surface_power_deriv(bc03, z, rho, beta, lamb)
    term4 = np.exp(-A_f(lamb, z, rho, beta))
    term5 = np.sin(beta)
    result = prefactor * term1 * term2 * term3 * term4 * term5

    if type(tau) == np.ndarray:
        np.putmask(result, tau == 0, 0)

    return result


# eqn A7, part 2: do the integral for given wavelength
def i_sca(lamb, bc03):
    print("check:", lamb)

    """
    # plot integrand for I_sca:
    # vars: theta, tau, R, zs
    # plot vs. theta (for given tau, R, zs)
    theta_vals = np.linspace(0, 2 * np.pi, 50)
    tau_temp = 0.1
    R_temp = 1
    zs_temp = 1.5
    integrand_vals = i_sca_integrand(theta_vals, tau_temp, R_temp, zs_temp, lamb, bc03_f)
    plt.plot(theta_vals, integrand_vals)
    plt.title('I_sca integrand vs. theta')
    plt.show()

    # plot vs. tau (for given theta, R, zs)
    tau_vals = np.linspace(0, tau_f(lamb, 0)-.00001, 300)
    theta_temp = np.pi
    R_temp = 1
    zs_temp = 1.5
    integrand_vals = i_sca_integrand(theta_temp, tau_vals, R_temp, zs_temp, lamb, bc03_f)
    plt.plot(tau_vals, integrand_vals)
    plt.title('I_sca integrand vs. tau')
    plt.show()

    # plot vs. z_s (for given tau, theta, R)
    zs_vals = np.linspace(-100, 1000, 50)
    tau_temp = 0.1
    theta_temp = np.pi
    R_temp = 1
    integrand_vals = i_sca_integrand(theta_temp, tau_temp, R_temp, zs_vals, lamb, bc03_f)
    plt.plot(zs_vals, integrand_vals)
    plt.title('I_sca integrand vs. z_s')
    plt.show()

    # plot vs. R (for given tau theta, zs)
    R_vals = np.linspace(.1, 1000, 50)
    tau_temp = 0.21  # temp to get z = zs: tau_f(lamb, 1.5)
    theta_temp = np.pi
    zs_temp = 1.5
    integrand_vals = i_sca_integrand(theta_temp, tau_temp, R_vals, zs_temp, lamb, bc03_f)
    plt.plot(R_vals, integrand_vals, '.')
    plt.title('I_sca integrand vs. R')
    plt.show()
    """
    """
    # plot vs. lambda (for given tau, theta, R, zs)
    R_temp = 1
    lamb_vals = np.linspace(4000, 10000, 1200)
    tau_temp = tau_f(lamb_vals, 1.5)  # temp to get z = zs; was 0.1
    theta_temp = np.pi
    zs_temp = 1.5
    integrand_vals = i_sca_integrand(theta_temp, tau_temp, R_temp, zs_temp, lamb_vals, bc03_f)
    plt.plot(lamb_vals, integrand_vals / bc03_f(lamb_vals))
    plt.title('I_sca integrand vs. lambda')
    plt.show()
    """

    n_theta = 30
    n_rho = 27
    n_tau = 28
    n_beta = 29
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
    rho_min = 0
    rho_max = 500
    h_rho = (rho_max - rho_min) / 2 / n_rho
    rho_grid = np.linspace(rho_min + h_rho, rho_max - h_rho, n_rho)

    taus, betas, rhos, thetas = np.meshgrid(tau_grid, beta_grid, rho_grid, theta_grid, indexing='ij')
    ww = i_sca_integrand(thetas, taus, rhos, betas, lamb, bc03)

    # sum over grid
    theta_div = (theta_max - theta_min) / (n_theta - 1)
    rho_div = (rho_max - rho_min) / (n_rho - 1)
    tau_div = (tau_max - tau_min) / (n_tau - 1)
    beta_div = (beta_max - beta_min) / (n_beta - 1)
    result = np.sum(ww) * theta_div * rho_div * tau_div * beta_div
    return result


# solar metallicity: .012
# (use .008 for now)
# index 2: t=5e9, 12 Gyr, Z=.008
# 5th to last: similar but t9e9

# loop through the bc03 spectra
for p in paths[2:3]:
    # load and interpolate the bc03 spectra
    a = np.loadtxt(p)
    wav = a[:, 0]  # angstroms
    bc03 = a[:, 1]
    bc03_f = interp1d(wav, bc03, kind='cubic')

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
    rho_grid = np.linspace(0, 500, num_div)
    lamb_grid = np.linspace(100, 10**6, num_div)

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
    rho_plot = i_tir_integrand(6000, 1, rho_grid, np.pi / 4, bc03_f)
    print(rho_plot)
    plt.plot(rho_grid, rho_plot)
    plt.title('I_TIR integrand vs. rho')
    plt.show()
    """

    # compute total infrared radiation
    # (0 to inf throws an error)
    # (range of lambda for dust model is 10^-4 to 10^4 microns, or 1 to 10^8 angstroms)
    # (note: range of lambda for bc03 is 91 to 1.6e6)

    I_tir = i_tir()
    print("I_tir result")
    print(I_tir)

    # convert to 100 micron (there is an associated uncertainty)
    nu_I_nu_100 = .52 * I_tir  # units of angstroms * sigma

    print("starting integral")

    #wavelength_partial = wavelength[(wavelength > 4600) & (wavelength < 5600)]  #5380 to 5480
    wavelength_partial = [3842, 3920, 4450, 4996, 6480, 6660, 8875, 9553]
    # wavelength_partial = wavelength
    # try to use broadcasting for this?
    i_sca_array = np.array([i_sca(lamb, bc03_f) for lamb in wavelength_partial])  # units of sigma * parsecs
    i_lam_array = i_sca_array * wavelength_partial  # units of sigma * parsecs * angstroms
    alphas = i_lam_array / nu_I_nu_100

    print("Wavelength and alphas")
    print(wavelength_partial)
    print(alphas)

    plt.plot(wavelength_partial, alphas / bc03_f(wavelength_partial))
    plt.plot(wav, bc03_f(wav))
    plt.show()

    # save the result
    save_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
    save_path += p.split('/')[-1].rsplit('.spec')[0] + '.npy'
    print(save_path)
    np.save(save_path, alphas)
