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
# (maybe not needed any more)
def z_of_tau(lamb, tau):
    cross_section = dust_cross_f(lamb)  # cm^2
    return sig_dust * np.sqrt(2) * special.erfinv(
        1 - 2 * tau / (density_prefactor * cross_section * sig_dust * np.sqrt(2 * np.pi)))


# A2
def A_f(lamb, z, rho, beta):
    term1 = tau_f(lamb, z) - tau_f(lamb, z - rho * np.cos(beta))
    term2 = np.cos(beta)
    result = np.abs(term1 / term2)

    if type(z) == np.ndarray or type(z_s) == np.ndarray:
        np.putmask(result, np.cos(beta) == 0, density_prefactor * dust_cross_f(lamb) * rho * np.sin(beta) * np.exp(-z ** 2 / 2 / sig_dust ** 2))
    elif np.cos(beta) == 0:
        return density_prefactor * dust_cross_f(lamb) * rho * np.sin(beta) * np.exp(-z ** 2 / 2 / sig_dust ** 2)

    return result


# eqn A5
def cos_xi(z, rho, theta, beta):
    numer = rho * np.cos(beta)**2 - z * np.sin(beta) * np.cos(theta) / np.tan(b)
    denom_sqr = z**2 / np.tan(b)**2 + rho**2 * np.cos(beta)**2
    return numer / np.sqrt(denom_sqr)


# scale heights 300 pc and 1350 pc
a_300 = 0.9
a_1350 = 0.1


# surface power density (arbitrary units; any normalization factor will cancel out later)
# it's a function of height z_s
# def surface_power_fn(bc03, z_s, lamb):  # z_s in pc
#     return bc03(lamb) * (a_300 * np.exp(-z_s / 300) + a_1350 * np.exp(-z_s / 1350))


# I found the derivative analytically
def surface_power_deriv(bc03, z, rho, beta, lamb):
    return bc03(lamb) * \
           (-a_300 * np.exp(-(z - rho *np.cos(beta)) / 300) / 300
            - a_1350 * np.exp(-(z - rho *np.cos(beta)) / 1350) / 1350)


# eqn A4: integrand
# units: angstroms * units of sigma after integration
def i_tir_integrand(lamb, z, rho, beta, bc03):
    prefactor = 1 / (8 * np.pi) / np.sin(np.abs(b))
    term1 = 1 - dust_albedo_f(lamb)
    term2 = surface_power_deriv(bc03, z, rho, beta, lamb)
    term3 = np.exp(-A_f(lamb, z, rho, beta))
    term4 = np.sin(beta)

    # apply neg. sign:
    term3 = -term3

    result = prefactor * term1 * term2 * term3 * term4
    return result

"""
# eqn A7, part 1: just the integrand.
# note that we integrate over z_s also.
# units: parsecs * units of sigma
def i_sca_integrand(theta, tau, rho, beta, lamb, bc03):
    prefactor = (1 / np.sin(np.abs(b))) * dust_albedo_f(lamb)
    z = z_of_tau(lamb, tau)

    # check for tau = 0
    if type(tau) != np.ndarray and tau == 0:
        return 0  # this happens when z = infinity
    z[tau == 0] = 1  # just so we don't have z = inf in cos_xi

    term1 = np.exp((-1 / np.sin(np.abs(b))) * (tau_f(lamb, 0) - tau))
    term2 = R
    term3 = henyey(cos_xi(z, rho, theta, beta), dust_cos_f(lamb))
    term4_exp = np.exp(-A_f(lamb, z, rho, beta))
    term4 = -surface_power_deriv(bc03, z, rho, beta, lamb) * term4_exp  # TEMP: added - and switched to deriv
    term5 = 4 * np.pi * ((z - z_s) ** 2 + R ** 2)
    result = prefactor * term1 * term2 * term3 * term4 / term5

    if type(tau) == np.ndarray:
        np.putmask(result, tau == 0, 0)

    return result
"""

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

    num_div = 50
    x_min = 0
    x_max = 2 * np.pi
    x = np.linspace(x_min, x_max, num_div)  # theta grid, 0 to 2pi
    y_min = 0  # .01  # get better estimate of lower bound. was using .01.
    y_max = tau_f(lamb, 0)  # -.00001
    y = np.linspace(y_min, y_max, num_div)  # tau grid, 0 to tau(0)
    z_min = .01
    z_max = 1000
    z = np.linspace(z_min, z_max, num_div)  # R grid
    v_min = -100
    v_max = 1000
    v = np.linspace(v_min, v_max, num_div)  # z_s grid, -inf to inf  (peak depends on tau)

    xx, yy, zz, vv = np.meshgrid(x, y, z, v)
    ww = i_sca_integrand(xx, yy, zz, vv, lamb, bc03)

    # try basic addition on the grid:
    x_div = (x_max - x_min) / (num_div - 1)
    y_div = (y_max - y_min) / (num_div - 1)
    z_div = (z_max - z_min) / (num_div - 1)
    v_div = (v_max - v_min) / (num_div - 1)
    result = np.sum(ww) * x_div * y_div * z_div * v_div
    return result

    # print("number of nans")
    # print(np.count_nonzero(np.isnan(ww)))

    # inner = [integrate.simps(ww_x, x) for ww_x in ww]
    # middle1 = [integrate.simps(ww_y, y) for ww_y in inner]
    # middle2 = [integrate.simps(ww_z, z) for ww_z in middle1]
    # result = integrate.simps(middle2, v)
    # return result


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

    # testing A4:
    # i_tir = i_tir_integrand(6000, 2, bc03_f)
    # print("I tir test")
    # print(i_tir)
    # exit(0)

    # compute total infrared radiation
    # (0 to inf throws an error)
    # (range of lambda for dust model is 10^-4 to 10^4 microns, or 1 to 10^8 angstroms)
    # (note: range of lambda for bc03 is 91 to 1.6e6)

    """
    # plot integrand for I_TIR:
    # vars: z, lambda, z_s
    # plot vs. lambda (for given z, zs)
    lambda_vals = np.linspace(91, 6*10**4, 50)
    z_s_temp = 1.5
    plt.plot(lambda_vals, i_tir_integrand(lambda_vals, z_s_temp, bc03_f))
    plt.title('I_TIR integrand vs. wavelength')
    plt.show()

    # plot vs. z_s (for given lambda)
    zs_vals = np.linspace(-3, 70, 50)
    lambda_temp = 6000
    z_temp = 1.5
    plt.plot(zs_vals, i_tir_integrand(lambda_temp, zs_vals, bc03_f))
    plt.title('I_TIR integrand vs. z_s')
    plt.show()
    """

    num_div = 50
    lamb_min = 100
    lamb_max = 10 ** 6  # 6*10**4
    lamb_grid = np.linspace(lamb_min, lamb_max, num_div)  # lambda grid, 91 to 10**4 A  # min, max, n
    tau_min = 0
    tau_max = tau_f(lamb, 0)
    tau_grid = np.linspace(y_min, y_max, num_div)  # z_s grid, -inf to inf
    beta_min =
    beta_max =
    beta_grid =
    rho_min = 0
    rho_max = 1000
    rho_grid =
    yy, zz = np.meshgrid(y, z)
    ww = i_tir_integrand(yy, zz, bc03_f)
    ww.shape

    # integrate on simple grid
    y_div = (y_max - y_min) / (num_div - 1)
    z_div = (z_max - z_min) / (num_div - 1)
    I_tir = np.sum(ww) * y_div * z_div  # units of angstroms * sigma

    # inner = [integrate.simps(ww_y, y) for ww_y in ww]
    # I_tir = integrate.simps(inner, z)

    print("I_tir result")
    print(I_tir)

    # convert to 100 micron (there is an associated uncertainty)
    nu_I_nu_100 = .52 * I_tir  # units of angstroms * sigma

    exit(0)

    print("starting integral")

    # wavelength_partial = wavelength[(wavelength > 4600) & (wavelength < 5600)]  #5380 to 5480
    wavelength_partial = [3842, 3920, 4450, 4996, 6480, 6660, 8875, 9553]
    # wavelength_partial = wavelength
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
