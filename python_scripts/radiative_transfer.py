# TODO
# check units again, for ex. cross section is cm^2
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
dust_wav = np.flip(dust_model[:, 0].flatten()) * 10**4  # convert to angstroms from um
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
def henyey(mu, g):
    return (1 - g**2) / (1 + g - 2*g*mu)**(3/2) / 2

# UNITS: LEFT OFF HERE

sig = 250  # parsecs (from density eqn)
# prefactor: calculate based on tau(z=0) = 0.15. Use V band,
V_band_wav = 5510  # A
density_prefactor = 0.3 / dust_cross_f(V_band_wav) / sig / np.sqrt(2 * np.pi)


# eqn A1
# density: rho ~ exp(-z^2 / 2 / sig^2).
# I took sigma out of the integral and integrated analytically.
# units of cross section should cancel out b/c it's also in the density prefactor.
def tau_f(lamb, z):
    cross_section = dust_cross_f(lamb)
    tau_lambda = density_prefactor * cross_section * sig * np.sqrt(2 * np.pi) / 2 \
        * (1 - special.erf(z / sig / np.sqrt(2)))
    return tau_lambda


# reverse: find z given tau
# (see paper calculations)
# (maybe not needed any more)
def z_of_tau(lamb, tau):
    cross_section = dust_cross_f(lamb)
    return sig * np.sqrt(2) * special.erfinv(1 - 2*tau/(density_prefactor * cross_section * sig * np.sqrt(2*np.pi)))


# A2
def A_f(lamb, z, z_s, R):
    term1 = np.abs(tau_f(lamb, z) - tau_f(lamb, z_s))
    term2 = np.sqrt((z-z_s)**2 + R**2)
    term3 = np.abs(z - z_s)
    result = term1 * term2 / term3
    return result


# eqn A5
def cos_xi(z, R, theta, z_s):
    numer = (z - z_s)**2 - R * z * np.cos(theta) / np.tan(b)
    term1 = z**2 / np.tan(b)**2 + (z - z_s)**2
    term2 = R**2 + (z - z_s)**2
    denom = np.sqrt(term1 * term2)
    result = numer / denom
    print("cos xi")
    print(z)
    print(z_s)
    print(result)
    if type(z) == np.ndarray or type(z_s) == np.ndarray:
        result[z == z_s] = -np.cos(theta)  # avoid singularity; see notes
    elif z == z_s:
        return -np.cos(theta)
    return result


# scale heights 300 pc and 1350 pc
a_300 = 0.9
a_1350 = 0.1


# surface power density (any normalization factor will cancel out later)
# it's a function of height z_s
def surface_power_fn(bc03, z_s, lamb):  # z_s in pc
    return bc03(lamb) * (a_300 * np.exp(-z_s / 300) + a_1350 * np.exp(-z_s / 1350))


# I found the derivative analytically
def surface_power_deriv(bc03, z_s, lamb):
    return bc03(lamb) * (-a_300 * np.exp(-z_s / 300) / 300 - a_1350 * np.exp(-z_s / 1350) / 1350)


# eqn A4: integrand
def i_tir_integrand(lamb, z_s, bc03):
    prefactor = 1 / (8*np.pi) / np.sin(np.abs(b))
    term1 = 1 - dust_albedo_f(lamb)
    term2 = surface_power_deriv(bc03, z_s, lamb)
    tau_zs = tau_f(lamb, z_s)
    tau_0 = tau_f(lamb, 0)
    term3 = 2 + tau_zs*special.expi(tau_zs) - np.exp(tau_zs) \
            + (tau_0 - tau_zs)*special.expi(tau_0 - tau_zs) - np.exp(tau_0 - tau_zs)
    result = prefactor * term1 * term2 * term3
    return result

# eqn A7, part 1: just the integrand.
# note that we need to integrate over z_s also.
def i_sca_integrand(theta, tau, R, z_s, lamb, bc03):
    prefactor = (1 / np.sin(np.abs(b))) * dust_albedo_f(lamb)
    z = z_of_tau(lamb, tau)
    term1 = np.exp((-1 / np.sin(np.abs(b))) * (tau_f(lamb, 0) - tau))
    term2 = R
    term3 = henyey(cos_xi(z, R, theta, z_s), dust_cos_f(lamb))
    term4_exp = np.exp(-A_f(lamb, z, z_s, R))
    if type(term4_exp) == np.ndarray:
        term4_exp[z == z_s] = np.exp(-R*density_prefactor*dust_cross_f(lamb)) / R**2
    elif z == z_s:
        term4_exp = 0
    term4 = surface_power_fn(bc03, z_s, lamb) * term4_exp
    term5 = 4 * np.pi * ((z - z_s)**2 + R**2)
    result = prefactor * term1 * term2 * term3 * term4 / term5
    return result


# eqn A7, part 2: do the integral for given wavelength
def i_sca(lamb, bc03):

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
    tau_vals = np.linspace(0.0001, tau_f(lamb, 0)-.0001, 50)
    theta_temp = np.pi
    R_temp = 1
    zs_temp = 1.5
    integrand_vals = i_sca_integrand(theta_temp, tau_vals, R_temp, zs_temp, lamb, bc03_f)
    plt.plot(tau_vals, integrand_vals)
    plt.title('I_sca integrand vs. tau')
    plt.show()

    # plot vs. z_s (for given tau, theta, R)
    zs_vals = np.linspace(50, 100, 50)
    tau_temp = 0.1
    theta_temp = np.pi
    R_temp = 1
    integrand_vals = i_sca_integrand(theta_temp, tau_temp, R_temp, zs_vals, lamb, bc03_f)
    plt.plot(zs_vals, integrand_vals)
    plt.title('I_sca integrand vs. z_s')
    plt.show()

    # plot vs. R (for given tau theta, zs)
    R_vals = np.linspace(0.01, 1000, 50)
    tau_temp = 0.1
    theta_temp = np.pi
    zs_temp = 1.5
    integrand_vals = i_sca_integrand(theta_temp, tau_temp, R_vals, zs_temp, lamb, bc03_f)
    plt.plot(R_vals, integrand_vals)
    plt.title('I_sca integrand vs. R')
    plt.show()

    x = np.linspace(0, 2*np.pi, 50)  # theta grid, 0 to 2pi
    y = np.linspace(.0001, tau_f(lamb, 0)-.0001, 50)  # tau grid, 0 to tau(0)
    z = np.linspace(.01, 1000, 50)  # R grid
    v = np.linspace(-3, 1000, 50)  # z_s grid, -inf to inf  (peak depends on tau)

    xx, yy, zz, vv = np.meshgrid(x, y, z, v)
    ww = i_sca_integrand(xx, yy, zz, vv, lamb, bc03)
    # print("number of nans")
    # print(np.count_nonzero(np.isnan(ww)))
    ww.shape
    inner = [integrate.simps(ww_x, x) for ww_x in ww]
    middle1 = [integrate.simps(ww_y, y) for ww_y in inner]
    middle2 = [integrate.simps(ww_z, z) for ww_z in middle1]
    result = integrate.simps(middle2, v)
    #print("ww")
    #print(ww)
    #print("inner")
    #print(inner)
    return result


# loop through the bc03 spectra
for p in paths[:1]:
    # load and interpolate the bc03 spectra
    a = np.loadtxt(p)
    wav = a[:, 0]  # angstroms
    bc03 = a[:, 1]
    bc03_f = interp1d(wav, bc03, kind='cubic')

    # compute total infrared radiation
    # (integrate from 91 to 10^4 bc 0 to inf throws an error)
    # (range of lambda for dust model is 10^-4 to 10^4)
    # (note: range of lambda for bc03 is 91 to 1.6e6)

    # plot integrand for I_TIR:
    # vars: z, lambda, z_s
    # plot vs. lambda (for given z, zs)
    lambda_vals = np.linspace(100, 10 ** 4, 50)
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


    y = np.linspace(100, 10**4, 50)  # min, max, n   # lambda grid, 91 to 10**4 A
    z = np.linspace(-3, 70, 50)  # z_s grid, -inf to inf
    yy, zz = np.meshgrid(y, z)
    ww = i_tir_integrand(yy, zz, bc03_f)
    ww.shape
    inner = [integrate.simps(ww_y, y) for ww_y in ww]
    I_tir = integrate.simps(inner, z)

    print("I_tir result")
    print(I_tir)

    # convert to 100 micron (there is an associated uncertainty)
    nu_I_nu_100 = .52 * I_tir

    # REVIEWED THROUGH HERE

    print("starting integral")

    # this will take a long time
    # remember to multiply by lambda (TODO)
    I_sca = i_sca(6000, bc03_f)
    print(I_sca)
    wavelength_partial = wavelength[1000:1001]

    i_sca_array = np.array([i_sca(lamb, bc03_f) for lamb in wavelength_partial])

    i_lam_array = i_sca_array * wavelength_partial
    alphas = i_lam_array / nu_I_nu_100

    print(alphas)

    plt.plot(wavelength_partial, alphas)
    plt.show()

    # save the result
    save_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
    save_path += p.split('/')[-1].rsplit('.spec')[0] + '.npy'
    # np.save(save_path, result)

    # don't forget to multiply by lambda
