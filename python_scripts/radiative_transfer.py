#TODO
# check my scratch work and code (units, argument order)
# look at the estimated errors of the integrals
#   other checks: try MCMC method?
# (which latitude (b)?)
# phase function normalization
# integrals: can I just stop at the wavelength bounds? or best to estimate what the integral is outside?
# sanity checks?


# take the BC03 model spectra and apply radiative transfer calculations.
# based on the appendix of BD12.

import glob
import numpy as np
from scipy.interpolate import interp1d
import scipy.special as special
from scipy import integrate

# load some files
wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')
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
dust_cross = np.flip(dust_model[:, 3].flatten())
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
dust_cross_f = interp1d(dust_wav, dust_cross, kind='cubic')

# choose a latitude
b = 40 * np.pi / 180  # 40 degrees -> radians


# phase function
def henyey(mu, g):
    return (1 - g**2) / (1 + g - 2*g*mu)**(3/2) / 2


density_prefactor = 0.15
sig = 250  # parsecs (from density eqn)


# eqn A1
# density: rho ~ exp(-z^2 / 2 / sig^2).
# I took sigma out of the integral and integrated analytically.
def tau_f(lamb, z):
    cross_section = dust_cross_f(lamb)
    tau_lambda = density_prefactor * cross_section * sig * np.sqrt(2 * np.pi) / 2 \
        * (1 - special.erf(z / sig / np.sqrt(2)))
    return tau_lambda


# reverse: find z given tau
# (see paper calculations)
def z_of_tau(lamb, tau):
    cross_section = dust_cross_f(lamb)
    return sig * np.sqrt(2) * special.erfinv(1 - 2*tau/(density_prefactor * cross_section * sig * np.sqrt(2*np.pi)))


# A2
def A_f(lamb, z, z_s, R):
    term1 = np.abs(tau_f(lamb, z) - tau_f(lamb, z_s))
    term2 = np.sqrt((z-z_s)**2 + R**2)
    term3 = np.abs(z - z_s)
    return term1 * term2 / term3


# eqn A5
def cos_xi(tau, R, theta, z_s, lamb):
    z = z_of_tau(lamb, tau)
    numer = (z - z_s)**2 - R * z * np.cos(theta) / np.tan(b)
    term1 = z**2 / np.tan(b)**2 + (z - z_s)**2
    term2 = R**2 + (z - z_s)**2
    denom = np.sqrt(term1 * term2)
    return numer / denom


# eqn A4, part 1 (just the inner integral)
def i_tir_inner(lamb, z_s):
    def integrand(tau):
        return special.expi(np.abs(tau - tau_f(lamb, z_s)))
    result, err = integrate.quad(integrand, 0, tau_f(lamb, 0))
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


# eqn A4, part 2 (calculate the integrand, but don't do the outer 2 integrals)
def i_tir_integrand(lamb, z_s, bc03):
    prefactor = 1 / (8*np.pi) / np.sin(np.abs(b))
    term1 = 1 - dust_albedo_f(lamb)
    term2 = surface_power_deriv(bc03, z_s, lamb)
    term3 = i_tir_inner(lamb, z_s)
    return prefactor * term1 * term2 * term3


# eqn A7, part 1: just the integrand.
# note that we need to integrate over z_s also.
def i_sca_integrand(tau, R, theta, z_s, lamb, bc03):
    z = z_of_tau(lamb, tau)

    # prefactor is accounted for in wrapper function
    term1 = np.exp((-1 / np.sin(np.abs(b))) * (tau_f(lamb, 0) - tau))
    term2 = R
    term3 = henyey(cos_xi(tau, R, theta, z_s, lamb), dust_cos_f(lamb))
    term4 = surface_power_fn(bc03, z_s, lamb) * np.exp(-A_f(lamb, z, z_s, R))
    term5 = 4 * np.pi * ((z - z_s)**2 + R**2)
    return term1 * term2 * term3 * term4 / term5


# eqn A7, part 2: do the integral for given wavelength
def i_sca(lamb, bc03):
    result, err = integrate.nquad(i_sca_integrand, [[0, tau_f(lamb, 0)],
                                                   [0, np.inf],
                                                   [0, 2*np.pi],
                                                   [-np.inf, np.inf]], args=(lamb, bc03))
    prefactor = (1 / np.sin(np.abs(b))) * dust_albedo_f(lamb)
    return prefactor * result


# loop through the bc03 spectra
for p in paths[:1]:
    # load and interpolate the bc03 spectra
    a = np.loadtxt(p)
    wav = a[:, 0]  # angstroms
    bc03 = a[:, 1]
    bc03_f = interp1d(wav, bc03, kind='cubic')

    """
    # compute total infrared radiation
    # (integrate from 91 to 10^4 bc 0 to inf throws an error)
    # (range of lambda for dust model is 10^-4 to 10^4)
    # (note: range of lambda for bc03 is 91 to 1.6e6)
    I_tir, I_tir_err = integrate.nquad(i_tir_integrand, [[91, 10**4],
                                                         [-np.inf, np.inf]], args=(bc03_f,))
    # convert to 100 micron (there is an associated uncertainty)
    nu_I_nu_100 = .52 * I_tir
    """


    # REVIEWED THROUGH HERE

    print("tau from z")
    tau_test = tau_f(6000, 1)
    print(z_of_tau(6000, tau_test))

    temp = i_sca_integrand(tau=1, R=1, theta=np.pi, z_s=1, lamb=6000, bc03=bc03_f)
    print("success")
    print(temp)
    z = z_of_tau(lamb=6000, tau=1)
    print(z)

    cross_section = dust_cross_f(6000)
    print(cross_section)
    tau = 1
    erf_arg = 1 - 2 * tau / (density_prefactor * cross_section * sig * np.sqrt(2 * np.pi))
    print(erf_arg)
    erf_result = special.erfinv(erf_arg)
    print(erf_result)

    # sig * np.sqrt(2) * special.erfinv(1 - 2 * tau / (density_prefactor * cross_section * sig * np.sqrt(2 * np.pi)))

    #term1 = np.exp((-1 / np.sin(np.abs(b))) * (tau_f(lamb, 0) - tau))
    #term2 = R
    #term3 = henyey(cos_xi(tau, R, theta, z_s, lamb), dust_cos_f(lamb))
    #term4 = surface_power_fn(bc03, z_s, lamb) * np.exp(-A_f(lamb, z, z_s, R))
    #term5 = 4 * np.pi * ((z - z_s) ** 2 + R ** 2)
    #result = term1 * term2 * term3 * term4 / term5
    #print(result)

    print("starting integral")
    test_lambda = 6000
    i_sca(6000, bc03_f)
    exit(0)

    # this will take a long time
    # result_spectrum = np.array([i_sca(lamb, bc03) for lamb in wavelength[:1]]) / nu_I_nu_100

    # save the result
    save_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
    save_path += p.split('/')[-1].rsplit('.spec')[0] + '.npy'
    # np.save(save_path, result)

    # don't forget to multiply by lambda
