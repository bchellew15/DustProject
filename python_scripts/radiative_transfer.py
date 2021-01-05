# tues: start 3:30

#TODO
# change instances of b to |b| where necessary

# take the BC03 model spectra and apply radiative transfer calculations.
# based on the appendix of BD12.

import glob
import numpy as np
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import scipy.special as special

# load some files
wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')
paths = glob.glob('/Users/blakechellew/Documents/DustProject/BrandtFiles/bc03/*.spec')

# load dust models
dust_model_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/kext_albedo_WD_MW_3.1_60_D03.all'
# column 0: wavelength (microns)
# column 1: albedo
# column 2: <cos>
# column 3: cross section (cm^2)
# column 5: <cos^2>

dust_model = np.loadtxt(dust_model_path, skiprows=80, usecols=(0, 1, 2, 3, 4, 5))

# flip because wavelength goes largest to smallest
dust_wav = np.flip(dust_model[:, 0].flatten()) * 10**4  # convert to angstroms from um
dust_albedo = np.flip(dust_model[:, 1].flatten())
dust_cross = np.flip(dust_model[:, 3].flatten())
# sort by wavelength; some out of order:
sorted_idx = np.argsort(dust_wav)
dust_wav = dust_wav[sorted_idx]
dust_albedo = dust_albedo[sorted_idx]
dust_cross = dust_cross[sorted_idx]
#remove duplicate wavelength
bad_idx = 125
dust_wav = np.delete(dust_wav, bad_idx)
dust_albedo = np.delete(dust_albedo, bad_idx)
dust_cross = np.delete(dust_cross, bad_idx)

# interpolate dust model
dust_albedo_f = interp1d(dust_wav, dust_albedo, kind='cubic')
dust_albedo_interp = np.array([dust_albedo_f(w) for w in wavelength])
dust_cross_f = interp1d(dust_wav, dust_cross, kind='cubic')
dust_cross_interp = np.array([dust_cross_f(w) for w in wavelength])

# choose latitude
b = 40 * np.pi / 180  # 40 degrees -> radians

# phase function
def henyey(mu):
    g = 1
    return (1 - g**2) / (1 + g - 2*g*mu)**(3/2) / 2

# A1
# density: rho ~ exp(-z^2 / 2 / sig^2).
# I took sigma out of the integral and integrated analytically.
def tau_f(lamb, z):
    cross_section = dust_cross_f(lamb)
    sig = 250  # parsecs (from density eqn)
    density_prefactor = 0.15
    tau_lambda = density_prefactor * cross_section * sig * np.sqrt(2 * np.pi) / 2 \
        * (1 - special.erf(z / sig / np.sqrt(2)))
    return tau_lambda

# A2
def A_f(lamb, z, z_s, R):
    term1 = np.abs(tau_f(lamb, z) - tau_f(lamb, z_s))
    term2 = np.sqrt((z-z_s)**2 + R**2)
    term3 = np.abs(z - z_s)
    return term1 * term2 / term3

for p in paths[:1]:
    a = np.loadtxt(p)
    wav = a[:, 0]  # I think angstroms
    bc03 = a[:, 1]

    # interpolate the bc03 spectra
    bc03_f = interp1d(wav, bc03, kind='cubic')
    bc03_interp = np.array([bc03_f(w) for w in wavelength])

    # surface power density
    # scale heights 300 pc and 1350 pc
    # it's a function of height z_s
    a_300 = 1
    a_1350 = 1
    def surface_power_fn(z_s, lamb):  # z_s in pc
        return bc03_f(lamb) * (a_300 * np.exp(-z_s / 300) + a_1350 * np.exp(-z_s / 1350))

    def surface_power_deriv(z_s, lamb):
        return bc03_f(lamb) * (-a_300 * np.exp(-z_s / 300) / 300 - a_1350 * np.exp(-z_s / 1350) / 1350)

    # eqn A4
    def I_tir_f():
        prefactor = 1 / (8*np.pi) / np.sin(np.abs(b))
        term1 = integrate.quad(1 - dust_albedo_f, 0, np.inf)
        term2 = integrate.quad(lambda z_s: surface_power_deriv(z_s, ))


    # apply radiative transfer
    # eqn A7
    (1 / np.sin(b)) * dust_albedo_interp
    # result = ...

    # save the result
    save_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
    save_path += p.split('/')[-1].rsplit('.spec')[0] + '.npy'
    # np.save(save_path, result)