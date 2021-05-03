# brandt's c code for radiative transfer outputs to a .dat file.
# here I read the data and plot the alphas.

# I've moved a copy to the dust_scripts folder so I can back up to github.

import numpy as np
import matplotlib.pyplot as plt
import os
from generate_plots import generate_binned_alphas

dust_model = 'wd'  # wd

a_300 = 0.9
a_1350 = 0.1
sig_star = 1.2  # hard coded in the c code

# set up grid of f values (which transform into z_s)
n_f = 30
f_min, f_max = 0, 1
h_f = (f_max - f_min) / n_f / 2
f_grid = np.linspace(f_min + h_f, f_max - h_f, n_f)

# set up arrays to store total integral
# and get the wavelengths
# column 0: lambdas
# column 1: tau0
# column 2: sca integral
# column 3: tir integral
# column 4: albedo
if dust_model == 'zd':
    test_data = np.loadtxt('dustsheet_0.15_40_0.5.dat')
elif dust_model == 'wd':
    test_data = np.loadtxt('dustsheet_0.15_40_0.55_wd01.dat')
total_int = np.zeros(test_data.shape[0])
total_int_denom = np.zeros(test_data.shape[0])
lambdas = test_data[:, 0]

# get bc03 spectrum
import glob
from scipy.interpolate import interp1d
paths = glob.glob('/Users/blakechellew/Documents/DustProject/BrandtFiles/bc03/*.spec')
for p in paths[22:23]:
    print("verify path name:", p)
    a = np.loadtxt(p)
    wav = a[:, 0]  # angstroms
    bc03 = a[:, 1]
    bc03_f = interp1d(wav, bc03, kind='cubic')

# integrate
# maybe run "make" first? in case I forgot.
for f in f_grid:
    if dust_model == 'zd':
        load_path = 'dustsheet_0.15_40_' + str(f) + '.dat'
        if not os.path.isfile(load_path):
            os.system('./sheetint 0.15 40 ' + str(f))
    elif dust_model == 'wd':
        load_path = 'dustsheet_0.15_40_' + str(f) + '_wd01.dat'
        if not os.path.isfile(load_path):
            os.system('./sheetint 0.15 40 ' + str(f) + ' 2')

    data = np.loadtxt(load_path)
    intval = data[:, 2]
    intval_denom = data[:, 3]
    if f < .5:
        transform_factor = sig_star / f
        z_s = sig_star * np.log(2*f)
    else:
        transform_factor = sig_star / (1-f)
        z_s = -sig_star * np.log(2*(1-f))
    print("z_s:", z_s)

    # add to the total integral:
    surface_power_factor = a_300 * np.exp(-np.abs(z_s)/300) + a_1350 * np.exp(-np.abs(z_s)/1350)
    df = (f_max - f_min) / n_f
    dz = 2 * sig_star * np.exp(np.abs(z_s) / sig_star) * df
    print("integral result:", intval)
    print("integral denom:", intval_denom)
    print("surface_power_factor:", surface_power_factor)
    total_int += intval * surface_power_factor * dz
    total_int_denom += intval_denom * surface_power_factor * dz

# add up I_tir:
I_tir = np.sum(total_int_denom * bc03_f(lambdas))
nu_I_nu = I_tir * 0.52

# plot various things
# TODO: need to plot against BOSS wavelengths with interpolation and binning (to see where the peak really is)
# TODO: seems like the z_s range might not be large enough. only goes up to 6 pc I think.
# TODO: actually should multiply by derivative of surface_power_factor (I think...)
# TODO: reproduce the results of the other dust model also

# load BOSS wavelengths, prep for interpolation
boss_lambdas = np.load('/Users/blakechellew/Documents/DustProject/alphas_and_stds/wavelength_boss.npy')  # angstroms
total_int_f = interp1d(lambdas, total_int, kind='cubic')

if dust_model == 'zd':
    bd12_factor = 0.52
elif dust_model == 'wd':
    bd12_factor = 0.49
plot1 = lambdas * bc03_f(lambdas)
plot2 = total_int / nu_I_nu
plot3 = bd12_factor * boss_lambdas * total_int_f(boss_lambdas) / nu_I_nu * bc03_f(boss_lambdas)

# bin the main plot
lambdas_boss_bin, plot3_bin, stds = generate_binned_alphas([plot3], [np.ones(plot3.shape)], boss_lambdas)
plot3_bin = plot3_bin[0]

plot1_mean = np.mean(plot1)
plot2_mean = np.mean(plot2)
plot3_mean = np.mean(plot3)
print("plot1", plot1)
print("plot2", plot2)
print("plot3", plot3)
print("plot 1 mean:", plot1_mean)
print("plot 2 mean:", plot2_mean)
print("plot 3 mean:", plot3_mean)
plt.plot(lambdas, plot1 * plot3_mean / plot1_mean, label='bc03 * wav')
plt.plot(lambdas, plot2 * plot3_mean / plot2_mean, label='just the integral')
plt.plot(lambdas_boss_bin, plot3_bin, step='pre', label='predicted alphas')
plt.legend()
plt.xlim(3800, 9200)
plt.show()