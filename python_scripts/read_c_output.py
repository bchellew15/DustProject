# brandt's c code for radiative transfer outputs to a .dat file.
# here I read the data and plot the alphas.

# I've moved a copy to the dust_scripts folder so I can back up to github.

# TODO: update the binning
# TODO: cross section should be differet in V band for each dust model...
# TODO: seems like the z_s range might not be large enough. only goes up to 6 pc I think...
# TODO: do I even need to use the SDSS wavelengths at all? maybe it would be better to use a uniform grid
# (and anyway, the real values are weighted and these aren't... should I weight these?)

import numpy as np
import matplotlib.pyplot as plt
import os
from generate_plots import generate_binned_alphas
from astropy.io import fits

dust_model = 'zd'  # wd, zd
folder_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/brandt_radiative/integrals/'
boss = False
regen_files = True

a_300 = 1  # 0.9
a_1350 = 0  # 0.1
sig_star = 1.2  # 1.2  # hard coded in the c code
path_num = 2  # see radiative_transfer_v2. 22 is z=.02, 2 is z=.08
sig_star_1 = sig_star  # 300
sig_star_2 = 1350

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
    file_length = 334
elif dust_model == 'wd':
    file_length = 203
total_int = np.zeros(file_length)
total_int_denom = np.zeros(file_length)

# get bc03 spectrum
import glob
from scipy.interpolate import interp1d
paths = glob.glob('/Users/blakechellew/Documents/DustProject/BrandtFiles/bc03/*.spec')
for p in paths[path_num:path_num+1]:
    print("verify path name:", p)
    a = np.loadtxt(p)
    wav = a[:, 0]  # angstroms
    bc03 = a[:, 1]
    bc03_f = interp1d(wav, bc03, kind='cubic')

lambdas = None  # get the wavelengths from one of the files in the loop below

# integrate
# but run "make" first in case I forgot
os.system('make')

for f in f_grid:
    if dust_model == 'zd':
        load_path = folder_path + 'dustsheet_0.15_40_' + str(f) + '.dat'
        if regen_files or not os.path.isfile(load_path):
            os.system(folder_path + 'sheetint 0.15 40 ' + str(f))
    elif dust_model == 'wd':
        load_path = folder_path + 'dustsheet_0.15_40_' + str(f) + '_wd01.dat'
        print("load path:", load_path)
        if regen_files or not os.path.isfile(load_path):
            print("running")
            os.system(folder_path + 'sheetint 0.15 40 ' + str(f) + ' 2')

    print("about to load data")
    data = np.loadtxt(load_path)
    intval = data[:, 2]
    intval_denom = data[:, 3]
    if f < .5:
        # transform_factor = sig_star / f  # this is included as part of dz
        z_s = sig_star * np.log(2*f)
    else:
        # transform_factor = sig_star / (1-f)  # this is included as part of dz
        z_s = -sig_star * np.log(2*(1-f))
    print("z_s:", z_s)

    # add to the total integral:
    surface_power_factor = a_300 * np.exp(-np.abs(z_s)/sig_star_1) + a_1350 * np.exp(-np.abs(z_s)/sig_star_2)
    df = (f_max - f_min) / n_f
    dz = 2 * sig_star * np.exp(np.abs(z_s) / sig_star) * df
    print("integral result:", intval)
    print("integral denom:", intval_denom)
    print("surface_power_factor:", surface_power_factor)
    total_int += intval * surface_power_factor * dz
    total_int_denom += intval_denom * surface_power_factor * dz

    if lambdas is None:
        lambdas = data[:, 0]


# add up I_tir:
d_lambdas = -np.diff(lambdas) # because they're decreasing
d_lambdas = np.append(d_lambdas, d_lambdas[-1])
assert total_int_denom.ndim == 1
I_tir = np.sum(total_int_denom * bc03_f(lambdas) * d_lambdas)
nu_I_nu = I_tir * 0.52

# load survey wavelengths, prep for interpolation
if boss:
    survey_lambdas = np.load('/Users/blakechellew/Documents/DustProject/alphas_and_stds/wavelength_boss.npy')  # angstroms
else:
    hdulist_direc = '/Users/blakechellew/Documents/DustProject/BrandtFiles/'
    hdulist_sdsswav = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
    survey_lambdas = np.array(hdulist_sdsswav[1].data)
total_int_f = interp1d(lambdas, total_int, kind='cubic')

if dust_model == 'zd':
    bd12_factor = 0.52
elif dust_model == 'wd':
    bd12_factor = 0.49
plot1 = lambdas * bc03_f(lambdas)
plot2 = total_int / nu_I_nu
plot3 = bd12_factor * survey_lambdas * total_int_f(survey_lambdas) / nu_I_nu * bc03_f(survey_lambdas)

# bin the main plot
lambdas_bin, plot3_bin, stds = generate_binned_alphas([plot3], [np.ones(plot3.shape)], survey_lambdas)
plot3_bin = plot3_bin[0]

# load bd12 plot for comparison
if dust_model == 'wd':
    bd12_plot = np.loadtxt('/Users/blakechellew/Documents/DustProject/alphas_and_stds/bd12_fig3_green_052621.csv',
                            delimiter=",")
else:
    bd12_plot = np.loadtxt('/Users/blakechellew/Documents/DustProject/alphas_and_stds/bd12_fig3_blue_052621.csv',
                        delimiter=",")

plot1_mean = np.mean(plot1)
plot2_mean = np.mean(plot2)
plot3_mean = np.mean(plot3)
print("plot1", plot1)
print("plot2", plot2)
print("plot3", plot3)
print("plot 1 mean:", plot1_mean)
print("plot 2 mean:", plot2_mean)
print("plot 3 mean:", plot3_mean)
plt.plot(bd12_plot[:, 0], bd12_plot[:, 1], 'green', label='bd12 wd01', drawstyle='steps-mid')
plt.plot(lambdas, plot1 * plot3_mean / plot1_mean, label='bc03 * wav')
plt.plot(lambdas, plot2 * plot3_mean / plot2_mean, label='just the integral')
plt.plot(lambdas_bin, plot3_bin, 'red', drawstyle='steps', label='predicted alphas')
plt.plot(lambdas, total_int_f(lambdas), 'purple', drawstyle='steps', label='Just the numerator')
plt.legend()
plt.xlim(3800, 9200)
plt.ylim(0, .28)
plt.show()

wd_path = folder_path + 'alphas_wd.npy'
zd_path = folder_path + 'alphas_zd.npy'
if os.path.isfile(wd_path) and os.path.isfile(zd_path):
    alphas_wd, alphas_zd = np.load(wd_path), np.load(zd_path)
    lambdas_bin, alphas_wd_bin, _ = generate_binned_alphas([alphas_wd], [np.ones(alphas_wd.shape)], survey_lambdas)
    _, alphas_zd_bin, _ = generate_binned_alphas([alphas_zd], [np.ones(alphas_zd.shape)], survey_lambdas)
    alphas_wd_bin, alphas_zd_bin = alphas_wd_bin[0], alphas_zd_bin[0]

    plt.plot(lambdas_bin, alphas_wd_bin, drawstyle='steps', label='wd01 alphas')
    plt.plot(lambdas_bin, alphas_zd_bin, drawstyle='steps', label='zda04 alphas')
    plt.plot(bd12_plot[:, 0], bd12_plot[:, 1], 'green', label='bd12 wd01', drawstyle='steps-mid')
    plt.legend()
    plt.xlim(3800, 9200)
    plt.ylim(0, 0.28)
    plt.show()

# save alphas
if dust_model == 'wd':
    np.save(wd_path, plot3)
else:
    np.save(zd_path, plot3)