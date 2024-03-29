# check convergence and sensitivity to different params
# also make plots and calculations for ERE

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from generate_plots import generate_binned_alphas
from astropy.io import fits
import pickle
from scipy.interpolate import interp1d

# paramters
save = False
save_4000 = False
load = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/convergence_test/'
# set bin width = 10 for finding the peak and width. else 50.
bin_width = 50
bootstrap = True
cst_coeff = .3  # .4 is 18.4% CST. .3 is 14.5%. .5 is 22%.
cst_coeff_wd = 0.6 # 0.8
# set wav_clipped range = 5000 to 8000 for ratios, bc that's the red band.
# for integrating the ERE use 4500 to 8500. (or 8000 for north? maybe not)
wav_clipped_min = 5000  # 4500
wav_clipped_max = 8000  # 8500
bias_factor = 2.1
third_panel = False
ere_bootstrap = False

# matplotlib options
matplotlib.rcParams['axes.labelsize'] = 'large'
matplotlib.rcParams['xtick.labelsize'] = 'large'
matplotlib.rcParams['ytick.labelsize'] = 'large'
matplotlib.rcParams['legend.fontsize'] = 'large'
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['ytick.left'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['ytick.major.size'] = 6.0
matplotlib.rcParams['ytick.minor.size'] = 3.0
matplotlib.rcParams['xtick.minor.size'] = 3.0
matplotlib.rcParams['xtick.major.size'] = 6.0
matplotlib.rcParams['axes.unicode_minus'] = False

##########################
# Load alphas, apply correction factors, etc.
##########################

wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')
hdulist_sdsswav = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
wavelength_sdss = np.array(hdulist_sdsswav[1].data)

loadkey = '111521'
mask_ecliptic = False

# load BOSS spectra north/south
if mask_ecliptic:
    alphas_norad = [np.load('../alphas_and_stds/alphas_boss_iris_2d_' + loadkey + '_10.npy'),
                    np.load('../alphas_and_stds/alphas_north_mask_ecliptic.npy'),
                    np.load('../alphas_and_stds/alphas_south_mask_ecliptic.npy')]
    alpha_stds_norad = [np.load('../alphas_and_stds/alpha_stds_boss_iris_2d_' + loadkey + '_10.npy'),
                        np.load('../alphas_and_stds/alpha_stds_north_mask_ecliptic.npy'),
                        np.load('../alphas_and_stds/alpha_stds_south_mask_ecliptic.npy')]
else:
    alphas_norad = [np.load('../alphas_and_stds/alphas_boss_iris_2d_' + loadkey + '_10.npy'),
                    np.load('../alphas_and_stds/alphas_boss_iris_2d_north_' + loadkey + '.npy'),
                    np.load('../alphas_and_stds/alphas_boss_iris_2d_south_' + loadkey + '.npy')]
    alpha_stds_norad = [np.load('../alphas_and_stds/alpha_stds_boss_iris_2d_' + loadkey + '_10.npy'),
                        np.load('../alphas_and_stds/alpha_stds_boss_iris_2d_north_' + loadkey + '.npy'),
                        np.load('../alphas_and_stds/alpha_stds_boss_iris_2d_south_' + loadkey + '.npy')]

# scaling factors
boss_fluxfactor = 1.38
sdss_fluxfactor = 1.38

# apply correction factor
correction_factors = [np.load('../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_10_smooth.npy'),
                      np.load('../alphas_and_stds/correction_factor_boss_iris_north_' + loadkey + '_smooth.npy'),
                      np.load('../alphas_and_stds/correction_factor_boss_iris_south_' + loadkey + '_smooth.npy')]
correction_factor_sdss = np.load('../alphas_and_stds/correction_factor_sdss_iris_smooth.npy')
lambdas_bin, binned_corrections, _ = generate_binned_alphas(correction_factors, 3 * [np.ones(len(correction_factors[0]))],
                                                  wavelength, boss=True, bin_width=bin_width)
_, binned_sdss_correction, _ = generate_binned_alphas([correction_factor_sdss], [np.ones(len(correction_factor_sdss))],
                                                      wavelength_sdss, boss=False, bin_width=bin_width)
binned_sdss_correction = binned_sdss_correction[0]

alphas_norad = [a / boss_fluxfactor * corr for a, corr in zip(alphas_norad, correction_factors)]
# bin
alphas_norad_bin_wav, alphas_norad_bin, _ = generate_binned_alphas(alphas_norad, alpha_stds_norad, wavelength, bin_offset=0, bin_width=bin_width)
alphas_boss_bin = alphas_norad_bin[0]
alphas_north_bin = alphas_norad_bin[1]
alphas_south_bin = alphas_norad_bin[2]

# bin sdss also
alphas_sdss = [np.load('../alphas_and_stds/alphas_sdss_iris_2d_102019.npy') * correction_factor_sdss / sdss_fluxfactor]  # just for the shape
sdss_bin_wav, sdss_bin, _ = generate_binned_alphas(alphas_sdss, [np.ones(alphas_sdss[0].shape)], wavelength_sdss, bin_offset=0, bin_width=bin_width)
sdss_bin = sdss_bin[0]

if bootstrap:
    # load bootstrap errors and bin them
    with open('../alphas_and_stds/bootstrap_binned_stds_boss' + loadkey + '.p', 'rb') as pf:
        bootstrap_binned_stds = pickle.load(pf)
    with open('../alphas_and_stds/bootstrap_binned_upper_boss' + loadkey + '.p', 'rb') as pf:
        bootstrap_binned_upper = pickle.load(pf)
    with open('../alphas_and_stds/bootstrap_binned_lower_boss' + loadkey + '.p', 'rb') as pf:
        bootstrap_binned_lower = pickle.load(pf)
    with open('../alphas_and_stds/bootstrap_binned_stds_sdss' + loadkey + '.p', 'rb') as pf:
        bootstrap_binned_stds_sdss = pickle.load(pf)
    if mask_ecliptic:
        bootstrap_binned_upper[0] = np.load('../alphas_and_stds/bootstrap_binned_upper_north_mask_ecliptic.npy')  # north
        bootstrap_binned_lower[0] = np.load('../alphas_and_stds/bootstrap_binned_lower_north_mask_ecliptic.npy')
        bootstrap_binned_upper[1] = np.load('../alphas_and_stds/bootstrap_binned_upper_south_mask_ecliptic.npy')  # south
        bootstrap_binned_lower[1] = np.load('../alphas_and_stds/bootstrap_binned_lower_south_mask_ecliptic.npy')

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

##########################
# comparisons of different rad. transfer params
##########################

default = np.load(load + 'defaults.npy') # t539, z=.02
no_1350 = np.load(load + 'no_1350.npy')
metal_008 = np.load(load + 'metallicity_008.npy')
b20 = np.load(load + 'b20.npy')
b30 = np.load(load + 'b30.npy')
b50 = np.load(load + 'b50.npy')
b60 = np.load(load + 'b60.npy')
tau005 = np.load(load + 'tau005.npy')
tau010 = np.load(load + 'tau010.npy')
zd = np.load(load + 'zd.npy')
zd_tau005 = np.load(load + 'zd_005.npy')
zd_tau010 = np.load(load + 'zd_010.npy')
test_converge = np.load(load + 'test_converge.npy')
alt_transform = np.load(load + 'alt_transform.npy')
dense_grid = np.load(load + 'dense_grid.npy')
denser_grid = np.load(load + 'denser_grid.npy')
grid_100s = np.load(load + 't5e9_grid_100_wd.npy')  # from ahab
grid_150s = np.load(load + 't5e9_grid_150_wd.npy')  # from ahab

mean_default = np.mean(default[(wavelength > 4200) & (wavelength < 5000)])
mean_no1350 = np.mean(no_1350[(wavelength > 4200) & (wavelength < 5000)])
mean_metal008 = np.mean(metal_008[(wavelength > 4200) & (wavelength < 5000)])
mean_b20 = np.mean(b20[(wavelength > 4200) & (wavelength < 5000)])
mean_b30 = np.mean(b30[(wavelength > 4200) & (wavelength < 5000)])
mean_b50 = np.mean(b50[(wavelength > 4200) & (wavelength < 5000)])
mean_b60 = np.mean(b60[(wavelength > 4200) & (wavelength < 5000)])
mean_tau005 = np.mean(tau005[(wavelength > 4200) & (wavelength < 5000)])
mean_tau010 = np.mean(tau010[(wavelength > 4200) & (wavelength < 5000)])
mean_zd = np.mean(zd[(wavelength > 4200) & (wavelength < 5000)])
mean_zd_tau005 = np.mean(zd_tau005[(wavelength > 4200) & (wavelength < 5000)])
mean_zd_tau010 = np.mean(zd_tau010[(wavelength > 4200) & (wavelength < 5000)])
mean_test = np.mean(test_converge[(wavelength > 4200) & (wavelength < 5000)])
mean_alt = np.mean(alt_transform[(wavelength > 4200) & (wavelength < 5000)])# this is with bigger wavelength range
if bootstrap:
    mean_north = np.average(alphas_north_bin[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)], weights=1/bootstrap_binned_stds_north[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)]**2)
    mean_south = np.average(alphas_south_bin[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)], weights=1/bootstrap_binned_stds_south[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)]**2)
else:
    mean_north = np.average(alphas_north_bin[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)])
    mean_south = np.average(alphas_south_bin[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)])

mean_choice = mean_south

# default parameters:
plt.plot(wavelength, default * mean_choice / mean_default, label='default')

# change coefficients to 350 vs. 1350 components:
# plt.plot(wavelength, no_1350 * mean_choice / mean_no1350, label='no 1350 component')

# change metallicity:
plt.plot(wavelength, metal_008 * mean_choice / mean_metal008, label='metal .008')

# check b values:
# plt.plot(wavelength, b20 * mean_choice / mean_b20, label='b = 20')
# plt.plot(wavelength, b30 * mean_choice / mean_b30, label='b = 30')
# plt.plot(wavelength, b50 * mean_choice / mean_b50, label='b = 50')
# plt.plot(wavelength, b60 * mean_choice / mean_b60, label='b = 60')

# check tau values for wd and zd:
# plt.plot(wavelength, tau005 * mean_choice / mean_tau005, label='tau = 0.05')
# plt.plot(wavelength, tau010 * mean_choice / mean_tau010, label='tau = 0.10')
# plt.plot(wavelength, zd * mean_choice / mean_zd, label='zd')
# plt.plot(wavelength, zd_tau005 * mean_choice / mean_zd_tau005, label='zd, tau=0.05')
# plt.plot(wavelength, zd_tau010 * mean_choice / mean_zd_tau010, label='zd, tau=0.10')

# check for convergence. these all have a different scaling bc of wider wav range
# plt.plot(wavelength, test_converge * mean_choice / mean_default, label='test convergence')
# plt.plot(wavelength, alt_transform * mean_choice / mean_default, label='alternate transformation')
# plt.plot(wavelength, dense_grid * mean_choice / mean_default, label='dense_grid')
# plt.plot(wavelength, denser_grid * mean_choice / mean_default, label='denser_grid')
# plt.plot(wavelength, grid_100s * mean_choice / mean_default, label='all 100s, from ahab')
# plt.plot(wavelength, grid_150s * mean_choice / mean_default, label='some 150s, from ahab')

plt.plot(alphas_norad_bin_wav, alphas_south_bin, color='blue', label='BOSS South', drawstyle='steps-mid')
plt.plot(alphas_norad_bin_wav, alphas_north_bin, color='red', label='BOSS North', drawstyle='steps-mid')
if bootstrap:
    plt.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_north,
                     bootstrap_binned_upper_north, linewidth=0.0, color='red', alpha=0.2, step='mid')
    plt.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_south,
                     bootstrap_binned_upper_south, linewidth=0.0, color='blue', alpha=0.2, step='mid')

plt.ylim(0, 0.4)
plt.legend()
plt.show()

# some calculations of effect size:
idx = np.argmin(np.abs(wavelength - 9000))
ratio_metal = (default[idx] / mean_default) / (metal_008[idx] / mean_metal008)
print("metal ratio:", ratio_metal)
ratio_dust = (default[idx] / mean_default) / (zd[idx] / mean_zd)
print("dust ratio:", ratio_dust)
ratio_tau = (default[idx] / mean_default) / (tau005[idx] / mean_tau005)
print("tau ratio:", ratio_tau)
ratio_b = (default[idx] / mean_default) / (b60[idx] / mean_b60)
print("b ratio:", ratio_b)


##########################
# PLOT BOSS WITH MULTIPLE DUST MODELS.
# SEE WHAT SCALING LOOKS GOOD.
##########################
show_full = False

# load_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/'
load_path = '/Users/blakechellew/Documents/DustProject/BrandtFiles/radiative/convergence_test/'
wd_t5e9_filename = 't5e9_grid_150_wd.npy'
zd_t5e9_filename = 't5e9_grid_150_zd.npy'
wd_t9e9_filename = 't9e9_grid_150_wd.npy'
zd_t9e9_filename = 't9e9_grid_150_zd.npy'
wd_cst_filename = 'cst_grid_150_wd.npy'
zd_cst_filename = 'cst_grid_150_zd.npy'
wd01_t5e9_alphas = np.load(load_path + wd_t5e9_filename)
zd01_t5e9_alphas = np.load(load_path + zd_t5e9_filename)
wd01_t9e9_alphas = np.load(load_path + wd_t9e9_filename)
zd01_t9e9_alphas = np.load(load_path + zd_t9e9_filename)
wd01_cst_alphas = np.load(load_path + wd_cst_filename)
zd01_cst_alphas = np.load(load_path + zd_cst_filename)
wd01_combo = np.load(load_path + wd_t5e9_filename) + cst_coeff_wd * np.load(load_path + wd_cst_filename)
zd01_combo = np.load(load_path + zd_t5e9_filename) + cst_coeff * np.load(load_path + zd_cst_filename)
# bin the alphas
_, wd01s_bin, _ = generate_binned_alphas([wd01_t5e9_alphas, wd01_t9e9_alphas, wd01_cst_alphas, wd01_combo],
                                         [np.ones(wd01_t5e9_alphas.shape) for i in range(4)], wavelength,
                                         bin_offset=0, bin_width=bin_width)
_, zd01s_bin, _ = generate_binned_alphas([zd01_t5e9_alphas, zd01_t9e9_alphas, zd01_cst_alphas, zd01_combo],
                                         [np.ones(wd01_t5e9_alphas.shape) for i in range(4)], wavelength,
                                         bin_offset=0, bin_width=bin_width)

# load alphas from BD12:
bd12_alphas_wd = np.loadtxt('/Users/blakechellew/Documents/DustProject/alphas_and_stds/bd12_fig3_green_052621.csv',
                             delimiter=",")
bd12_alphas_zd = np.loadtxt('/Users/blakechellew/Documents/DustProject/alphas_and_stds/bd12_fig3_blue_052621.csv',
                             delimiter=",")

# calculate scaling factor: (4200 to 5000 for now)
# and just avg value for now
# and use the north spectrum
mean_wds = [np.mean(wd01_bin[(lambdas_bin > 4200) & (lambdas_bin < 5000)]) for wd01_bin in wd01s_bin]
mean_zds = [np.mean(zd01_bin[(lambdas_bin > 4200) & (lambdas_bin < 5000)]) for zd01_bin in zd01s_bin]
if bootstrap:
    mean_boss = np.average(
        alphas_boss_bin[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)],
        weights=1/bootstrap_binned_stds_boss[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)]**2)
else:
    mean_boss = np.average(alphas_boss_bin[(alphas_norad_bin_wav > 4200) & (alphas_norad_bin_wav < 5000)])
scale_factors_wd_boss = [mean_boss / mean_wd for mean_wd in mean_wds]
scale_factors_zd_boss = [mean_boss / mean_zd for mean_zd in mean_zds]
scale_factors_wd_north = [mean_north / mean_wd for mean_wd in mean_wds]
scale_factors_zd_north = [mean_north / mean_zd for mean_zd in mean_zds]
scale_factors_wd_south = [mean_south / mean_wd for mean_wd in mean_wds]
scale_factors_zd_south = [mean_south / mean_zd for mean_zd in mean_zds]
south_scale = mean_north / mean_south
south_scale_bias = round(mean_north / mean_south * bias_factor, 1)
print("south scale:", south_scale)  # I initially chose 1.15 by eye

# 2-panel plot
if third_panel:
    plt.figure(figsize=(10, 8))
else:
    plt.figure(figsize=(10, 4))
# order: dark blue, yellow, cyan, light red, green
colors = ['#4477AA', '#CCBB44', '#66CCEE', '#EE6677', '#228833']

if third_panel:
    ax1 = plt.subplot(2, 2, 1)
else:
    ax1 = plt.subplot(1, 2, 1)
ax1.set_title('Correlation Spectra and BC03 Models')
# plt.plot(alphas_norad_bin_wav, alphas_boss_bin, 'orange', label='BOSS alphas (observed)', drawstyle='steps-mid')
ax1.plot(alphas_norad_bin_wav, bias_factor * south_scale * alphas_south_bin, colors[3], label=r'BOSS South  [$C=$' + str(south_scale_bias) + ']',
         drawstyle='steps-mid', linewidth=1.3)
ax1.plot(alphas_norad_bin_wav, bias_factor * alphas_north_bin, colors[0], label=r'BOSS North  [$C=$' + str(bias_factor) + ']',
         drawstyle='steps-mid', linewidth=1.3)
if show_full:
    ax1.plot(alphas_norad_bin_wav, alphas_boss_bin, 'k', label='BOSS overall',
             drawstyle='steps-mid')

# plot alphas from BD12 (manually found scaling factor for time reasons)
# the original scaling factor was 0.52, so now it's 0.52 * 1.15 = 0.6
# it's ZD, but different metallicity
# ax1.plot(bd12_alphas_zd[:, 0], 1.15*bd12_alphas_zd[:, 1], 'k', label='BD12 t5e9 (x 0.60)', drawstyle='steps-post')


scale_factors_wd = scale_factors_wd_north
scale_factors_zd = scale_factors_zd_north

# ax1.plot(lambdas_bin, wd01s_bin[0] * scale_factors_wd[0], 'k',
#          label='WD01 t5e9 (x ' + str(scale_factors_wd[0]) + ')', drawstyle='steps-mid')
"""
ax1.plot(lambdas_bin, zd01s_bin[0] * scale_factors_zd[0], 'blue',
         label='ZDA04 t5e9 (x ' + str(scale_factors_zd[0]) + ')', drawstyle='steps-mid')
ax1.plot(lambdas_bin, wd01s_bin[1] * scale_factors_wd[1], 'orange',
         label='WD01 t9e9 (x ' + str(scale_factors_wd[1]) + ')', drawstyle='steps-mid')
ax1.plot(lambdas_bin, zd01s_bin[1] * scale_factors_zd[1], 'yellow',
         label='ZDA04 t9e9 (x ' + str(scale_factors_zd[1]) + ')', drawstyle='steps-mid')
ax1.plot(lambdas_bin, wd01s_bin[2] * scale_factors_wd[2], 'pink',
         label='WD01 cst (x ' + str(scale_factors_wd[2]) + ')', drawstyle='steps-mid')
ax1.plot(lambdas_bin, zd01s_bin[2] * scale_factors_zd[2], 'cyan',
         label='ZDA04 cst (x ' + str(scale_factors_zd[2]) + ')', drawstyle='steps-mid')
ax1.plot(lambdas_bin, wd01s_bin[3] * scale_factors_wd[3], 'grey',
         label='WD01 combo (x ' + str(scale_factors_wd[3]) + ')', drawstyle='steps-mid')
ax1.plot(lambdas_bin, zd01s_bin[3] * scale_factors_zd[3], 'brown',
         label='ZDA04 combo (x ' + str(scale_factors_zd[3]) + ')', drawstyle='steps-mid')
"""

# plot some models, scaled to north and south:
ax1.plot(lambdas_bin, bias_factor * zd01s_bin[0] * scale_factors_zd_north[0], colors[1], linewidth=1.7,
         label=r't5e9  [$\times$' + str(round(bias_factor * scale_factors_zd_north[0], 2)) + ']', drawstyle='steps-mid')
# ax1.plot(lambdas_bin, zd01s_bin[0] * scale_factors_zd_south[0], colors[1],
#          label=r'ZDA04 t5e9 ($\times$ ' + str(scale_factors_zd_south[0]) + ')', drawstyle='steps-mid')
ax1.plot(lambdas_bin, bias_factor * zd01s_bin[3] * scale_factors_zd_north[3], colors[4], linewidth=1.7,
         label=r't5e9 + cst  [$\times$' + str(round(bias_factor * scale_factors_zd_north[3], 2)) + ']', drawstyle='steps-mid')


if bootstrap:
    ax1.fill_between(alphas_norad_bin_wav, bias_factor * bootstrap_binned_lower_north,
                     bias_factor * bootstrap_binned_upper_north, linewidth=0.0, color=colors[0], alpha=0.2,
                     step='mid')
    ax1.fill_between(alphas_norad_bin_wav, bias_factor * south_scale * bootstrap_binned_lower_south,
                     bias_factor * south_scale * bootstrap_binned_upper_south, linewidth=0.0, color=colors[3], alpha=0.2,
                     step='mid')
    if show_full:
        ax1.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_boss,
                         bootstrap_binned_upper_boss, linewidth=0.0, color='k', alpha=0.2,
                         step='mid')

ax1.set_ylim(0, 0.55)  # was .26 before bias factor
ax1.set_xlim(3800, 10000)
ax1.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.47, 0))
ax1.set_ylabel(r"$C\alpha_\lambda^{\prime}$ = $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
ax1.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")

ax1.xaxis.set_major_locator(MultipleLocator(1000))
ax1.xaxis.set_minor_locator(MultipleLocator(200))
ax1.yaxis.set_major_locator(MultipleLocator(0.1))
ax1.yaxis.set_minor_locator(MultipleLocator(0.02))

if third_panel:
    ax2 = plt.subplot(2, 2, 2)
else:
    ax2 = plt.subplot(1, 2, 2)

# define functions to get everything to line up
zd_fns = [interp1d(lambdas_bin, zd01_bin) for zd01_bin in zd01s_bin]
wd_fns = [interp1d(lambdas_bin, wd01_bin) for wd01_bin in wd01s_bin]
north_fn = interp1d(alphas_norad_bin_wav, alphas_north_bin)
south_fn = interp1d(alphas_norad_bin_wav, alphas_south_bin)
boss_fn = interp1d(alphas_norad_bin_wav, alphas_boss_bin)
sdss_fn = interp1d(sdss_bin_wav, sdss_bin)
wav_clipped = lambdas_bin[(alphas_norad_bin_wav > wav_clipped_min) & (alphas_norad_bin_wav < wav_clipped_max)]

if bootstrap:
    bootstrap_binned_upper_north_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_upper_north)
    bootstrap_binned_lower_north_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_lower_north)
    bootstrap_binned_upper_south_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_upper_south)
    bootstrap_binned_lower_south_fn = interp1d(alphas_norad_bin_wav, bootstrap_binned_lower_south)

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

if bootstrap:
    # find median signal to noise for boss vs. sdss
    boss_signal_to_noise = boss_fn(wav_clipped) / boss_errors_fn(wav_clipped)
    sdss_signal_to_noise = sdss_fn(wav_clipped) / sdss_errors_fn(wav_clipped)
    print("median signal to noise:")
    print("BOSS:", np.median(boss_signal_to_noise))
    print("SDSS:", np.median(sdss_signal_to_noise))

# TEST: was using ZD dust model for the last one

ax2.plot(alphas_norad_bin_wav, south_scale * bias_factor * (south_fn(alphas_norad_bin_wav) - zd_fns[3](alphas_norad_bin_wav) * scale_factors_zd_south[3]), colors[4],
         label=r'South $-$ (t5e9 + cst)  [$C=$' + str(south_scale_bias) + ']',
         drawstyle='steps-mid', linewidth=1.3)
ax2.plot(alphas_norad_bin_wav, south_scale * bias_factor * (south_fn(alphas_norad_bin_wav) - zd_fns[0](alphas_norad_bin_wav) * scale_factors_zd_south[0]), colors[1],
         label=r'South $-$ t5e9  [$C=$' + str(south_scale_bias) + ']', drawstyle='steps-mid', linewidth=1.3)
ax2.plot(alphas_norad_bin_wav, bias_factor * (north_fn(alphas_norad_bin_wav) - zd_fns[0](alphas_norad_bin_wav) * scale_factors_zd_north[0]), colors[0],
         label=r'North $-$ t5e9  [$C=$' + str(bias_factor) + ']', drawstyle='steps-mid', linewidth=1.7)
# TEST:
# ax2.plot(alphas_norad_bin_wav, south_scale * (south_fn(alphas_norad_bin_wav) - wd_fns[3](alphas_norad_bin_wav) * scale_factors_wd_south[3]), colors[3],
#          label='south - wd (combo) (x ' + str(south_scale) + ')', drawstyle='steps-mid')

if bootstrap:
    ax2.fill_between(alphas_norad_bin_wav,
                     bias_factor * (bootstrap_binned_lower_north_fn(alphas_norad_bin_wav) - zd_fns[0](alphas_norad_bin_wav) * scale_factors_zd_north[0]),
                     bias_factor * (bootstrap_binned_upper_north_fn(alphas_norad_bin_wav) - zd_fns[0](alphas_norad_bin_wav) * scale_factors_zd_north[0]),
                     linewidth=0.0, color=colors[2], alpha=0.2, step='mid')
    ax2.fill_between(alphas_norad_bin_wav,
                     south_scale * bias_factor * (bootstrap_binned_lower_south_fn(alphas_norad_bin_wav) - zd_fns[3](alphas_norad_bin_wav) * scale_factors_zd_south[3]),
                     south_scale * bias_factor * (bootstrap_binned_upper_south_fn(alphas_norad_bin_wav) - zd_fns[3](alphas_norad_bin_wav) * scale_factors_zd_south[3]),
                     linewidth=0.0, color=colors[4], alpha=0.2, step='mid')
    ax2.fill_between(alphas_norad_bin_wav,
                     south_scale * bias_factor * (bootstrap_binned_lower_south_fn(alphas_norad_bin_wav) - zd_fns[0](alphas_norad_bin_wav) * scale_factors_zd_south[0]),
                     south_scale * bias_factor * (bootstrap_binned_upper_south_fn(alphas_norad_bin_wav) - zd_fns[0](alphas_norad_bin_wav) * scale_factors_zd_south[0]),
                     linewidth=0.0, color=colors[1], alpha=0.2, step='mid')

ax2.set_ylabel(r"$C\left(\alpha^{\prime} - \alpha_{\mathrm{model}}^{\prime}\right)$")
ax2.set_title(r"Observations $-$ Model  (ERE Peak)")
ax2.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.47, 0))
ax2.hlines(0, 3650, 10200, color='grey')
ax2.set_ylim(-0.149, 0.17) # was -.07 to .085. TEST: scale to avoid overlapping the legend
# ax2.set_xlim(4000, 9000)
ax2.set_xlim(3800, 10000)
ax2.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")

ax2.xaxis.set_major_locator(MultipleLocator(1000))
ax2.xaxis.set_minor_locator(MultipleLocator(200))
ax2.yaxis.set_major_locator(MultipleLocator(0.05))
ax2.yaxis.set_minor_locator(MultipleLocator(0.01))

if third_panel:
    # add a 3rd panel
    ax3 = plt.subplot(2, 2, 3)

    ax3.set_title('Correlation Spectra and BC03 Models')
    # plt.plot(alphas_norad_bin_wav, alphas_boss_bin, 'orange', label='BOSS alphas (observed)', drawstyle='steps-mid')
    ax3.plot(alphas_norad_bin_wav, south_scale * alphas_south_bin, colors[3], label=r'BOSS South  [$\times$' + str(south_scale) + ']',
             drawstyle='steps-mid', linewidth=1.3)
    ax3.plot(alphas_norad_bin_wav, alphas_north_bin, colors[0], label='BOSS North',
             drawstyle='steps-mid', linewidth=1.3)

    # plot some models, scaled to north and south:
    ax3.plot(lambdas_bin, wd01s_bin[0] * scale_factors_wd_north[0], colors[1], linewidth=1.7,
             label=r't5e9  [$\times$' + str(scale_factors_wd_north[0]) + ']', drawstyle='steps-mid')
    # ax1.plot(lambdas_bin, zd01s_bin[0] * scale_factors_zd_south[0], colors[1],
    #          label=r'ZDA04 t5e9 ($\times$ ' + str(scale_factors_zd_south[0]) + ')', drawstyle='steps-mid')
    ax3.plot(lambdas_bin, wd01s_bin[3] * scale_factors_wd_north[3], colors[4], linewidth=1.7,
             label=r't5e9 + cst  [$\times$' + str(scale_factors_wd_north[3]) + ']', drawstyle='steps-mid')


    if bootstrap:
        ax3.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_north,
                         bootstrap_binned_upper_north, linewidth=0.0, color=colors[0], alpha=0.2,
                         step='mid')
        ax3.fill_between(alphas_norad_bin_wav, south_scale * bootstrap_binned_lower_south,
                         south_scale * bootstrap_binned_upper_south, linewidth=0.0, color=colors[3], alpha=0.2,
                      step='mid')
        if show_full:
            ax3.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_boss,
                             bootstrap_binned_upper_boss, linewidth=0.0, color='k', alpha=0.2,
                             step='mid')

    ax3.set_ylim(0, 0.26)
    ax3.set_xlim(3800, 10000)
    ax3.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.47, 0))
    ax3.set_ylabel(r"$\alpha_\lambda^{\prime}$ = $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
    ax3.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")

    ax3.xaxis.set_major_locator(MultipleLocator(1000))
    ax3.xaxis.set_minor_locator(MultipleLocator(200))
    ax3.yaxis.set_major_locator(MultipleLocator(0.1))
    ax3.yaxis.set_minor_locator(MultipleLocator(0.02))

# save and show
plt.tight_layout()
if save:
    plt.savefig('/Users/blakechellew/Documents/DustProject/paper_figures/ere_2plot_ecliptic_022322.pdf', bbox_inches='tight')
plt.show()

#######################

# bin in 10 A bins:
bin_10_wav, alphas_norad_bin, _ = generate_binned_alphas(alphas_norad, alpha_stds_norad, wavelength, bin_offset=0, bin_width=10)
alphas_north_bin10 = alphas_norad_bin[1]
alphas_south_bin10 = alphas_norad_bin[2]
_, zd01s_bin_10, _ = generate_binned_alphas([zd01_t5e9_alphas, zd01_t9e9_alphas, zd01_cst_alphas, zd01_combo],
                                         [np.ones(wd01_t5e9_alphas.shape) for i in range(4)], wavelength,
                                         bin_offset=0, bin_width=10)

# plot the 4000 A break up close, see what it looks like
plt.figure(figsize=(5, 7))
ax1 = plt.subplot(1, 1, 1)

colors = ['#4477AA', '#CCBB44', '#66CCEE', '#EE6677', '#228833']
ax1.set_title(r'$4000\,\mathrm{\AA}$ Break')
# [2] is south, [1] is north
ax1.plot(bin_10_wav, bias_factor * south_scale * alphas_north_bin10, colors[3], label=r'BOSS South  [$C=$' + str(south_scale_bias) + ']', drawstyle='steps-mid', linewidth=1.3, zorder=7)
ax1.plot(bin_10_wav, bias_factor * alphas_south_bin10, colors[0], label=r'BOSS North  [$C=$' + str(bias_factor) + ']', drawstyle='steps-mid', linewidth=1.3, zorder=7)


# plot some models, scaled to north and south:
ax1.plot(wavelength, bias_factor * zd01_t5e9_alphas * scale_factors_zd_north[0], colors[1], linewidth=1.9,
         label=r't5e9  [$\times$' + str(round(bias_factor * scale_factors_zd_north[0], 2)) + ']', drawstyle='steps-mid', zorder=5)
ax1.plot(wavelength, bias_factor * zd01_combo * scale_factors_zd_north[3], colors[4], linewidth=1.7,
         label=r't5e9 + cst  [$\times$' + str(round(bias_factor * scale_factors_zd_north[3], 2)) + ']', drawstyle='steps-mid', zorder=5)


# binning the models:
"""
# plot some models, scaled to north and south:
ax1.plot(bin_10_wav, bias_factor * zd01s_bin_10[0] * scale_factors_zd_north[0], colors[1], linewidth=1.9,
         label=r't5e9  [$\times$' + str(round(bias_factor * scale_factors_zd_north[0], 2)) + ']', drawstyle='steps-mid', zorder=5)
ax1.plot(bin_10_wav, bias_factor * zd01s_bin_10[3] * scale_factors_zd_north[3], colors[4], linewidth=1.7,
         label=r't5e9 + cst  [$\times$' + str(round(bias_factor * scale_factors_zd_north[3], 2)) + ']', drawstyle='steps-mid', zorder=5)
"""

# vertical line
ax1.axvline(x=4000, color='k', linewidth=1, linestyle='--')

"""
if bootstrap:
    plt.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_north,
                     bootstrap_binned_upper_north, linewidth=0.0, color=colors[0], alpha=0.2,
                     step='mid')
    plt.fill_between(alphas_norad_bin_wav, south_scale * bootstrap_binned_lower_south,
                     south_scale * bootstrap_binned_upper_south, linewidth=0.0, color=colors[3], alpha=0.2,
                     step='mid')
    if show_full:
        plt.fill_between(alphas_norad_bin_wav, bootstrap_binned_lower_boss,
                         bootstrap_binned_upper_boss, linewidth=0.0, color='k', alpha=0.2,
                         step='mid')
"""

ax1.set_ylim(0, 0.38)
ax1.set_xlim(3850, 4150)
ax1.legend(loc='lower right', facecolor='white', framealpha=1, edgecolor='white')  # bbox_to_anchor=(0.47, 0)
ax1.set_ylabel(r"$C\alpha_\lambda^{\prime}$ = $\lambda I_{\lambda}$ / $\nu I_\nu$ (100 $\mathrm{\mu}$m)")
ax1.set_xlabel(r"Wavelength ($\mathrm{\AA}$)")

ax1.xaxis.set_major_locator(MultipleLocator(100))
ax1.xaxis.set_minor_locator(MultipleLocator(20))
ax1.yaxis.set_major_locator(MultipleLocator(0.1))
ax1.yaxis.set_minor_locator(MultipleLocator(0.02))
plt.tight_layout()

if save_4000:
    plt.savefig('/Users/blakechellew/Documents/DustProject/paper_figures/break_4000_122021.pdf', bbox_inches='tight')
plt.show()


############################

# some ERE calculations:

# weighted i100 should be an average:
i100_weighted = np.load('/Users/blakechellew/Documents/DustProject/alphas_and_stds/avg_i100_boss_iris_' + loadkey + '_smooth.npy')[0]
i100_weighted_north = np.load('/Users/blakechellew/Documents/DustProject/alphas_and_stds/avg_i100_boss_iris_north_' + loadkey + '_smooth.npy')[0]
i100_weighted_south = np.load('/Users/blakechellew/Documents/DustProject/alphas_and_stds/avg_i100_boss_iris_south_' + loadkey + '_smooth.npy')[0]
print("weighted i100:")
print("north:", i100_weighted_north)
print("south:", i100_weighted_south)
# OVERRIDE: regular average
i100_weighted_north = 1.397
i100_weighted_south = 2.572

# integrate south - north:
integrand_south_minus_north = north_south_diff_fn(wav_clipped)
integrand_south_minus_north /= wav_clipped
integral_south_minus_north = 50 * np.nansum(integrand_south_minus_north)
print("Integral of south - north:", integral_south_minus_north)
if bootstrap:
    integral_south_minus_north_err = np.sqrt(np.sum(north_south_errors_fn(wav_clipped) ** 2))
    print("+/-", integral_south_minus_north_err)

if ere_bootstrap:
    # bootstrapping for ERE numbers:

    if ere_bootstrap:
        # really I should add the correction factor also, but won't make much difference
        bootstrap_alphas_boss = [np.load('../alphas_and_stds/bootstrap_alphas_boss_iris_2d_north_111521.npy') / boss_fluxfactor * correction_factors[1],
                                 np.load('../alphas_and_stds/bootstrap_alphas_boss_iris_2d_south_111521.npy') / boss_fluxfactor * correction_factors[2]]
        bootstrap_alpha_stds_boss = [np.load('../alphas_and_stds/bootstrap_alpha_stds_boss_iris_2d_north_111521.npy') / boss_fluxfactor * correction_factors[1],
                                     np.load('../alphas_and_stds/bootstrap_alpha_stds_boss_iris_2d_south_111521.npy') / boss_fluxfactor* correction_factors[2]]
        wav_boss_binned, bootstrap_binned_alphas_boss, _ = generate_binned_alphas(bootstrap_alphas_boss,
                                                                                  bootstrap_alpha_stds_boss, wavelength,
                                                                                  boss=True)
        bootstrap_binned_north, bootstrap_binned_south = bootstrap_binned_alphas_boss

    north_eres = []
    south_eres = []
    for i in range(1000):
        north_fn_boot = interp1d(wav_boss_binned, bootstrap_binned_north[i])
        south_fn_boot = interp1d(wav_boss_binned, bootstrap_binned_south[i])

        integrand_north_minus_model = bias_factor * (
                    north_fn_boot(wav_clipped) - zd_fns[0](wav_clipped) * scale_factors_zd_north[0])
        integrand_north_minus_model /= wav_clipped  # to get a unitless integral
        integral_north_minus_model = 50 * np.nansum(
            integrand_north_minus_model * i100_weighted_north * (3 * 10 ** 12))  # multiply by bin width 50 A
        integral_north_minus_model *= 10 ** -17  # unit conversion (to erg / s / cm^2 / sr)
        integrand_south_minus_combo = bias_factor * south_scale * (
                    south_fn_boot(wav_clipped) - zd_fns[3](wav_clipped) * scale_factors_zd_south[3])
        integrand_south_minus_combo /= wav_clipped
        integral_south_minus_combo = 50 * np.nansum(integrand_south_minus_combo * i100_weighted_south * (3 * 10 ** 12))
        integral_south_minus_combo *= 10 ** -17  # unit conversion (to erg / s / cm^2 / sr)
        north_eres.append(integral_north_minus_model)
        south_eres.append(integral_south_minus_combo)
    print("PERCENTILES")
    print(np.percentile(north_eres, (16, 50, 84)))
    print(np.percentile(south_eres, (16, 50, 84)))


# integrate north - model and south - model (assuming no model errors)
integrand_north_minus_model = bias_factor * (north_fn(wav_clipped) - zd_fns[0](wav_clipped) * scale_factors_zd_north[0])
integrand_north_minus_model /= wav_clipped  # to get a unitless integral
integral_north_minus_model = 50 * np.nansum(
    integrand_north_minus_model * i100_weighted_north * (3 * 10 ** 12))  # multiply by bin width 50 A
integral_north_minus_model *= 10 ** -17  # unit conversion (to erg / s / cm^2 / sr)
integrand_south_minus_combo = bias_factor * south_scale * (south_fn(wav_clipped) - zd_fns[3](wav_clipped) * scale_factors_zd_south[3])
integrand_south_minus_combo /= wav_clipped
integral_south_minus_combo = 50 * np.nansum(integrand_south_minus_combo * i100_weighted_south * (3 * 10 ** 12))
integral_south_minus_combo *= 10 ** -17  # unit conversion (to erg / s / cm^2 / sr)
integrand_south_minus_model = bias_factor * south_scale * (south_fn(wav_clipped) - zd_fns[0](wav_clipped) * scale_factors_zd_south[0])
integrand_south_minus_model /= wav_clipped
integral_south_minus_model = 50 * np.nansum(integrand_south_minus_model * i100_weighted_south * (3 * 10 ** 12))
integral_south_minus_model *= 10 ** -17  # unit conversion (to erg / s / cm^2 / sr)
integrand_south_minus_wd = bias_factor * south_scale * (south_fn(wav_clipped) - wd_fns[0](wav_clipped) * scale_factors_wd_south[0])
integrand_south_minus_wd /= wav_clipped
integral_south_minus_wd = 50 * np.nansum(integrand_south_minus_wd * i100_weighted_south * (3 * 10 ** 12))
integral_south_minus_wd *= 10 ** -17  # unit conversion (to erg / s / cm^2 / sr)
integrand_south_minus_wd_combo = bias_factor * south_scale * (south_fn(wav_clipped) - wd_fns[3](wav_clipped) * scale_factors_wd_south[3])
integrand_south_minus_wd_combo /= wav_clipped
integral_south_minus_wd_combo = 50 * np.nansum(integrand_south_minus_wd_combo * i100_weighted_south * (3 * 10 ** 12))
integral_south_minus_wd_combo *= 10 ** -17  # unit conversion (to erg / s / cm^2 / sr)
# multiplying by (nu*I_nu)_100 should give and actual value
# units here are MJy / sr (bc i100 is MJy / sr, and the dlambda units cancel with division by wav)
# convert:
# start: (MJy / sr) * Hz
# conversion factor: x10^-20 to get erg / s / cm^2 / sr

print("Integral of south - t5e9 zd:", integral_south_minus_model)
print("Integral of north - model:", integral_north_minus_model)
print("Integral of south - combo:", integral_south_minus_combo)
print("Integral of south - t5e9 wd:", integral_south_minus_wd)
print("Integral of south - wd combo:", integral_south_minus_wd_combo)
if bootstrap:
    integral_south_minus_model_err = np.sqrt(np.sum(south_errors_fn(wav_clipped) ** 2))
    integral_north_minus_model_err = np.sqrt(np.sum(north_errors_fn(wav_clipped) ** 2))
    print("error south: +/-", integral_south_minus_model_err)
    print("error north: +/-", integral_north_minus_model_err)

# total flux:
total_integrand_south = bias_factor * south_scale * south_fn(wav_clipped)
total_integrand_south /= wav_clipped
total_flux_south = np.nansum(50 * total_integrand_south * i100_weighted_south * (3 * 10 ** 12))
total_flux_south *= 10 ** -17  # unit conversion
total_integrand_north = bias_factor * north_fn(wav_clipped)
total_integrand_north /= wav_clipped
total_flux_north = np.nansum(50 * total_integrand_north * i100_weighted_north * (3 * 10 ** 12))
total_flux_north *= 10 ** -17  # unit conversion
print("Total flux south:", total_flux_south)
print("Total flux north:", total_flux_north)

# flux ratio:
south_ERE_ratio = integral_south_minus_model / (total_flux_south - integral_south_minus_model)
north_ERE_ratio = integral_north_minus_model / (total_flux_north - integral_north_minus_model)
south_ERE_ratio_combo = integral_south_minus_combo / (total_flux_south - integral_south_minus_combo)
south_ERE_ratio_wd = integral_south_minus_wd / (total_flux_south - integral_south_minus_wd)
south_ERE_ratio_wd_combo = integral_south_minus_wd_combo / (total_flux_south - integral_south_minus_wd_combo)
print("ratio ERE to scattered, south:", south_ERE_ratio)
print("ratio ERE to scattered, north:", north_ERE_ratio)
print("ratio ERE to scattered, south combo:", south_ERE_ratio_combo)
print("ratio ERE to scattered, south wd:", south_ERE_ratio_wd)
print("ratio ERE to scattered, south wd:", south_ERE_ratio_wd_combo)

# compare flux for CST vs t5e9
# NOTE: I'm adding them AFTER applying rad transfer
# also this is the whole range from 3600 to 10400
total_flux_t5e9 = np.nansum(50 * zd01s_bin[0] / alphas_norad_bin_wav)  # 50 not needed bc ratio
total_flux_cst = cst_coeff * np.nansum(50 * zd01s_bin[2] / alphas_norad_bin_wav)
t5e9_cst_ratio = total_flux_t5e9 / total_flux_cst
print("t5e9 cst ratio:", t5e9_cst_ratio)

if bootstrap:
    # comparing spectra (z-values):
    avg_north = np.mean(north_fn(wav_clipped))
    avg_south = np.mean(south_fn(wav_clipped))
    avg_north_err = np.median(north_errors_fn(wav_clipped))
    avg_south_err = np.median(south_errors_fn(wav_clipped))
    north_vs_south_err = np.sqrt(avg_north_err ** 2 + avg_south_err ** 2)
    z_north_vs_south = (avg_south - avg_north) / north_vs_south_err
    print("z north vs south:", z_north_vs_south)
    avg_sdss = np.mean(sdss_fn(wav_clipped))
    avg_boss = np.mean(boss_fn(wav_clipped))
    avg_sdss_err = np.median(bootstrap_binned_stds_sdss_fn(wav_clipped))
    sdss_vs_north_err = np.sqrt(avg_sdss_err ** 2 + avg_north_err ** 2)
    sdss_vs_south_err = np.sqrt(avg_sdss_err ** 2 + avg_south_err ** 2)
    z_sdss_vs_north = (avg_north - avg_sdss) / sdss_vs_north_err
    z_sdss_vs_south = (avg_south - avg_sdss) / sdss_vs_south_err
    print("z SDSS vs north:", z_sdss_vs_north)
    print("z SDSS vs south:", z_sdss_vs_south)
    print("avgs:")
    print("north:", avg_north)
    print("south:", avg_south)
    print("sdss:", avg_sdss)
    print("boss:", avg_boss)

# get ratio of ERE to total IR intensity
# Draine and Li 2007 dust model: (nu I_nu)_100 = (0.52 +/- .05) I_TIR
# from above: use i100_weighted_north and south
nu = 2.998 * 10**12
I_TIR_north = (nu * i100_weighted_north) / 0.52  # units of MJy / sr * Hz
I_TIR_south = (nu * i100_weighted_south) / 0.52
print("I_tir north:", I_TIR_north)
print("I_tir south:", I_TIR_south)
# the factor of 10 ** -17 converts the intensity back to units of (integral alpha/lambda) * (MJy/sr) * (Hz)
ratio_ERE_TIR_north = (integral_north_minus_model / (10 ** -17)) / I_TIR_north
ratio_ERE_TIR_south = (integral_south_minus_combo / (10 ** -17)) / I_TIR_south
print("ratio ERE to TIR north:", ratio_ERE_TIR_north)
print("ratio ERE to TIR south:", ratio_ERE_TIR_south)

# energy efficienct:
# gamma = UV / TIR
# efficiency = ERE / UV = ERE / (gamma * TIR)
gamma1_t5e9 = 0.510
gamma2_t5e9 = 0.416
gamma1_cst = 0.672
gamma2_cst = 0.547
gamma1_combo = cst_coeff * gamma1_cst + (1-cst_coeff) * gamma1_t5e9
gamma2_combo = cst_coeff * gamma2_cst + (1-cst_coeff) * gamma2_t5e9
print("combo gammas", gamma1_combo, gamma2_combo)
eff1_t5e9_north = ratio_ERE_TIR_north / gamma1_t5e9
eff1_t5e9_south = ratio_ERE_TIR_south / gamma1_t5e9
eff2_t5e9_north = ratio_ERE_TIR_north / gamma2_t5e9
eff2_t5e9_south = ratio_ERE_TIR_south / gamma2_t5e9
eff1_combo_north = ratio_ERE_TIR_north / gamma1_combo
eff1_combo_south = ratio_ERE_TIR_south / gamma1_combo
eff2_combo_north = ratio_ERE_TIR_north / gamma2_combo
eff2_combo_south = ratio_ERE_TIR_south / gamma2_combo
print("efficiencies: bound1, bound2")
print("t5e9 north:", eff1_t5e9_north, eff2_t5e9_north)
print("t5e9 south:", eff1_t5e9_south, eff2_t5e9_south)
print("combo north:", eff1_combo_north, eff2_combo_north)
print("combo south:", eff1_combo_south, eff2_combo_south)

# Szomoru calculations:
# A_v = 0.053 * (i100 / 1 MJy/sr)
Av_1 = 0.27
Av_2 = 0.64
i100_sz_1 = Av_1 / 0.053  # = 5.1 MJy/sr
i100_sz_2 = Av_2 / 0.053  # = 12.1 MJy/sr
# now convert to I_TIR:   nu * I_nu / 0.52
I_TIR_sz1 = nu * i100_sz_1 / 0.52 * 10**-17  # 10**-17 converts to same intensity units as Szomoru
I_TIR_sz2 = nu * i100_sz_2 / 0.52 * 10**-17
# now compare to Szomoru ERE value
I_ERE_sz = 1.2 * 10**-5  # (erg/s/cm^2/sr)
sz_ratio_1 = I_ERE_sz / I_TIR_sz1
sz_ratio_2 = I_ERE_sz / I_TIR_sz2
print("Szomoru I_ERE / I_TIR:")
print(sz_ratio_1)
print(sz_ratio_2)

# find the peaks and quartiles (don't care about the scaling)
# first the edges
def find_quartiles(integrand, lower_edge=5000, upper_edge=8000):
    #  = wav_clipped[np.nonzero(integrand < 0)]
    # lower_edge = 6500 - np.min(np.abs(wavs_subzero[wavs_subzero < 6500] - 6500))
    # upper_edge = 6500 + np.min(np.abs(wavs_subzero[wavs_subzero > 6500] - 6500))
    # print("lower edge:", lower_edge)
    # print("upper edge:", upper_edge)
    peak_integrand = integrand[(wav_clipped > lower_edge) & (wav_clipped < upper_edge)]
    peak_wavs = wav_clipped[(wav_clipped > lower_edge) & (wav_clipped < upper_edge)]
    total_integrand = np.nansum(peak_integrand)
    print("total ERE integrand:", total_integrand)  # intermediate step, units not adjusted
    partial_sum = 0
    quartile_1 = None
    quartile_2 = None
    quartile_3 = None
    for w, e in zip(peak_wavs, peak_integrand):
        if np.isnan(e):
            continue
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

print("Quartiles: south - t5e9")
find_quartiles(integrand_south_minus_model)
print("Quartiles: south - (t5e9 + cst)")
find_quartiles(integrand_south_minus_combo)
print("Quartiles: north - t5e9")
find_quartiles(integrand_north_minus_model)