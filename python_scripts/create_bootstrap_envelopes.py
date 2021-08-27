# run this on ahab.
# create the bootstrap envolopes here so the plotting code can run faster.
# need to divide by flux conversion factor somewhere else.

import numpy as np
from generate_plots import generate_binned_alphas
from astropy.io import fits
import pickle

loadkey = '012720'
alpha_direc_boot = '../data/'
hdulist_direc = '../data/'

# load wavelengths
wavelength_boss = np.load('../alphas_and_stds/wavelength_boss.npy')
hdulist = fits.open(hdulist_direc + 'SDSS_allskyspec.fits')
wavelength_sdss = np.array(hdulist[1].data)

#############################################################

# now for the threshold plots

bootstrap_alphas_thresh_1d = [np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_10.npy'),
                              np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_15.npy'),
                              np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_20.npy'),
                              np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_25.npy'),
                              np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_1d_' + loadkey + '_30.npy')]
bootstrap_alpha_stds_thresh_1d = [np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_10.npy'),
                                  np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_15.npy'),
                                  np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_20.npy'),
                                  np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_25.npy'),
                                  np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_1d_' + loadkey + '_30.npy')]
bootstrap_alphas_thresh_2d = [np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_10.npy'),
                              np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_15.npy'),
                              np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_20.npy'),
                              np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_25.npy'),
                              np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_30.npy')]
bootstrap_alpha_stds_thresh_2d = [np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_10.npy'),
                                  np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_15.npy'),
                                  np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_20.npy'),
                                  np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_25.npy'),
                                  np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_30.npy')]
print("loaded threshold alphas")

_, bootstrap_binned_alphas_thresh_1d, _ = generate_binned_alphas(bootstrap_alphas_thresh_1d, bootstrap_alpha_stds_thresh_1d, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_thresh_2d, _ = generate_binned_alphas(bootstrap_alphas_thresh_2d, bootstrap_alpha_stds_thresh_2d, wavelength_boss, boss=True)
print("done binning")

# bootstrap percentiles
bootstrap_binned_lower_thresh_1d = [np.percentile(b, 16, axis=0) for b in bootstrap_binned_alphas_thresh_1d]  # 68 percent confidence interval
bootstrap_binned_upper_thresh_1d = [np.percentile(b, 84, axis=0) for b in bootstrap_binned_alphas_thresh_1d]
bootstrap_binned_stds_thresh_1d = [(bootstrap_binned_upper_thresh_1d[i] - bootstrap_binned_lower_thresh_1d[i]) / 2 for i in range(len(bootstrap_binned_lower_thresh_1d))]
bootstrap_binned_lower_thresh_2d = [np.percentile(b, 16, axis=0) for b in bootstrap_binned_alphas_thresh_2d]  # 68 percent confidence interval
bootstrap_binned_upper_thresh_2d = [np.percentile(b, 84, axis=0) for b in bootstrap_binned_alphas_thresh_2d]
bootstrap_binned_stds_thresh_2d = [(bootstrap_binned_upper_thresh_2d[i] - bootstrap_binned_lower_thresh_2d[i]) / 2 for i in range(len(bootstrap_binned_lower_thresh_2d))]
print("done with percentiles")

# save
with open(alpha_direc_boot + 'bootstrap_binned_lower_thresh_1d' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_lower_thresh_1d, fp)
with open(alpha_direc_boot + 'bootstrap_binned_upper_thresh_1d' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_upper_thresh_1d, fp)
with open(alpha_direc_boot + 'bootstrap_binned_stds_thresh_1d' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_stds_thresh_1d, fp)
with open(alpha_direc_boot + 'bootstrap_binned_lower_thresh_2d' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_lower_thresh_2d, fp)
with open(alpha_direc_boot + 'bootstrap_binned_upper_thresh_2d' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_upper_thresh_2d, fp)
with open(alpha_direc_boot + 'bootstrap_binned_stds_thresh_2d' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_stds_thresh_2d, fp)
print("done saving")

#########################################################

# load bootstrap stuff:
bootstrap_alphas_boss = [np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_north_' + loadkey + '.npy'),
                         np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_south_' + loadkey + '.npy'),
                         np.load(alpha_direc_boot + 'bootstrap_alphas_boss_iris_2d_' + loadkey + '_10.npy')]
bootstrap_alpha_stds_boss = [np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_north_' + loadkey + '.npy'),
                             np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_south_' + loadkey + '.npy'),
                             np.load(alpha_direc_boot + 'bootstrap_alpha_stds_boss_iris_2d_' + loadkey + '_10.npy')]

bootstrap_alphas_sdss = [np.load(alpha_direc_boot + 'bootstrap_alphas_sdss_1d_' + loadkey + '.npy'),
                         np.load(alpha_direc_boot + 'bootstrap_alphas_sdss_2d_' + loadkey + '.npy'),
                         np.load(alpha_direc_boot + 'bootstrap_alphas_sdss_iris_2d_' + loadkey + '.npy'),
                         np.load(alpha_direc_boot + 'bootstrap_alphas_sdss_iris_1d_' + loadkey + '.npy')]
bootstrap_alpha_stds_sdss = [np.load(alpha_direc_boot + 'bootstrap_alpha_stds_sdss_1d_' + loadkey + '.npy'),
                             np.load(alpha_direc_boot + 'bootstrap_alpha_stds_sdss_2d_' + loadkey + '.npy'),
                             np.load(alpha_direc_boot + 'bootstrap_alpha_stds_sdss_iris_2d_' + loadkey + '.npy'),
                             np.load(alpha_direc_boot + 'bootstrap_alpha_stds_sdss_iris_1d_' + loadkey + '.npy')]
print("loaded main alphas")

# BOSS:
bootstrap_lower_boss = [np.percentile(b, 16, axis=0) for b in bootstrap_alphas_boss]
bootstrap_upper_boss = [np.percentile(b, 84, axis=0) for b in bootstrap_alphas_boss]
bootstrap_stds_boss = [(bootstrap_upper_boss[i] - bootstrap_lower_boss[i]) / 2 for i in range(len(bootstrap_lower_boss))]
for i in range(len(bootstrap_stds_boss)):
    assert bootstrap_stds_boss[i].ndim == 1

print("check 1")

# SDSS:
bootstrap_lower_sdss = [np.percentile(b, 16, axis=0) for b in bootstrap_alphas_sdss]
bootstrap_upper_sdss = [np.percentile(b, 84, axis=0) for b in bootstrap_alphas_sdss]
bootstrap_stds_sdss = [(bootstrap_upper_sdss[i] - bootstrap_lower_sdss[i]) / 2 for i in range(len(bootstrap_lower_sdss))]
for i in range(len(bootstrap_stds_boss)):
    assert bootstrap_stds_boss[i].ndim == 1

print("check 2")

# binning:
wav_boss_binned, bootstrap_binned_alphas_boss, _ = generate_binned_alphas(bootstrap_alphas_boss, bootstrap_alpha_stds_boss, wavelength_boss, boss=True)
wav_sdss_binned, bootstrap_binned_alphas_sdss, _ = generate_binned_alphas(bootstrap_alphas_sdss, bootstrap_alpha_stds_sdss, wavelength_sdss, boss=False)
# bin the boss alphas according to sdss wavelengths
_, bootstrap_binned_alphas_boss_to_sdss, _ = generate_binned_alphas([bootstrap_alphas_boss[-1]], [bootstrap_alpha_stds_boss[-1]], wavelength_sdss, wavelength_boss, boss=True)

# percentiles for binned plots:

# BOSS
bootstrap_binned_lower_boss = [np.percentile(b, 16, axis=0) for b in bootstrap_binned_alphas_boss] #68 percent confidence interval
bootstrap_binned_upper_boss = [np.percentile(b, 84, axis=0) for b in bootstrap_binned_alphas_boss]
bootstrap_binned_stds_boss = [(bootstrap_binned_upper_boss[i] - bootstrap_binned_lower_boss[i])/2 for i in range(len(bootstrap_binned_lower_boss))]

# SDSS
bootstrap_binned_lower_sdss = [np.percentile(b, 16, axis=0) for b in bootstrap_binned_alphas_sdss] #68 percent confidence interval
bootstrap_binned_upper_sdss = [np.percentile(b, 84, axis=0) for b in bootstrap_binned_alphas_sdss]
bootstrap_binned_stds_sdss = [(bootstrap_binned_upper_sdss[i] - bootstrap_binned_lower_sdss[i])/2 for i in range(len(bootstrap_binned_lower_sdss))]

# boss to sdss:
bootstrap_binned_lower_boss_to_sdss = np.percentile(bootstrap_binned_alphas_boss_to_sdss[0], 16, axis=0)  #68 percent confidence interval
bootstrap_binned_upper_boss_to_sdss = np.percentile(bootstrap_binned_alphas_boss_to_sdss[0], 84, axis=0)
bootstrap_binned_stds_boss_to_sdss = (bootstrap_binned_upper_boss_to_sdss - bootstrap_binned_lower_boss_to_sdss)/2
assert bootstrap_binned_stds_boss_to_sdss.ndim == 1
print("done with binning")

# save unbinned
with open(alpha_direc_boot + 'bootstrap_lower_boss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_lower_boss, fp)
with open(alpha_direc_boot + 'bootstrap_upper_boss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_upper_boss, fp)
with open(alpha_direc_boot + 'bootstrap_stds_boss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_stds_boss, fp)
with open(alpha_direc_boot + 'bootstrap_lower_sdss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_lower_sdss, fp)
with open(alpha_direc_boot + 'bootstrap_upper_sdss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_upper_sdss, fp)
with open(alpha_direc_boot + 'bootstrap_stds_sdss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_stds_sdss, fp)
# save binned
with open(alpha_direc_boot + 'bootstrap_binned_lower_boss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_lower_boss, fp)
with open(alpha_direc_boot + 'bootstrap_binned_upper_boss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_upper_boss, fp)
with open(alpha_direc_boot + 'bootstrap_binned_stds_boss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_stds_boss, fp)
with open(alpha_direc_boot + 'bootstrap_binned_lower_sdss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_lower_sdss, fp)
with open(alpha_direc_boot + 'bootstrap_binned_upper_sdss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_upper_sdss, fp)
with open(alpha_direc_boot + 'bootstrap_binned_stds_sdss' + loadkey + '.p', 'wb') as fp:
    pickle.dump(bootstrap_binned_stds_sdss, fp)
# save boss to sdss
np.save(alpha_direc_boot + 'bootstrap_binned_lower_boss_to_sdss' + loadkey + '.npy', bootstrap_binned_lower_boss_to_sdss)
np.save(alpha_direc_boot + 'bootstrap_binned_upper_boss_to_sdss' + loadkey + '.npy', bootstrap_binned_upper_boss_to_sdss)
np.save(alpha_direc_boot + 'bootstrap_binned_stds_boss_to_sdss' + loadkey + '.npy', bootstrap_binned_stds_boss_to_sdss)
print("done saving")

###############################################

# misc calculations:
# averages and percentiles for some of these
select_north = bootstrap_binned_alphas_boss[0][:, (wav_boss_binned > 5000) & (wav_boss_binned < 8000)]
select_south = bootstrap_binned_alphas_boss[1][:, (wav_boss_binned > 5000) & (wav_boss_binned < 8000)]
select_boss = bootstrap_binned_alphas_boss[2][:, (wav_boss_binned > 5000) & (wav_boss_binned < 8000)]
select_sdss = bootstrap_binned_alphas_sdss[2][:, (wav_sdss_binned > 5000) & (wav_sdss_binned < 8000)]
means_north = np.mean(select_north, axis=1).flatten()
means_south = np.mean(select_south, axis=1).flatten()
means_sdss = np.mean(select_sdss, axis=1).flatten()
means_boss = np.mean(select_boss, axis=1).flatten()
means_north_lower = np.percentile(means_north, 16)
means_north_upper = np.percentile(means_north, 84)
means_north_err = (means_north_upper - means_north_lower) / 2
means_south_lower = np.percentile(means_south, 16)
means_south_upper = np.percentile(means_south, 84)
means_south_err = (means_south_upper - means_south_lower) / 2
means_sdss_lower = np.percentile(means_sdss, 16)
means_sdss_upper = np.percentile(means_sdss, 84)
means_sdss_err = (means_sdss_upper - means_sdss_lower) / 2
means_boss_lower = np.percentile(means_boss, 16)
means_boss_upper = np.percentile(means_boss, 84)
means_boss_err = (means_boss_upper - means_boss_lower) / 2
print("means north err:", means_north_err)
print("means south err:", means_south_err)
print("means sdss err:", means_sdss_err)
print("means boss err:", means_boss_err)
# 50th percentiles:
means_south_50th = np.percentile(means_south, 50)
means_north_50th = np.percentile(means_north, 50)
print("50th percentiles:", means_south_50th, means_north_50th)
# try from 50th to inside:
means_south_diff_50 = means_south_50th - means_south_lower
means_north_diff_50 = means_north_upper - means_north_50th
print("means from 50th: (south, north)", means_south_diff_50, means_north_diff_50)
# compare to average bootstrap error:
bootstrap_stds_north = bootstrap_binned_stds_boss[0][(wav_boss_binned > 5000) & (wav_boss_binned < 8000)]
bootstrap_stds_south = bootstrap_binned_stds_boss[1][(wav_boss_binned > 5000) & (wav_boss_binned < 8000)]
bootstrap_mean_north = np.mean(bootstrap_stds_north)
bootstrap_mean_south = np.mean(bootstrap_stds_south)
print("avg bootstrap stds:")
print("north:", bootstrap_mean_north)
print("south:", bootstrap_mean_south)


