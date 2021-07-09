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

_, bootstrap_binned_alphas_thresh_1d, _ = generate_binned_alphas(bootstrap_alphas_thresh_1d, bootstrap_alpha_stds_thresh_1d, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_thresh_2d, _ = generate_binned_alphas(bootstrap_alphas_thresh_2d, bootstrap_alpha_stds_thresh_2d, wavelength_boss, boss=True)

# bootstrap percentiles
bootstrap_binned_lower_thresh_1d = [np.nanpercentile(b, 16, axis=0) for b in bootstrap_binned_alphas_thresh_1d]  # 68 percent confidence interval
bootstrap_binned_upper_thresh_1d = [np.nanpercentile(b, 84, axis=0) for b in bootstrap_binned_alphas_thresh_1d]
bootstrap_binned_stds_thresh_1d = [(bootstrap_binned_upper_thresh_1d[i] - bootstrap_binned_lower_thresh_1d[i]) / 2 for i in range(len(bootstrap_binned_lower_thresh_1d))]
bootstrap_binned_lower_thresh_2d = [np.nanpercentile(b, 16, axis=0) for b in bootstrap_binned_alphas_thresh_2d]  # 68 percent confidence interval
bootstrap_binned_upper_thresh_2d = [np.nanpercentile(b, 84, axis=0) for b in bootstrap_binned_alphas_thresh_2d]
bootstrap_binned_stds_thresh_2d = [(bootstrap_binned_upper_thresh_2d[i] - bootstrap_binned_lower_thresh_2d[i]) / 2 for i in range(len(bootstrap_binned_lower_thresh_2d))]

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

# BOSS:
bootstrap_lower_boss = [np.nanpercentile(b, 16, axis=0) for b in bootstrap_alphas_boss]
bootstrap_upper_boss = [np.nanpercentile(b, 84, axis=0) for b in bootstrap_alphas_boss]
bootstrap_stds_boss = [(bootstrap_upper_boss[i] - bootstrap_lower_boss[i]) / 2 for i in range(len(bootstrap_lower_boss))]
for i in range(len(bootstrap_stds_boss)):
    assert bootstrap_stds_boss[i].ndim == 1

print("check 1")

# SDSS:
bootstrap_lower_sdss = [np.nanpercentile(b, 16, axis=0) for b in bootstrap_alphas_sdss]
bootstrap_upper_sdss = [np.nanpercentile(b, 84, axis=0) for b in bootstrap_alphas_sdss]
bootstrap_stds_sdss = [(bootstrap_upper_sdss[i] - bootstrap_lower_sdss[i]) / 2 for i in range(len(bootstrap_lower_sdss))]
for i in range(len(bootstrap_stds_boss)):
    assert bootstrap_stds_boss[i].ndim == 1

print("check 2")

# binning:
_, bootstrap_binned_alphas_boss, _ = generate_binned_alphas(bootstrap_alphas_boss, bootstrap_alpha_stds_boss, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_sdss, _ = generate_binned_alphas(bootstrap_alphas_sdss, bootstrap_alpha_stds_sdss, wavelength_sdss, boss=False)
# bin the boss alphas according to sdss wavelengths
_, bootstrap_binned_alphas_boss_to_sdss, _ = generate_binned_alphas([bootstrap_alphas_boss[-1]], [bootstrap_alpha_stds_boss[-1]], wavelength_sdss, wavelength_boss, boss=True)

# percentiles for binned plots:

# BOSS
bootstrap_binned_lower_boss = [np.nanpercentile(b, 16, axis=0) for b in bootstrap_binned_alphas_boss] #68 percent confidence interval
bootstrap_binned_upper_boss = [np.nanpercentile(b, 84, axis=0) for b in bootstrap_binned_alphas_boss]
bootstrap_binned_stds_boss = [(bootstrap_binned_upper_boss[i] - bootstrap_binned_lower_boss[i])/2 for i in range(len(bootstrap_binned_lower_boss))]

# SDSS
bootstrap_binned_lower_sdss = [np.nanpercentile(b, 16, axis=0) for b in bootstrap_binned_alphas_sdss] #68 percent confidence interval
bootstrap_binned_upper_sdss = [np.nanpercentile(b, 84, axis=0) for b in bootstrap_binned_alphas_sdss]
bootstrap_binned_stds_sdss = [(bootstrap_binned_upper_sdss[i] - bootstrap_binned_lower_sdss[i])/2 for i in range(len(bootstrap_binned_lower_sdss))]

# boss to sdss:
bootstrap_binned_lower_boss_to_sdss = np.nanpercentile(bootstrap_binned_alphas_boss_to_sdss[0], 16, axis=0)  #68 percent confidence interval
bootstrap_binned_upper_boss_to_sdss = np.nanpercentile(bootstrap_binned_alphas_boss_to_sdss[0], 84, axis=0)
bootstrap_binned_stds_boss_to_sdss = (bootstrap_binned_upper_boss_to_sdss - bootstrap_binned_lower_boss_to_sdss)/2
assert bootstrap_binned_stds_boss_to_sdss.ndim == 1

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

###############################################