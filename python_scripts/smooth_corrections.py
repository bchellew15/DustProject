# quick script to smooth out the correction factors

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from astropy.io import fits

wav_boss = np.load('../alphas_and_stds/wavelength_boss.npy')
hdulist = fits.open( '../data/SDSS_allskyspec.fits')
wav_sdss = np.array(hdulist[1].data)
loadkey = '111521'

"""
# Smooth the i100
i100_files = [
    '../alphas_and_stds/avg_i100_boss_iris.npy',
    '../alphas_and_stds/avg_i100_boss_iris_north.npy',
    '../alphas_and_stds/avg_i100_boss_iris_south.npy',
]
i100s = [np.load(f) for f in i100_files]

for i, i100 in enumerate(i100s):
    i100_smooth = ndimage.gaussian_filter(i100, 30)
    i100_mean_smooth = np.copy(i100_smooth)
    i100_mean_smooth[:] = np.mean(i100_mean_smooth[2890:3010])  # wavelengths 6900 to 7090
    plt.plot(wav_boss, i100)
    plt.plot(wav_boss, i100_smooth)
    plt.plot(wav_boss, i100_mean_smooth)
    plt.show()
    np.save(i100_files[i].split('.npy')[0] + '_smooth.npy', i100_mean_smooth)

exit(0)
"""

######################
# Now smooth the correction factors
##########################

# BOSS
files_boss = [
    '../alphas_and_stds/correction_factor_boss_iris_north_' + loadkey + '.npy',
    '../alphas_and_stds/correction_factor_boss_iris_south_' + loadkey + '.npy',
    '../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_10.npy',
    '../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_15.npy',
    '../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_20.npy',
    '../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_25.npy',
    '../alphas_and_stds/correction_factor_boss_iris_2d_' + loadkey + '_30.npy'
]
correction_factors_boss = [np.load(file) for file in files_boss]

# SDSS
files_sdss = [
    '../alphas_and_stds/correction_factor_sdss_sfd.npy',
    '../alphas_and_stds/correction_factor_sdss_iris.npy'
]
correction_factors_sdss = [np.load(file) for file in files_sdss]

for i, corr in enumerate(correction_factors_boss):
    corr_smooth = ndimage.gaussian_filter(corr, 30)
    plt.plot(wav_boss, corr)
    plt.plot(wav_boss, corr_smooth)
    plt.show()
    np.save(files_boss[i].split('.npy')[0] + '_smooth.npy', corr_smooth)

"""
for i, corr in enumerate(correction_factors_sdss):
    corr_smooth = ndimage.gaussian_filter(corr, 30)
    plt.plot(wav_sdss, corr)
    plt.plot(wav_sdss, corr_smooth)
    plt.show()
    np.save(files_sdss[i].split('.npy')[0] + '_smooth.npy', corr_smooth)
"""

