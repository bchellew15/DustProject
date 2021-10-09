# TODO: incorporate uncertainties in the average

import os.path
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 4:
    print("Usage: run_bias_factor_code.py [recalculate: 0, 1] [loc: 0, 1, 2] [savekey]")
    exit(0)
recalculate = int(sys.argv[1])
location = int(sys.argv[2])
save_key = sys.argv[3]  # originally "biasfactor"

# min_thresholds = [0, 1]
min_thresholds = np.linspace(0, 5, 25)
save_path = '../alphas_and_stds/alphas_' + save_key + '.npy'
save_num_fibers = '../alphas_and_stds/num_fibers_' + save_key + '.npy'
# if that works try more:
# min_thresholds = [0, 0.33, 0.66, 1, 1.33, 1.66, 2, 2.33, 2.66, 3, 3.33, 3.66, 4]

if recalculate:
    # make sure there isn't a save file yet
    if os.path.isfile(save_path) or os.path.isfile(save_num_fibers):
        print("ERROR: file already exists.")
        exit(0)

    for min_thresh in min_thresholds:

        # run the code
        print("running:", min_thresh)
        os.system('python generate_alphas_biasfactor.py iris_2d 1 ' + save_key + ' ' + str(min_thresh) + ' 10 ' + str(location) + ' 0')

# make some plots:

# first just the alphas
boss_wav = np.load('../alphas_and_stds/wavelength_boss.npy')
threshold_alphas = np.load(save_path)
num_fibers = np.load(save_num_fibers)
N = threshold_alphas.shape[0]

# plot the first few:
# for i in range(4):
#     a = threshold_alphas[i]
#     plt.plot(boss_wav, a)
#     plt.ylim(0, 1)
#     plt.show()

# then plot the region 6600-6700 vs. threshold
averages = np.zeros(N)
for i in range(N):
    a = threshold_alphas[i]
    a_in_range = a[(boss_wav > 6600) & (boss_wav < 6700)]
    avg = np.mean(a_in_range)
    averages[i] = avg
plt.plot(min_thresholds, averages, label='average alphas')
plt.plot(min_thresholds, np.log(num_fibers) / 10, label='log(num fibers) / 10')  # scale to plot both at once
plt.legend()
plt.title("Average alphas vs. minimum excess 100 micron intensity")
plt.show()


