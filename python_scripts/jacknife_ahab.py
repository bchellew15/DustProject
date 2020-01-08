#read input from files with xxsig and yxsig sums
#jacknife to leave out each plate one-by-one

import numpy as np
import sys
import os

#command line args:
if len(sys.argv) != 2:
    print("Usage: jacknife_ahab.py [save_key: sdss_1d, sdss_2d, boss_1d, boss_2d]")
    exit(0)
save_key = sys.argv[1] #should be same as input file

#load 
yxsig_jacknife = np.load('yx_bootstrap_' + save_key + '.npy')
xxsig_jacknife = np.load('xx_bootstrap_' + save_key + '.npy')

save_path1 = '../data/jacknife_alphas_' + save_key + '.npy'
save_path2 = '../data/jacknife_alpha_stds_' + save_key + '.npy'

#check sizes:
#print("check sizes")
#print(yxsig_bootstrap.shape)
#print(xxsig_bootstrap.shape)

for i in range(yxsig_jacknife.shape[0]):
    jacknife_indices = [j for j in np.arange(yxsig_jacknife.shape[0]) if j != i]
    #all the plates except the one at index i

    sums1 = np.sum(yxsig_jacknife[jacknife_indices], axis=0)
    sums2 = np.sum(xxsig_jacknife[jacknife_indices], axis=0)

    alphas = np.divide(sums1, sums2)
    alpha_stds = np.sqrt(1/sums2)


    if os.path.isfile(save_path1) and os.path.isfile(save_path2): 
        jacknife_alphas = np.load(save_path1)
        jacknife_stds = np.load(save_path2)
        np.save(save_path1, np.append(jacknife_alphas, alphas.reshape(1, len(alphas)), axis=0))
        np.save(save_path2, np.append(jacknife_stds, alpha_stds.reshape(1, len(alphas)), axis=0))
    elif not os.path.isfile(save_path1) and not os.path.isfile(save_path2):
        np.save(save_path1, alphas.reshape(1, len(alphas)))
        np.save(save_path2, alpha_stds.reshape(1, len(alphas)))
    else:
        print("error: one file exists and the other doesn't")
        exit(0)
