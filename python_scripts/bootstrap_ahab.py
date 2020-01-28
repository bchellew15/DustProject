#read input from files with xxsig and yxsig sums
#bootstrap, using random samples of those

import numpy as np
import sys
import os

#command line args:
if len(sys.argv) != 3:
    print("Usage: bootstrap_ahab.py [save_key: sdss_1d, sdss_2d, boss_1d, boss_2d] [num_iter]")
    exit(0)
save_key = sys.argv[1] #should be same as input file
num_iter = int(sys.argv[2]) #number of bootstrap samples to add to the file

#load 
yxsig_bootstrap = np.load('yx_bootstrap_' + save_key + '.npy')
xxsig_bootstrap = np.load('xx_bootstrap_' + save_key + '.npy')

print(yxsig_bootstrap)
print(xxsig_bootstrap)

save_path1 = '../data/bootstrap_alphas_' + save_key + '.npy'
save_path2 = '../data/bootstrap_alpha_stds_' + save_key + '.npy'

#print("check sizes")
#print(yxsig_bootstrap.shape)
#print(xxsig_bootstrap.shape)

for i in range(num_iter):
    bootstrap_indices = np.random.choice(yxsig_bootstrap.shape[0], yxsig_bootstrap.shape[0])
    #bootstrap_indices = np.array(range(yxsig_bootstrap.shape[0]))

    sums1 = np.sum(yxsig_bootstrap[bootstrap_indices], axis=0)
    sums2 = np.sum(xxsig_bootstrap[bootstrap_indices], axis=0)

    alphas = np.divide(sums1, sums2)
    alpha_stds = np.sqrt(1/sums2)

    if os.path.isfile(save_path1) and os.path.isfile(save_path2): 
        bootstrap_alphas = np.load(save_path1)
        bootstrap_stds = np.load(save_path2)
        np.save(save_path1, np.append(bootstrap_alphas, alphas.reshape(1, len(alphas)), axis=0))
        np.save(save_path2, np.append(bootstrap_stds, alpha_stds.reshape(1, len(alphas)), axis=0))
    elif not os.path.isfile(save_path1) and not os.path.isfile(save_path2):
        np.save(save_path1, alphas.reshape(1, len(alphas)))
        np.save(save_path2, alpha_stds.reshape(1, len(alphas)))
    else:
        print("error: one file exists and the other doesn't")
        exit(0)
