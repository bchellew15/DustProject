#loop over other scripts
#goal is to generate all the alphas and bootstrap samples needed for the paper
#input a unique key (a date, for example) that will be used for the filenames

import os
import sys

#savekey is a unique key that will be used for this set of filenames
if len(sys.argv) != 2:
    print("Usage: script_looper.py [savekey]")
    exit(0)
savekey = sys.argv[1]

num_samples = 10000 #10,000 bootstrap samples
thresholds = ['10', '15', '20', '25', '30']
threshold_labels = ['10', '15', '20', '25', '30'] #in case there is a decimal point

#generate the alphas:
"""
os.system('python generate_alphas_ahab.py 1d 0 sdss_1d_' + savekey + ' 10 0 0')
os.system('python generate_alphas_ahab.py iris_1d 0 sdss_iris_1d_' + savekey + ' 10 0 0')
os.system('python generate_alphas_ahab.py 2d 0 sdss_2d_' + savekey + ' 10 0 0')
os.system('python generate_alphas_ahab.py iris_2d 0 sdss_iris_2d_' + savekey + ' 10 0 0')
"""
# for i in range(len(thresholds)):
    # os.system('python generate_alphas_ahab.py iris_1d 1 boss_iris_1d_' + savekey + '_' + threshold_labels[i] + ' ' + thresholds[i] + ' 0 0 0')
    # os.system('python generate_alphas_ahab.py iris_2d 1 boss_iris_2d_' + savekey + '_' + threshold_labels[i] + ' ' + thresholds[i] + ' 0 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 boss_iris_2d_north_' + savekey + ' 10 1 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 boss_iris_2d_south_' + savekey + ' 10 2 0')

#intermediate bootstrap step for all of these:
"""
os.system('python generate_alphas_ahab.py 1d 0 sdss_1d_' + savekey + ' 10 0 1')
os.system('python generate_alphas_ahab.py iris_1d 0 sdss_iris_1d_' + savekey + ' 10 0 1')
os.system('python generate_alphas_ahab.py 2d 0 sdss_2d_' + savekey + ' 10 0 1')
os.system('python generate_alphas_ahab.py iris_2d 0 sdss_iris_2d_' + savekey + ' 10 0 1')
"""
for i in range(len(thresholds)):
    # os.system('python generate_alphas_ahab.py iris_1d 1 boss_iris_1d_' + savekey + '_' + threshold_labels[i] + ' ' + thresholds[i] + ' 0 1 0')
    os.system('python generate_alphas_ahab.py iris_2d 1 boss_iris_2d_' + savekey + '_' + threshold_labels[i] + ' ' + thresholds[i] + ' 0 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 boss_iris_2d_north_' + savekey + ' 10 1 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 boss_iris_2d_south_' + savekey + ' 10 2 1 0')


#now calculate bootstrap samples:
# os.system('python bootstrap_ahab.py sdss_1d_' + savekey + ' 10000')
# os.system('python bootstrap_ahab.py sdss_iris_1d_' + savekey + ' 10000')
# os.system('python bootstrap_ahab.py sdss_2d_' + savekey + ' 10000')
# os.system('python bootstrap_ahab.py sdss_iris_2d_' + savekey + ' 10000')
for i in range(len(thresholds)):
    # os.system('python bootstrap_ahab.py boss_iris_1d_' + savekey + '_' + threshold_labels[i] + ' 10000')
    os.system('python bootstrap_ahab.py boss_iris_2d_' + savekey + '_' + threshold_labels[i] + ' 10000')
os.system('python bootstrap_ahab.py boss_iris_2d_north_' + savekey + ' 10000')
os.system('python bootstrap_ahab.py boss_iris_2d_south_' + savekey + ' 10000')

# after this, need to calculate bootstrap envelopes

