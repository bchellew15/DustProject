#a script to loop over other scripts with various parameter combinations

import os

thresholds = ['5', '7.5', '10', '12.5', '15']
threshold_labels = ['5', '75', '10', '125', '15']

#first I'll do SFD SDSS at 5 thresholds
#then SFD BOSS at the same 5 thresholds
#what about IRIS data? (maybe iris 2d at some point)

for i in range(5):
    os.system('python reproduce_figs_boss.py ' + 'iris_1d' + ' ' + '1' + ' ' + 'boss_iris_1d_91119_' + threshold_labels[i] + ' ' + thresholds[i])
