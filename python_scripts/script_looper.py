#a script to loop over other scripts with various parameter combinations

import os

#thresholds = ['5', '7.5', '10', '12.5', '15']
#threshold_labels = ['5', '75', '10', '125', '15']
thresholds = ['20', '25', '30']
threshold_labels = ['20', '25', '30']

'''
#first I'll do SFD SDSS at 5 thresholds
#then SFD BOSS at the same 5 thresholds
#what about IRIS data? (maybe iris 2d at some point)

for i in range(len(thresholds)):
    os.system('python reproduce_figs_boss.py ' + 'iris_1d' + ' ' + '1' + ' ' + 'boss_iris_1d_92719_' + threshold_labels[i] + ' ' + thresholds[i])
'''

'''
#bootstrapping (generate xx and yx):
commands = ['python reproduce_figs_ahab.py 2d 0 sdss_2d_102019 10 1',
            'python reproduce_figs_ahab.py iris_1d 0 sdss_iris_1d_102019 10 1',
            'python reproduce_figs_ahab.py iris 0 sdss_iris_2d_102019 10 1',
            'python reproduce_figs_ahab.py 1d 1 boss_1d_102019 10 1',
            'python reproduce_figs_ahab.py 2d 1 boss_2d_102019 10 1',
            'python reproduce_figs_ahab.py iris_1d 1 boss_iris_1d_102019 10 1',
            'python reproduce_figs_ahab.py iris 1 boss_iris_2d_102019 10 1']

#for c in commands:
#    os.system(c)
'''



#bootstrapping (generate bootstrap samples from xx and yx)
commands = ['python bootstrap_ahab.py sdss_2d_102019 10000',
            'python bootstrap_ahab.py sdss_iris_1d_102019 10000',
            'python bootstrap_ahab.py sdss_iris_2d_102019 10000',
            'python bootstrap_ahab.py boss_1d_102019 10000',
            'python bootstrap_ahab.py boss_2d_102019 10000',
            'python bootstrap_ahab.py boss_iris_1d_102019 10000',
            'python bootstrap_ahab.py boss_iris_2d_102019 10000']

#for c in commands:
#    os.system(c)

os.system(commands[4])
print("command 4 done")
os.system(commands[5])
print("command 4 done")
os.system(commands[6])
print("command 4 done")

