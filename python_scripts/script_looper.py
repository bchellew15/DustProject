#loop over other scripts with various parameter combinations

#goal: input a date, and it generates all the relevant alpha files and bootstrap files
#set it up so can do each step individually if needed

import os

'''
#thresholds = ['5', '7.5', '10', '12.5', '15']
#threshold_labels = ['5', '75', '10', '125', '15']
thresholds = ['17.5', '22.5', '27.5', '30', '32.5']
threshold_labels = ['175', '225', '275', '30', '325']

#first I'll do SFD SDSS at 5 thresholds
#then SFD BOSS at the same 5 thresholds
#what about IRIS data? (maybe iris 2d at some point)

for i in range(len(thresholds)):
    os.system('python reproduce_figs_ahab.py ' + 'iris_1d' + ' ' + '1' + ' ' + 'boss_iris_1d_111119_' + threshold_labels[i] + ' ' + thresholds[i])
'''

'''
#bootstrapping (generate xx and yx):
commands = ['python reproduce_figs_ahab.py 2d 0 sdss_2d_102019 10 0 1',
            'python reproduce_figs_ahab.py iris_1d 0 sdss_iris_1d_102019 10 0 1',
            'python reproduce_figs_ahab.py iris 0 sdss_iris_2d_102019 10 0 1',
            'python reproduce_figs_ahab.py 1d 1 boss_1d_102019 10 0 1',
            'python reproduce_figs_ahab.py 2d 1 boss_2d_102019 10 0 1',
            'python reproduce_figs_ahab.py iris_1d 1 boss_iris_1d_102019 10 0 1',
            'python reproduce_figs_ahab.py iris 1 boss_iris_2d_102019 10 0 1']

#for north and south BOSS:
commands = ['python reproduce_figs_ahab.py iris 1 north_111719 10 1 1',
            'python reproduce_figs_ahab.py iris 1 south_111719 10 2 1']

for c in commands:
    os.system(c)
'''


#bootstrapping (generate bootstrap samples from xx and yx)
commands = ['python bootstrap_ahab.py sdss_1d_101819 9000',
            'python bootstrap_ahab.py sdss_2d_102019 9000',
            'python bootstrap_ahab.py sdss_iris_1d_102019 9000',
            'python bootstrap_ahab.py sdss_iris_2d_102019 9000',
            'python bootstrap_ahab.py boss_1d_102019 9000',
            'python bootstrap_ahab.py boss_2d_102019 9000',
            'python bootstrap_ahab.py boss_iris_1d_102019 9000',
            'python bootstrap_ahab.py boss_iris_2d_102019 9000']

commands = ['python bootstrap_ahab.py north_111719 500',
            'python bootstrap_ahab.py south_111719 500']

commands = ['python bootstrap_ahab.py south_111719 500']


for c in commands:
    os.system(c)

