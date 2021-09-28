# wrapper to run rad_transfer_v2.py a few times.

# 34 is t5e9, z=.02  (33 is z=.008)
# 37 is t9e9, 1 is const

import os

# os.system('python radiative_transfer_v2.py 34 1 t5e9_grid_150_wd.npy')
# os.system('python radiative_transfer_v2.py 37 1 t9e9_grid_150_wd.npy')
# os.system('python radiative_transfer_v2.py 1 1 cst_grid_150_wd.npy')
os.system('python radiative_transfer_v2.py 34 0 t5e9_grid_150_zd.npy')
os.system('python radiative_transfer_v2.py 37 0 t9e9_grid_150_zd.npy')
os.system('python radiative_transfer_v2.py 1 0 cst_grid_150_zd.npy')