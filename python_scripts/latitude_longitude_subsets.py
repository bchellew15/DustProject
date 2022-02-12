# run on ahab
# just a wrapper to help run generate_alphas for different lat / lon ranges

import os

# latitude cuts
os.system('python generate_alphas_ahab.py iris_2d 1 lat_a 10 5 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_b 10 6 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_c 10 7 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_d 10 8 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_e 10 9 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_f 10 10 0 0')

# longitude cuts
os.system('python generate_alphas_ahab.py iris_2d 1 lon_a 10 11 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_b 10 12 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_c 10 13 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_d 10 14 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_e 10 15 0 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_f 10 16 0 0')
