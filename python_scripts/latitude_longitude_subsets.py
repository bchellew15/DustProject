# run on ahab
# just a wrapper to help run generate_alphas for different lat / lon ranges

import os
from multiprocessing import Pool
from generate_plots import generate_binned_alphas

num_processes = 6

####################################

# # latitude cuts
# os.system('python generate_alphas_ahab.py iris_2d 1 lat_a 10 5 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lat_b 10 6 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lat_c 10 7 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lat_d 10 8 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lat_e 10 9 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lat_f 10 10 0 0')
#
# # longitude cuts
# os.system('python generate_alphas_ahab.py iris_2d 1 lon_a 10 11 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lon_b 10 12 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lon_c 10 13 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lon_d 10 14 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lon_e 10 15 0 0')
# os.system('python generate_alphas_ahab.py iris_2d 1 lon_f 10 16 0 0')

# then bootstrap
os.system('python generate_alphas_ahab.py iris_2d 1 lat_a 10 5 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_b 10 6 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_c 10 7 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_d 10 8 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_e 10 9 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lat_f 10 10 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_a 10 11 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_b 10 12 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_c 10 13 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_d 10 14 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_e 10 15 1 0')
os.system('python generate_alphas_ahab.py iris_2d 1 lon_f 10 16 1 0')

# then generate bootstrap samples
def run_bootstrapping(savekey):
    execute_string = 'python bootstrap_ahab.py ' + savekey + ' 2000'
    os.system(execute_string)

inputs = ['lat_a', 'lat_b', 'lat_c', 'lat_d', 'lat_e', 'lat_f', 'lon_a', 'lon_b', 'lon_c', 'lon_d', 'lon_e', 'lon_f']
with Pool(num_processes) as p:
    outputs = p.imap_unordered(run_bootstrapping, inputs)
    for output in outputs:
        print(output)

# generate envelopes
wavelength_boss = np.load('/home/blakechellew/alphas_and_stds/wavelength_boss.npy')
# load bootstrap alphas
bootstrap_alphas_lat_a = np.load('../data/bootstrap_alphas_lat_a.npy')
bootstrap_alphas_lat_b = np.load('../data/bootstrap_alphas_lat_b.npy')
bootstrap_alphas_lat_c = np.load('../data/bootstrap_alphas_lat_c.npy')
bootstrap_alphas_lat_d = np.load('../data/bootstrap_alphas_lat_d.npy')
bootstrap_alphas_lat_e = np.load('../data/bootstrap_alphas_lat_e.npy')
bootstrap_alphas_lat_f = np.load('../data/bootstrap_alphas_lat_f.npy')
bootstrap_alphas_lon_a = np.load('../data/bootstrap_alphas_lon_a.npy')
bootstrap_alphas_lon_b = np.load('../data/bootstrap_alphas_lon_b.npy')
bootstrap_alphas_lon_c = np.load('../data/bootstrap_alphas_lon_c.npy')
bootstrap_alphas_lon_d = np.load('../data/bootstrap_alphas_lon_d.npy')
bootstrap_alphas_lon_e = np.load('../data/bootstrap_alphas_lon_e.npy')
bootstrap_alphas_lon_f = np.load('../data/bootstrap_alphas_lon_f.npy')
# load bootstrap stds
bootstrap_alpha_stds_lat_a = np.load('../data/bootstrap_alpha_stds_lat_a.npy')
bootstrap_alpha_stds_lat_b = np.load('../data/bootstrap_alpha_stds_lat_b.npy')
bootstrap_alpha_stds_lat_c = np.load('../data/bootstrap_alpha_stds_lat_c.npy')
bootstrap_alpha_stds_lat_d = np.load('../data/bootstrap_alpha_stds_lat_d.npy')
bootstrap_alpha_stds_lat_e = np.load('../data/bootstrap_alpha_stds_lat_e.npy')
bootstrap_alpha_stds_lat_f = np.load('../data/bootstrap_alpha_stds_lat_f.npy')
bootstrap_alpha_stds_lon_a = np.load('../data/bootstrap_alpha_stds_lon_a.npy')
bootstrap_alpha_stds_lon_b = np.load('../data/bootstrap_alpha_stds_lon_b.npy')
bootstrap_alpha_stds_lon_c = np.load('../data/bootstrap_alpha_stds_lon_c.npy')
bootstrap_alpha_stds_lon_d = np.load('../data/bootstrap_alpha_stds_lon_d.npy')
bootstrap_alpha_stds_lon_e = np.load('../data/bootstrap_alpha_stds_lon_e.npy')
bootstrap_alpha_stds_lon_f = np.load('../data/bootstrap_alpha_stds_lon_f.npy')
# generate binned stds:
_, bootstrap_binned_alphas_lat_a, _ = generate_binned_alphas(bootstrap_alphas_lat_a, bootstrap_alpha_stds_lat_a, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lat_b, _ = generate_binned_alphas(bootstrap_alphas_lat_b, bootstrap_alpha_stds_lat_b, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lat_c, _ = generate_binned_alphas(bootstrap_alphas_lat_c, bootstrap_alpha_stds_lat_c, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lat_d, _ = generate_binned_alphas(bootstrap_alphas_lat_d, bootstrap_alpha_stds_lat_d, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lat_e, _ = generate_binned_alphas(bootstrap_alphas_lat_e, bootstrap_alpha_stds_lat_e, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lat_f, _ = generate_binned_alphas(bootstrap_alphas_lat_f, bootstrap_alpha_stds_lat_f, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lon_a, _ = generate_binned_alphas(bootstrap_alphas_lon_a, bootstrap_alpha_stds_lon_a, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lon_b, _ = generate_binned_alphas(bootstrap_alphas_lon_b, bootstrap_alpha_stds_lon_b, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lon_c, _ = generate_binned_alphas(bootstrap_alphas_lon_c, bootstrap_alpha_stds_lon_c, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lon_d, _ = generate_binned_alphas(bootstrap_alphas_lon_d, bootstrap_alpha_stds_lon_d, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lon_e, _ = generate_binned_alphas(bootstrap_alphas_lon_e, bootstrap_alpha_stds_lon_e, wavelength_boss, boss=True)
_, bootstrap_binned_alphas_lon_f, _ = generate_binned_alphas(bootstrap_alphas_lon_f, bootstrap_alpha_stds_lon_f, wavelength_boss, boss=True)
# percentiles: lower
bootstrap_binned_lower_lat_a = np.percentile(bootstrap_binned_alphas_lat_a, 16, axis=0)
bootstrap_binned_lower_lat_b = np.percentile(bootstrap_binned_alphas_lat_b, 16, axis=0)
bootstrap_binned_lower_lat_c = np.percentile(bootstrap_binned_alphas_lat_c, 16, axis=0)
bootstrap_binned_lower_lat_d = np.percentile(bootstrap_binned_alphas_lat_d, 16, axis=0)
bootstrap_binned_lower_lat_e = np.percentile(bootstrap_binned_alphas_lat_e, 16, axis=0)
bootstrap_binned_lower_lat_f = np.percentile(bootstrap_binned_alphas_lat_f, 16, axis=0)
bootstrap_binned_lower_lon_a = np.percentile(bootstrap_binned_alphas_lon_a, 16, axis=0)
bootstrap_binned_lower_lon_b = np.percentile(bootstrap_binned_alphas_lon_b, 16, axis=0)
bootstrap_binned_lower_lon_c = np.percentile(bootstrap_binned_alphas_lon_c, 16, axis=0)
bootstrap_binned_lower_lon_d = np.percentile(bootstrap_binned_alphas_lon_d, 16, axis=0)
bootstrap_binned_lower_lon_e = np.percentile(bootstrap_binned_alphas_lon_e, 16, axis=0)
bootstrap_binned_lower_lon_f = np.percentile(bootstrap_binned_alphas_lon_f, 16, axis=0)
# percentiles: upper
bootstrap_binned_upper_lat_a = np.percentile(bootstrap_binned_alphas_lat_a, 84, axis=0)
bootstrap_binned_upper_lat_b = np.percentile(bootstrap_binned_alphas_lat_b, 84, axis=0)
bootstrap_binned_upper_lat_c = np.percentile(bootstrap_binned_alphas_lat_c, 84, axis=0)
bootstrap_binned_upper_lat_d = np.percentile(bootstrap_binned_alphas_lat_d, 84, axis=0)
bootstrap_binned_upper_lat_e = np.percentile(bootstrap_binned_alphas_lat_e, 84, axis=0)
bootstrap_binned_upper_lat_f = np.percentile(bootstrap_binned_alphas_lat_f, 84, axis=0)
bootstrap_binned_upper_lon_a = np.percentile(bootstrap_binned_alphas_lon_a, 84, axis=0)
bootstrap_binned_upper_lon_b = np.percentile(bootstrap_binned_alphas_lon_b, 84, axis=0)
bootstrap_binned_upper_lon_c = np.percentile(bootstrap_binned_alphas_lon_c, 84, axis=0)
bootstrap_binned_upper_lon_d = np.percentile(bootstrap_binned_alphas_lon_d, 84, axis=0)
bootstrap_binned_upper_lon_e = np.percentile(bootstrap_binned_alphas_lon_e, 84, axis=0)
bootstrap_binned_upper_lon_f = np.percentile(bootstrap_binned_alphas_lon_f, 84, axis=0)
# save lower
np.save('../data/bootstrap_binned_lower_lat_a.npy', bootstrap_binned_lower_lat_a)
np.save('../data/bootstrap_binned_lower_lat_b.npy', bootstrap_binned_lower_lat_b)
np.save('../data/bootstrap_binned_lower_lat_c.npy', bootstrap_binned_lower_lat_c)
np.save('../data/bootstrap_binned_lower_lat_d.npy', bootstrap_binned_lower_lat_d)
np.save('../data/bootstrap_binned_lower_lat_e.npy', bootstrap_binned_lower_lat_e)
np.save('../data/bootstrap_binned_lower_lat_f.npy', bootstrap_binned_lower_lat_f)
np.save('../data/bootstrap_binned_lower_lon_a.npy', bootstrap_binned_lower_lon_a)
np.save('../data/bootstrap_binned_lower_lon_b.npy', bootstrap_binned_lower_lon_b)
np.save('../data/bootstrap_binned_lower_lon_c.npy', bootstrap_binned_lower_lon_c)
np.save('../data/bootstrap_binned_lower_lon_d.npy', bootstrap_binned_lower_lon_d)
np.save('../data/bootstrap_binned_lower_lon_e.npy', bootstrap_binned_lower_lon_e)
np.save('../data/bootstrap_binned_lower_lon_f.npy', bootstrap_binned_lower_lon_f)
# save upper
np.save('../data/bootstrap_binned_upper_lat_a.npy', bootstrap_binned_upper_lat_a)
np.save('../data/bootstrap_binned_upper_lat_b.npy', bootstrap_binned_upper_lat_b)
np.save('../data/bootstrap_binned_upper_lat_c.npy', bootstrap_binned_upper_lat_c)
np.save('../data/bootstrap_binned_upper_lat_d.npy', bootstrap_binned_upper_lat_d)
np.save('../data/bootstrap_binned_upper_lat_e.npy', bootstrap_binned_upper_lat_e)
np.save('../data/bootstrap_binned_upper_lat_f.npy', bootstrap_binned_upper_lat_f)
np.save('../data/bootstrap_binned_upper_lon_a.npy', bootstrap_binned_upper_lon_a)
np.save('../data/bootstrap_binned_upper_lon_b.npy', bootstrap_binned_upper_lon_b)
np.save('../data/bootstrap_binned_upper_lon_c.npy', bootstrap_binned_upper_lon_c)
np.save('../data/bootstrap_binned_upper_lon_d.npy', bootstrap_binned_upper_lon_d)
np.save('../data/bootstrap_binned_upper_lon_e.npy', bootstrap_binned_upper_lon_e)
np.save('../data/bootstrap_binned_upper_lon_f.npy', bootstrap_binned_upper_lon_f)
