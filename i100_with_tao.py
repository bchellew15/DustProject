from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

#multiply i100 by correction factor that takes optical depth into account
#note: it does not use updated values of i100

#load data
fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
i100 = np.array(fiberinfo[:, 4])  # 100 micron intensity (MJy/sr)

#convert i100 using optical depth:
taos = np.load("/Users/blakechellew/Documents/DustProject/taos.npy")
i100_factor = np.divide(1 - np.exp(np.negative(taos)), taos)
i100 = i100.reshape(len(i100), 1)
i100 = i100 * i100_factor

np.save("/Users/blakechellew/Documents/DustProject/i100_tao.npy", i100)
