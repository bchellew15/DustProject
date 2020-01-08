#look at stacked spectra
#(an array of spectra where each row is a different spectrum)

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


#jacknife: plot the first few
a = np.load("../alphas_and_stds/jacknife_alphas_south_111719.npy")

print(a.shape)

wavelength = np.load('../alphas_and_stds/wavelength_boss.npy')

plt.plot(wavelength, a[0], drawstyle='steps')
plt.xlim(5840, 5855)
plt.ylim(0, 1)
plt.show()
#plt.savefig('test_fig.png')
