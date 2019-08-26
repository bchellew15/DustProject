#jiaoyue code: no tao

'''
my comments:

masking: main issue; 0 when taking sums, but does not mask when taking average of x. I deleted the masked stuff.
doesn't affect result, but -1 indexing is unnecessary for y = lambda*ilambda
ivar needs to be unit-converted AND divided by wavelength^2

'''

from astropy.io import fits
import numpy as np

fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
pn = fiberinfo[:, 0]
l = fiberinfo[:, 2]
b = fiberinfo[:, 3]
i100 = fiberinfo[:, 4]

hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
plate = hdulist[0].data
wavelength = hdulist[1].data
flambda = hdulist[2].data
ivar = hdulist[3].data

ilambda=np.array(flambda)

#mask the high-intensity values: (blake method)
mask_indices = np.arange(len(i100))[i100>10]
l = np.delete(l, mask_indices)
b = np.delete(b, mask_indices)
i100 = np.delete(i100, mask_indices)
plate = np.delete(plate, mask_indices)
pn = np.delete(pn, mask_indices)
ilambda = np.delete(ilambda, mask_indices, 0)
ivar = np.delete(ivar, mask_indices, 0)

'''
#intensity = flux per sr
for i in range(91847):
    if i100[i]>=10:
        ivar[i]=0
#exclude fibers >100MJy
'''

#multiply by nu:
i100=[i*3*10**12 for i in i100]

from collections import defaultdict
import collections
#plateN,i100
dict=defaultdict(list)
for i, j in zip(pn, i100):
    dict[i].append(j)
mean_dict={} #avg i100 over plate
for k,v in dict.items():
    mean_dict[k]=[sum(v)/float(len(v))]*len(v)#to have same len as i100 in each plate
#i100,avg_i100,inverse variances lists after taking out bad fibers

odict=collections.OrderedDict(sorted(dict.items()))
omean_dict=collections.OrderedDict(sorted(mean_dict.items()))
#keep plate orders
import itertools
mean=list(itertools.chain(*omean_dict.values()))

excess=[a-b for a,b in zip(i100,mean)]
x=[i for i in excess]#times 100micron in frequency

#blake add:
ivar /= np.array([pow(w, 2) for w in wavelength])
#ivar = pow(1.619*10**(-10), 2) * ivar #no need for this if converting x only

#y=lambda*ilambda
wavelengthMatrix=np.zeros((len(i100),len(wavelength)))
y=np.zeros((len(i100),len(wavelength)))
for i in range(len(i100)):
    wavelengthMatrix[i-1]=wavelength
    y[i-1]=wavelengthMatrix[i-1]*ilambda[i-1]
    
#make x and ivar array
xa=np.array(x)
xa *= 1.619*10**(-10) #unit conversion

ivara=np.array(ivar)
#y*ivar
yi=np.zeros((len(i100),len(wavelength)))
for i in range(len(i100)):
    yi[i-1]=y[i-1]*ivar[i-1]

#print("ivar: ", ivara)

'''
#save x, y, ivar
np.save('/Users/blakechellew/Documents/DustProject/jiaoyue_x.npy', xa)
np.save('/Users/blakechellew/Documents/DustProject/jiaoyue_y.npy', y)
np.save('/Users/blakechellew/Documents/DustProject/jiaoyue_ivar.npy', ivara)
'''
    
#sum over fibers for each wavelength    
top=np.dot(yi.T,xa)
bottom=np.dot(ivara.T,xa**2)
alpha=top/bottom

print("alphas before bin: ")
for i in range(100):
    print(alpha[i])

np.save('/Users/blakechellew/Documents/DustProject/jiaoyue_alphas.npy', alpha)
    


#blake add:
stds = np.sqrt(1/bottom)
print("stds: ", stds)

from scipy import stats
from matplotlib import pyplot as plt
import fileinput

bin_means, bin_edges, binnumber = stats.binned_statistic(wavelength, alpha,statistic='mean', bins=108)
bin_width = (bin_edges[1] - bin_edges[0])
plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g')
plt.xlabel('wavelength(A)')
plt.ylabel('Alpha')
plt.axis([3900, 9200,0,0.30]) #was 0.28, not 0.30
#plt.savefig('corrSpec.png', dpi=2000)

plt.show()
