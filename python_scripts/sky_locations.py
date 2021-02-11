#plot sky locations
#color coded according to density and weighted density (see Brandt paper)

import numpy as np
from matplotlib import pyplot as plt
import sys
import matplotlib #for lognorm
from astropy.io import fits
from mpl_toolkits.basemap import Basemap

if len(sys.argv) != 4:
    print("Usage: sky_locations [boss: 0, 1] [save (1), load (0)] [weighted: 0, 1]")
    exit(0)
        
boss = int(sys.argv[1])
save = int(sys.argv[2])
weighted = int(sys.argv[3])

#truncate colormap (https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib)
import matplotlib.colors as colors 

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
cmap = plt.get_cmap('Blues')
new_cmap = truncate_colormap(cmap, 0.2, 1)

if boss: 
    coords = np.loadtxt("/Users/blakechellew/Documents/DustProject/BrandtFiles/BOSS_locations_galactic.txt")

    i100_old = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/SFD_i100_at_BOSS_locations.txt")[:,2]
    i100 = np.load("/Volumes/TOSHIBA/Dust_Overflow/i100_tao_boss_iris.npy", mmap_mode='r')
    i100 = i100[:,2941] # around the halfway point in terms of wavelength

    hdulist = fits.open("/Volumes/TOSHIBA/Dust_Overflow/" + 'skyfibers_lam0.fits')
    plate = hdulist[6].data
    mjd = hdulist[5].data
    plate = 10000*plate + mjd%10000 #plate is now a unique identifier
    
else:
    coords = np.loadtxt('/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/infile.txt')
    fiberinfo = np.loadtxt('/Users/blakechellew/Documents/DustProject/BrandtFiles/fiberinfo_halpha.dat')
    i100_old = np.array(fiberinfo[:, 4])
    i100 = i100_old

    hdulist = fits.open('/Users/blakechellew/Documents/DustProject/BrandtFiles/SDSS_allskyspec.fits')
    plate = np.array(hdulist[0].data)

#mask large i100 values
#all masking is done based on 1D SFD i100
#"mask' is the indices of masked elements
threshold = 10
mask = np.where(i100_old>threshold)[0]
#mask plates with large averages
for p in np.unique(plate):
    if np.mean(i100_old[plate==p]) > threshold:
        mask = np.append(mask, np.where(plate==p)[0])
if boss:
    mask = np.append(mask, 3178) #data at this location is bad
    
longs = coords[:,0]
lats = coords[:,1]

if save and not weighted:
    densities = np.zeros(len(longs))
    for i in range(len(longs)):
        densities[i] = np.sum((longs > longs[i] - 1) * (longs < longs[i] + 1) * (lats > lats[i] - 1) * (lats < lats[i] + 1))

    #apply mask:
    densities = np.delete(densities, mask)
    
    if boss:
        np.save('../alphas_and_stds/boss_skyfiber_densities.npy', densities)
    else:
        np.save('../alphas_and_stds/sdss_skyfiber_densities.npy', densities)

    densities[densities>256] = 256 #max on scale
    
elif not save and not weighted:
    if boss:
        densities = np.load('../alphas_and_stds/boss_skyfiber_densities.npy')
    else:
        densities = np.load('../alphas_and_stds/sdss_skyfiber_densities.npy')
elif save and weighted:
    
    var_p = np.zeros(longs.shape[0]) #variance on a plate
    #print(var_p.shape)
    print(plate.shape)
    print(i100.shape)
    for p in np.unique(plate):
        i100_p = i100[plate == p]
        i100_sqr_p = np.power(i100_p, 2)
        var = np.mean(i100_sqr_p) - np.power(np.mean(i100_p), 2)
        var_p[plate==p] = var
        
    avg_var = np.mean(np.unique(var_p))
    #print("avg var:", avg_var)
    #print("median var:", np.median(np.unique(var_p)))
    #big_vars = var_p[var_p > avg_var]
    #print(len(big_vars))
    #print(len(var_p))

    densities = np.zeros(len(longs))
    for i in range(len(longs)):
        densities[i] = np.sum((longs > longs[i] - 1) * (longs < longs[i] + 1) * (lats > lats[i] - 1) * (lats < lats[i] + 1))
        densities[i] *= var_p[i] / avg_var
        #print("weight:", var_p[i] / avg_var)

    #apply mask:
    densities = np.delete(densities, mask)
    
    if boss:
        np.save('../alphas_and_stds/boss_skyfiber_densities_weighted.npy', densities)
    else:
        np.save('../alphas_and_stds/sdss_skyfiber_densities_weighted.npy', densities)

    densities[densities>256] = 256 #max on scale
    densities[densities<1] = 1 #min on scale
    
elif not save and weighted:
    if boss:
        densities = np.load('../alphas_and_stds/boss_skyfiber_densities_weighted.npy')
    else:
        densities = np.load('../alphas_and_stds/sdss_skyfiber_densities_weighted.npy')

'''
#convert to radians (for regular matplotlib, not basemap)
longs *= (np.pi/180)
lats *= (np.pi/180)

#convert to range -pi to pi:
for i in range(len(longs)):
    if longs[i] > np.pi:
        longs[i] -= 2*np.pi
'''

#mask coordinates:
longs = np.delete(longs, mask)
lats = np.delete(lats, mask)

#move darker points to the top:
longs_lats = [(long, lat) for _,long,lat in sorted(zip(densities,longs,lats))]
longs = [long for long,_ in longs_lats]
lats = [lat for _,lat in longs_lats]
densities = sorted(densities)

#invert coordinates because sky projection:
longs = np.subtract(360, longs)

#use basemap for projection
m = Basemap(projection='moll', lon_0=0) #TEST

#text (gridline labels)
text_xs = np.array([-60, -120, 0, 60, 120, 0, 0, 0])
text_ys = np.array([0, 0, 0, 0, 0, 30, -30, -60])
text_xs += 1
text_ys += 1
x, y = m(text_xs, text_ys)
plt.text(x[0], y[0], '60', fontsize=10)
plt.text(x[1], y[1], '120', fontsize=10)
plt.text(x[2], y[2], '0', fontsize=10)
plt.text(x[3], y[3], '-60', fontsize=10)
plt.text(x[4], y[4], '-120', fontsize=10)
plt.text(x[5], y[5], '30', fontsize=10)
plt.text(x[6], y[6], '-30', fontsize=10)
plt.text(x[7], y[7], '-60', fontsize=10)

#transform data coordinates
x, y = m(longs,lats)

m.drawmapboundary()
m.drawmeridians(np.arange(0,360,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')
m.drawparallels(np.arange(-90,90,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')

m.scatter(x, y, marker='.', s=1, c=densities, cmap=new_cmap, norm=matplotlib.colors.LogNorm())

if weighted:
    cb = m.colorbar(location='bottom', label=r'Sky Fiber Density (deg$^{-2}$)')
else:
    cb = m.colorbar(location='bottom', label=r'Weighted Sky Fiber Density (deg$^{-2}$)')
cb.set_ticks([1, 2, 4, 8, 16, 32, 64, 128, 256])
tick_labels = [1, 2, 4, 8, 16, 32, 64, 128, '256+']
if weighted:
    tick_labels[0] = '<1'
cb.set_ticklabels(tick_labels)

#plt.title("Sky Fiber Density")
plt.show()
