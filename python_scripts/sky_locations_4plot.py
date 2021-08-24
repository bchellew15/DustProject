#MODIFIED version of SKY_LOCATIONS.PY
#plot sky locations
#color coded according to density and weighted density (see Brandt paper)

import numpy as np
from matplotlib import pyplot as plt
import sys
import matplotlib #for lognorm
from astropy.io import fits
from mpl_toolkits.basemap import Basemap

'''
COMMAND LINE ARGS:
if save = 0, save the figure, don't display
otherwise, display but don't save
'''

#command line args:
if len(sys.argv) != 2:
    print("Usage: sky_locations_4plot.py [save: 0, 1]")
    exit(0)
save = int(sys.argv[1])

#truncate colormap (https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib)
import matplotlib.colors as colors 

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

trunc_map = True
if not trunc_map:
    #new_cmap = plt.get_cmap('YlGnBu')
    new_cmap = plt.get_cmap('magma_r')
    
else:
    #cmap = plt.get_cmap('cubehelix_r')
    #cmap = plt.get_cmap('BuPu')
    cmap = plt.get_cmap('PuBuGn')
    new_cmap = truncate_colormap(cmap, 0.1, 1)

#boss:
coords = np.loadtxt("/Users/blakechellew/Documents/DustProject/BrandtFiles/BOSS_locations_galactic.txt")
i100_old = np.loadtxt("/Users/blakechellew/Documents/DustProject/SFD_Maps/CodeC/SFD_i100_at_BOSS_locations.txt")[:,2]
i100 = np.load("/Volumes/TOSHIBA/Dust_Overflow/i100_tao_boss_iris.npy", mmap_mode='r')
i100 = i100[:,235] #around the halfway point
hdulist = fits.open("/Volumes/TOSHIBA/Dust_Overflow/" + 'skyfibers_lam0.fits')
plate = hdulist[6].data
mjd = hdulist[5].data
plate = 10000*plate + mjd%10000 #plate is now a unique identifier

#mask large i100 values
#all masking is done based on 1D SFD i100
#"mask' is the indices of masked elements
threshold = 10
mask = np.where(i100_old>threshold)[0]
#mask plates with large averages
for p in np.unique(plate):
    if np.mean(i100_old[plate==p]) > threshold:
        mask = np.append(mask, np.where(plate==p)[0])
mask = np.append(mask, 3178) #data at this location is bad


fig = plt.figure(figsize=(12, 8), dpi=200)
fig.subplots_adjust(hspace=0.1) # height spaces
fig.subplots_adjust(wspace=0.15) # width spaces
plt.tight_layout()

ax = fig.add_subplot(222)

#boss unweighted
longs = np.copy(coords[:,0])
lats = np.copy(coords[:,1])
densities = np.load('../alphas_and_stds/boss_skyfiber_densities.npy')
boss_colorbar_max = 256  # this is 256 * (9/4) to compare with SDSS
densities[densities > boss_colorbar_max] = boss_colorbar_max

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
plt.text(0, 1, 'BOSS', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=18)

#transform data coordinates
x, y = m(longs,lats)

m.drawmapboundary(linewidth=1, zorder=-1)
m.drawmeridians(np.arange(0,360,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')
m.drawparallels(np.arange(-90,90,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')

m.scatter(x, y, marker='.', s=1, c=densities, cmap=new_cmap, norm=matplotlib.colors.LogNorm())

cb = m.colorbar(location='bottom', label=r'Sky Fiber Density (deg$^{-2}$)')
cb.set_ticks([1, 2, 4, 8, 16, 32, 64, 128, 256])
tick_labels = [1, 2, 4, 8, 16, 32, 64, 128, '256+']
cb.set_ticklabels(tick_labels)


# subtract from 360 to invert the coordinates
x1, y1 = m(360-20.7, 70.8)
x2, y2 = m(360-279.3, 75.3)
x3, y3 = m(360-88.4, -26.1)
m.scatter(x1, y1, marker='.', s=10, c='r')
m.scatter(x2, y2, marker='.', s=10, c='r')
m.scatter(x3, y3, marker='.', s=10, c='r')

#plt.title("Sky Fiber Density")

ax = fig.add_subplot(224)

#boss weighted
longs = np.copy(coords[:,0])
lats = np.copy(coords[:,1])
densities = np.load('../alphas_and_stds/boss_skyfiber_densities_weighted.npy')
densities[densities > boss_colorbar_max] = boss_colorbar_max #max on scale
densities[densities<1] = 1 #min on scale

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
plt.text(0, 1, 'BOSS', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=18)

#transform data coordinates
x, y = m(longs,lats)

m.drawmapboundary(linewidth=1, zorder=-1)
m.drawmeridians(np.arange(0,360,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')
m.drawparallels(np.arange(-90,90,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')

m.scatter(x, y, marker='.', s=1, c=densities, cmap=new_cmap, norm=matplotlib.colors.LogNorm())

cb = m.colorbar(location='bottom', label=r'Weighted Sky Fiber Density (deg$^{-2}$)')
cb.set_ticks([1, 2, 4, 8, 16, 32, 64, 128, 256])
tick_labels = [1, 2, 4, 8, 16, 32, 64, 128, '256+']
tick_labels[0] = '<1'
cb.set_ticklabels(tick_labels)

m.scatter(x1, y1, marker='.', s=10, c='r')
m.scatter(x2, y2, marker='.', s=10, c='r')
m.scatter(x3, y3, marker='.', s=10, c='r')

#plt.title("Sky Fiber Density")


#sdss
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

#plt.subplot(2, 2, 1)
ax = fig.add_subplot(221)
        
#sdss unweighted
longs = np.copy(coords[:,0])
lats = np.copy(coords[:,1])
densities = np.load('../alphas_and_stds/sdss_skyfiber_densities.npy')
densities[densities>256] = 256 #max on scale

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
plt.text(0, 1, 'SDSS', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=18)

#transform data coordinates
x, y = m(longs,lats)

m.drawmapboundary(linewidth=1, zorder=-1)
m.drawmeridians(np.arange(0,360,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')
m.drawparallels(np.arange(-90,90,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')

m.scatter(x, y, marker='.', s=1, c=densities, cmap=new_cmap, norm=matplotlib.colors.LogNorm())

cb = m.colorbar(location='bottom', label=r'Sky Fiber Density (deg$^{-2}$)')
cb.set_ticks([1, 2, 4, 8, 16, 32, 64, 128, 256])
tick_labels = [1, 2, 4, 8, 16, 32, 64, 128, '256+']
cb.set_ticklabels(tick_labels)

#plt.title("Sky Fiber Density")

#plt.subplot(2, 2, 3)
ax = fig.add_subplot(223)

#sdss weighted
longs = np.copy(coords[:,0])
lats = np.copy(coords[:,1])
densities = np.load('../alphas_and_stds/sdss_skyfiber_densities_weighted.npy')
densities[densities>256] = 256 #max on scale
densities[densities<1] = 1 #min on scale

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
plt.text(0, 1, 'SDSS', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=18)

#transform data coordinates
x, y = m(longs,lats)

m.drawmapboundary(linewidth=1, zorder=-1)
m.drawmeridians(np.arange(0,360,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')
m.drawparallels(np.arange(-90,90,30), zorder=0, linewidth=0.5, dashes=[1, 0], color='gray')

m.scatter(x, y, marker='.', s=1, c=densities, cmap=new_cmap, norm=matplotlib.colors.LogNorm())

cb = m.colorbar(location='bottom', label=r'Weighted Sky Fiber Density (deg$^{-2}$)')
cb.set_ticks([1, 2, 4, 8, 16, 32, 64, 128, 256])
tick_labels = [1, 2, 4, 8, 16, 32, 64, 128, '256+']
tick_labels[0] = '<1'
cb.set_ticklabels(tick_labels)

#plt.title("Sky Fiber Density")

if save:
    plt.savefig('../paper_figures/skyfibers_4plot_082221.png', bbox_inches='tight')
else:
    plt.show()

