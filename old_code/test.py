from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt

# create figure and axes instances
fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# create polar stereographic Basemap instance.
m = Basemap(projection='stere',lon_0=0,lat_0=30.,lat_ts=45.,\
            width=10000000, height=4000000,
            rsphere=6371200.,resolution='l',area_thresh=10000)

# draw parallels.
parallels = np.arange(0.,90,5.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)

# draw meridians
merid_values = np.arange(0.,360.,10.)
meridians = m.drawmeridians(merid_values,labels=[0,0,0,1],fontsize=10)

plt.show()
