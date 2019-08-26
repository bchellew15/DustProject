import numpy as np
import idlwrap
import pydl
from astropy.coordinates import SkyCoord
from astropy.io import fits
import glob
import os
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
from time import time
import astropy.wcs as wcs #for ad2xy

#next steps:
#test functions independently
#look at output step by step
#better way to convert coordinates? (find examples)
#are there undefined values in the map?

#which file selection method (1 is mine, 0 is original)
file_sel = 1
#which coordinates (0 is SFD, 1 is IRIS example, 2 is BOSS)
coords = 2

start_idx_actual = 0
end_idx_actual = 1000

start_idx_sfd = 0
end_idx_sfd = 91847

#convert coordinates
#result is in degrees
def euler(alpha, delta):
    if coords == 0:
        c = SkyCoord(alpha, delta, frame='galactic', unit="deg")
        alpha = c.fk4.ra.deg
        delta = c.fk4.dec.deg
    elif coords == 1:
        c = SkyCoord(alpha, delta, frame='fk4', unit="deg")
        alpha = c.fk4.ra.deg
        delta = c.fk4.dec.deg
    #fk5 is J2000
    else:
        c = SkyCoord(alpha, delta, frame='fk5', unit="deg")
        alpha = c.fk4.ra.deg
        delta = c.fk4.dec.deg
    return alpha, delta

#works
def nan2undef(l, undef):
    l[np.isnan(l)] = undef

#implemented in IDL
def get_iris(num, direc, band):
    bd = str(band)
    iras_number= str(int(num))
    if num < 100:
        iras_number = '0' + iras_number
    if num < 10:
        iras_number = '0' + iras_number
    hcon = 0
    hnum = str(hcon)
    result = glob.glob(direc+'/?'+iras_number+'?'+bd+'?'+hnum+'.*')
    count = len(result)
    if count > 0:
        fits.setval(result[0], 'CDELT3', value=1) #TEST: modify header
        ifile = fits.open(result[0])
        imap = ifile[0].data[0]
        header = fits.open(result[0])[0].header
        #sxaddpar sets LONPOLE to 180, but already there.
        bad = np.where((imap < -5) + (imap == 0))[0]
        nbbad = len(bad)
        if nbbad > 0:
            map[bad] = -32768
    else:
        print("error: could not find data file for ISSA number", iras_number, "and band", bd)  

    return imap, header

#used formula for tangent projection:
#https://lambda.gsfc.nasa.gov/product/iras/coordproj.cfm
def my_ad2xy(alpha, delta, hi):
    #40 pixels per degree
    origin_x = hi['CRPIX1']
    origin_y = hi['CRPIX2']
    x_deg_pix = hi['CDELT1']
    y_deg_pix = hi['CDELT1']
    x_pix_deg = 1/x_deg_pix
    y_pix_deg = -1/y_deg_pix
    alpha0 = hi['CRVAL1']
    delta0 = hi['CRVAL2']

    #first atempt:
    alpha_res = alpha-alpha0
    alpha_pix = origin_x + x_pix_deg*alpha_res
    delta_res = delta-delta0
    delta_pix = origin_y + y_pix_deg*delta_res

    #better conversion:
    scale = x_pix_deg
    #convert to radians
    alpha = alpha*np.pi/180
    delta = delta*np.pi/180
    alpha0 = alpha0*np.pi/180
    delta0 = delta0*np.pi/180
    A = np.cos(delta) * np.cos(alpha - alpha0)
    F = scale * (180/np.pi)/(np.sin(delta0) * np.sin(delta) + A * np.cos(delta0))
    LINE = -F * (np.cos(delta0) * np.sin(delta) - A * np.sin(delta0))
    SAMPLE = -F * np.cos(delta) * np.sin(alpha - alpha0)
    alpha_pix = origin_x - SAMPLE
    delta_pix = origin_y + LINE
    
    return alpha_pix, delta_pix

#I did my own interpolation using scipy's interp2d
def mbilinear(x, y, array):
    #RectBivariateSpline supposed to be faster? (but only for grid...)
    #f = interpolate.RectBivariateSpline(np.arange(array.shape[0]), np.arange(array.shape[0]), array)
    f = interpolate.interp2d(np.arange(array.shape[0]), np.arange(array.shape[0]), array, bounds_error=False, fill_value=np.nan, kind='linear')
    return np.array([f(i, j) for i, j in zip(x, y)]).reshape(len(x))

def mosaique_iris(alpha, delta, direc):

    band = 4 #i100

    #the IDL code assumes square arrays, I think
    #get_cootype(astr)? I think not needed
    print("convert coordinates from Galactic to Celestial B1950")
    alpha, delta = euler(alpha, delta)
    
    nan2undef(alpha, undef=-32768)
    nan2undef(delta, undef=-32768)

    totsize = len(alpha)
    result = np.zeros(totsize)
    weight = np.zeros(totsize)

    catname = '/Users/blakechellew/Documents/DustProject/IRIS/irispro/info_issa_map4.txt'
    contents = np.loadtxt(catname)
    inum = contents[:,0]
    ramin = contents[:,1]
    ramax = contents[:,2]
    raavg = contents[:,3]
    demin = contents[:,4]
    demax = contents[:,5]
    deavg = contents[:,6]
    medval = contents[:,7]
    noise_key = contents[:,8]

    nb = len(inum)
    id_good = np.zeros(nb)
    print("Check for ISSA maps that intersect with the given header")
    
    ind = np.where((alpha!=-32768) * (delta!=-32768))[0]
    c1min = np.min(alpha[ind])
    c1max = np.max(alpha[ind])
    c2min = np.min(delta[ind])
    c2max = np.max(delta[ind])

    print("ra ranges from", c1min, "to", c1max)
    print("dec ranges from", c2min, "to", c2max)

    #my file selection code
    if file_sel == 1:
        for i in range(nb):
            if (c1min <= ramin[i] and c1max <= ramin[i]) or (c1min >= ramax[i] and c1max >= ramax[i]) \
               or (c2min <= demin[i] and c2max <= demin[i]) or (c2min >= demax[i] and c2max >= demax[i]):
                pass
            else:
                id_good[i] = 1
    #original file selection code
    else:
        for i in range(nb):
            if c1min >= ramin[i] and c1min <= ramax[i] and c2min >= demin[i] and c2min <= demax[i]:
                id_good[i] = 1
            if c1min >= ramin[i] and c1min <= ramax[i] and c2max >= demin[i] and c2max <= demax[i]:
                id_good[i] = 1
            if c1max >= ramin[i] and c1max <= ramax[i] and c2max >= demin[i] and c2max <= demax[i]:
                id_good[i] = 1
            if c1max >= ramin[i] and c1max <= ramax[i] and c2min >= demin[i] and c2min <= demax[i]:
                id_good[i] = 1
            if ramin[i] >= c1min and ramin[i] <= c1max and demin[i] >= c2min and demin[i] <= c2max:
                id_good[i] = 1
            if ramax[i] >= c1min and ramax[i] <= c1max and demin[i] >= c2min and demin[i] <= c2max:
                id_good[i] = 1
            if ramin[i] >= c1min and ramin[i] <= c1max and demax[i] >= c2min and demax[i] <= c2max:
                id_good[i] = 1
            if ramax[i] >= c1min and ramax[i] <= c1max and demax[i] >= c2min and demax[i] <= c2max:
                id_good[i] = 1
            
    ind = np.where(id_good > 0)[0]
    nbind = ind.shape[0]
    if nbind == 0:
        print("No ISSA map corresponds to the header given")
        return -1
    print(nbind, "Issa maps will be combined to produce the mosaic")
    print("PLATE NUMBERS: ", inum[ind])

    #extract the i100 values:
    for i in range(nbind):
        mapi, hi = get_iris(inum[ind[i]], direc=direc, band=band)
        print("file number: ", inum[ind[i]])
        w = wcs.WCS(hi) #, filei . . . 
        #filei.close()
        xi, yi, trash = w.all_world2pix(alpha, delta, np.zeros(len(alpha)), 0)#1) #origin 1 for FITS, zeros b/c axis3 has size 1
        
        '''
        #my own function:
        xi, yi = my_ad2xy(alpha, delta, hi)
        '''
        
        #xi, yi = ad2xy(alpha, delta, hi)
        tempo = mbilinear(xi, yi, mapi)
        tempo[np.isnan(tempo)] = -32768 #TESTING

        
        #uncomment for small number of datapoints
        #print("tempo:")
        #print(tempo)
        
        
        indw = np.where(tempo != -32768)[0]
        nbindw = indw.shape[0]
        if nbindw > 0:
              weight[indw] = weight[indw]+1
              result[indw] = result[indw] + tempo[indw]
            
        #print progress:
        if i % 10 == 0:
            print("progress: ", i)
            
    indw = np.where(weight > 0)[0] #nbindw, complement...
    mask = np.zeros(alpha.shape)
    mask[indw] = True
    indempty = mask==0
    nindempty = tempo[indempty].shape[0]
    if nindempty > 0:
        print("ERROR:", nindempty, "points failed")
    else:
        print("SUCCESS: values found for all points")
    nbindw = indw.shape[0]

    if nbindw > 0:
        result[indw] = result[indw] / weight[indw]
    if nindempty > 0:
        result[indempty] = -32768

    with open('IRIS_values_skyphot.dat', 'w') as f:
        for i in range(totsize):
            f.write(str(result[i]) + "\n")

    return result


#check against actual IRIS values:
if coords == 1:
    #these are FK4 coordinates
    actual_IRIS = np.loadtxt("/Users/blakechellew/Documents/DustProject/IRIS/brandt_iris/IRIS_values_skyphot.dat")
    #these are same coordinates but FK5
    #actual_IRIS2 = np.loadtxt("/Users/blakechellew/Documents/DustProject/IRIS/brandt_iris/skyphot_coords.dat")

    #print("shape: ", actual_IRIS.shape)

    start_idx = start_idx_actual
    end_idx = end_idx_actual

    #a, b either lat/long or not, and B1500 or J2000
    a = actual_IRIS[:,0]
    b = actual_IRIS[:,1]
    iris_i100 = actual_IRIS[:,2]

    #truncate:
    a = a[start_idx:end_idx]
    b = b[start_idx:end_idx]
    iris_i100 = iris_i100[start_idx:end_idx]

    result = mosaique_iris(a, b, '/Users/blakechellew/Documents/DustProject/IRIS/IRISNOHOLES_B4H0')

    from scipy.stats import linregress
    print(linregress(iris_i100, result))

    plt.plot(iris_i100, result, 'k.', markersize=1)
    x = np.arange(np.min(result), np.max(result), 0.1)
    plt.plot(x, x)
    plt.xlabel("iris i100")
    plt.ylabel("my code results")
    plt.show()

#check against SFD values:
elif coords == 0:
    a_b = np.load("/Users/blakechellew/Documents/DustProject/l_b_original.npy")

    start_idx = start_idx_sfd
    end_idx = end_idx_sfd

    a = a_b[start_idx:end_idx,0]
    b = a_b[start_idx:end_idx,1]    

    result = mosaique_iris(a, b, '/Users/blakechellew/Documents/DustProject/IRIS/IRISNOHOLES_B4H0')
    #np.save('iris_i100_at_sfd.npy', result)

    #check correlation
    i100_original = np.load("/Users/blakechellew/Documents/DustProject/i100_1d.npy")[start_idx:end_idx]
    
    
    '''
    #pick certain points:
    #lower branch:
    idx1 = np.where((result<1.9) * (result>1.7)*(i100_original<1.01)*(i100_original>0.95))[0]
    #upper branch:
    idx2 = np.where((result>2.275) * (result<2.295)*(i100_original<1.09)*(i100_original>1.0725))[0]
    print("indices of selected points:")
    print(idx1)
    #124  144 2928 2960
    print(idx2)
    #668  694  726 2132
    
    #print("original:")
    #print(i100_original[:50])
    '''

    from scipy.stats import linregress
    print(linregress(i100_original, result))
    
    plt.plot(i100_original, result, 'k.', markersize=1)
    #consider density plot
    x = np.arange(np.min(result), np.max(result), 0.1)
    plt.plot(x, x)
    plt.xlabel("SFD i100")
    plt.ylabel("IRIS i100")
    plt.show()

#get IRIS values at BOSS locations:
else:
    #load coordinates:
    a_b = np.loadtxt("/Users/blakechellew/Documents/DustProject/BrandtFiles/BOSS_locations.txt")

    a = a_b[:,0]
    b = a_b[:,1] 

    #get i100 values:
    result = mosaique_iris(a, b, '/Users/blakechellew/Documents/DustProject/IRIS/IRISNOHOLES_B4H0')

    #save i100 values:
    np.save('iris_i100_at_boss.npy', result)

    #check correlation with SFD BOSS values...





'''
#test mbilinear

mapi, hi = get_iris(100, '/Users/blakechellew/Documents/DustProject/IRIS/IRISNOHOLES_B4H0', band=4)
test = mbilinear([1000], [2], mapi)
print(test)
exit(0)
'''

