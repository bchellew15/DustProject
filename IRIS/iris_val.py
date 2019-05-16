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

'''
issues:
subtracting 1 in the ad2xy method


'''


#should work
#if select == 2, galactic to RA / DEC
#fk4 == 1 means convert to B1950 (0 means J2000)
#result is in degrees
def euler(alpha, delta, select, fk4):
    if select == 2:
        #convert galactic to celestial
        c = SkyCoord(alpha, delta, frame='galactic', unit="deg")
        alpha = c.fk4.ra.deg
        delta = c.fk4.dec.deg
    elif select == 3:
        #convert ecliptic to celestial
        print("not implemented yet")
    return alpha, delta

#not needed
#implemented in IDL
def mprecess():
    return True

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

#can't find
#I think not needed
def extast(hi):
    return True

#can't find
#source: acrider.wordpress.com/2015/02/04/revised-python-code-for-converting-ra-dec-to-x-y/
#added stuff for rotation from:
#http://danmoser.github.io/notes/gai_fits-imgs.html
#it mentioned a slight modification to the original code, a -1
#assuming no rotation so CROTA2 = 0
def ad2xy(alpha, delta, header):

    galaxy_RA = alpha
    galaxy_DEC = delta
    FITS_header = header

    if (FITS_header['CTYPE1'] != 'RA---TAN') or (FITS_header['CTYPE2'] != 'DEC--TAN') :
        print('ERROR: Wrong CTYPE1 or CTYPE2 in galaxy_X!')
        print('CTYPE1 = ', FITS_header['CTYPE1'])
        print( 'CTYPE2 = ', FITS_header['CTYPE2'])

    # These values allow translation from RA,DEC to X,Y and vice versa.
    crpix1 = FITS_header['CRPIX1'] - 1 # X of reference pixel
    crpix2 = FITS_header['CRPIX2'] - 1 # Y of reference pixel
    crval1 = FITS_header['CRVAL1'] # RA of reference pixel
    crval2 = FITS_header['CRVAL2'] # DEC of reference pixel

    #BC: I changed this
    delt1 = FITS_header['CDELT1']
    delt2 = FITS_header['CDELT2'] 
    cd1_1 = delt1*np.cos(0)
    cd1_2 = delt2*np.sin(0)
    cd2_1 = delt1*np.sin(0)
    cd2_2 = delt2*np.cos(0)

    # Find the X,Y values of the galaxy's RA and DEC.
    det = cd1_1 * cd2_2 - cd1_2 * cd2_1

    cdinv11 = cd2_2 / det
    cdinv12 = -cd1_2 / det
    cdinv21 = -cd2_1 / det
    cdinv22 = cd1_1 / det

    ra0 = crval1 / 180.0 * np.pi
    dec0 = crval2 / 180.0 * np.pi
    ra = galaxy_RA / 180.0 * np.pi
    dec = galaxy_DEC / 180.0 * np.pi

    bottom = np.sin(dec)*np.sin(dec0) + np.cos(dec)*np.cos(dec0)*np.cos(ra-ra0)

    xi = np.cos(dec) * np.sin(ra-ra0) / bottom
    eta = (np.sin(dec)*np.cos(dec0) - np.cos(dec)*np.sin(dec0)*np.cos(ra-ra0)) / bottom
    xi = xi * 180.0 / np.pi
    eta = eta * 180.0 / np.pi

    galaxy_X = cdinv11 * xi + cdinv12 * eta + crpix1
    galaxy_Y = cdinv21 * xi + cdinv22 * eta + crpix2

    return galaxy_X, galaxy_Y

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

#implemented in IDL
#I just did my own interpolation, might go back and update
def mbilinear(x, y, array):
    '''
    #I assume 1d x and y
    six = len(x)
    siy = len(y)
    siA = array.shape
    missing = -32768
    Nax = siA[0]
    Nay = siA[1]
    Nx = six
    Ny = siy
    output = np.zeros(Nx, Ny) #&output(*,*)=missing?

    min_array = np.argmin(array)
    ind = np.where(array != missing)[0]
    count = ind.shape[0]
    if count != 0:
        min_array = np.argmin(array[ind])
    #I don't understant the indbad line
    indbad = np.where((x < 0) + (x > Nax-1) or (y < 0) or (y > Nay-1))
    countbad = indbad.shape[0]
    inter_percent = 1*(Nx*Ny-countbad)/Nx/Ny*100

    for j in range(Ny):
        ind = np.where(x)[]
    '''
    
    #check indices
    #RectBivariateSpline supposed to be faster? (but only for grid...)
    #f = interpolate.RectBivariateSpline(np.arange(array.shape[0]), np.arange(array.shape[0]), array)
    f = interpolate.interp2d(np.arange(array.shape[0]), np.arange(array.shape[0]), array)
    return np.array([f(i, j) for i, j in zip(x, y)]).reshape(len(x))

def mosaique_iris(alpha, delta, direc):

    band = 4
    equinox = 2000

    #the IDL code assumes square arrays, I think
    totsize = len(alpha)
    
    #sysco = 1 #celestial
    sysco = 2 #galactic
    #get_cootype(astr)? I think not needed
    if sysco == 2:
        equinox = 1950
    if equinox == 1950:
        fk4 = 1
    else:
        fk4 = 0
    if sysco == 2:
        print("convert coordinates from Galactic to Celestial B1950")
        alpha, delta = euler(alpha, delta, select=2, fk4=fk4)
    elif sysco == 3:
        print("not implemented")
        #print("convert coordinates from Ecliptic to Celestial")
        #euler(alpha, delta, select=4, fk4=fk4) #how convert?
    if equinox == 2000:
        print("not implemented")
        #print("precess coordinates from J2000 to B1950")
        #mprecess(alpha, delta, 2000, 1950)
        
    nan2undef(alpha, undef=-32768)
    nan2undef(delta, undef=-32768)
        
    result = np.zeros(totsize)
    weight = np.zeros(totsize) #not sure

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
    
    for i in range(nb):
        if (c1min >= ramin[i] and c1min <= ramax[i] and c2min >= demin[i] and c2min <= demax[i]):
            id_good[i] = 1
        if (c1min >= ramin[i] and c1min <= ramax[i] and c2max >= demin[i] and c2max <= demax[i]):
            id_good[i] = 1
        if (c1max >= ramin[i] and c1max <= ramax[i] and c2max >= demin[i] and c2max <= demax[i]):
            id_good[i] = 1
        if (c1max >= ramin[i] and c1max <= ramax[i] and c2min >= demin[i] and c2min <= demax[i]):
            id_good[i] = 1
        if (ramin[i] >= c1min and ramin[i] <= c1max and demin[i] >= c2min and demin[i] <= c2max):
            id_good[i] = 1
        if (ramax[i] >= c1min and ramax[i] <= c1max and demin[i] >= c2min and demin[i] <= c2max):
            id_good[i] = 1
        if (ramin[i] >= c1min and ramin[i] <= c1max and demax[i] >= c2min and demax[i] <= c2max):
            id_good[i] = 1
        if (ramax[i] >= c1min and ramax[i] <= c1max and demax[i] >= c2min and demax[i] <= c2max):
            id_good[i] = 1

    ind = np.where(id_good > 0)[0]
    nbind = ind.shape[0]
    if nbind == 0:
        print("No ISSA map corresponds to the header given")
        return -1
    print(nbind, "Issa maps will be combined to produce the mosaic")
    print("PLATE NUMBERS: ", inum[ind])
    
    for i in range(nbind):
        mapi, hi = get_iris(inum[ind[i]], direc=direc, band=band)
        #print("header:")
        #print(hi)
        print("file number: ", inum[ind[i]])

        '''
        #TEMP
        #try: fits.setval(fits_file, 'CDELT3', value=1), inside get_iris
        w = wcs.WCS(hi) #, filei . . . 
        #filei.close()
        xi, yi, trash = w.all_world2pix(alpha, delta, np.zeros(len(alpha)), 1) #origin 1 for FITS, zeros b/c axis3 has size 1
        print(xi, yi)
        '''

        #my own function:
        xi, yi = my_ad2xy(alpha, delta, hi)
        
        #xi, yi = ad2xy(alpha, delta, hi)
        tempo = mbilinear(xi, yi, mapi)
        
        tempo[np.isnan(tempo)] = -32768 #TESTING
        
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
    nindempty = indempty.shape[0]
    nbindw = indw.shape[0]

    if nbindw > 0:
        result[indw] = result[indw] / weight[indw]
    if nindempty > 0:
        result[indempty] = -32768

    with open('IRIS_values_skyphot.dat', 'w') as f:
        for i in range(totsize):
            f.write(str(result[i]) + "\n")

    return result



start_idx = 0
end_idx = 60

a_b = np.load("/Users/blakechellew/Documents/DustProject/l_b_original.npy")
a = a_b[start_idx:end_idx,0]
b = a_b[start_idx:end_idx,1]
result = mosaique_iris(a, b, '/Users/blakechellew/Documents/DustProject/IRIS/IRISNOHOLES_B4H0')
#print(result)

i100_original = np.load("/Users/blakechellew/Documents/DustProject/i100_1d.npy")[start_idx:end_idx]

#print("original:")
#print(i100_original[:50])

from scipy.stats import linregress
print(linregress(i100_original, result))

plt.plot(i100_original, result, 'k.', markersize=1)
#consider density plot
x = np.arange(np.min(result), np.max(result), 0.1)
plt.plot(x, x)
plt.xlabel("original i100")
plt.ylabel("IRIS i100")
plt.show()
