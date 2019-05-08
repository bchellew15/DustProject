import numpy as np
import idlwrap
import pydl
from astropy.coordinates import SkyCoord
from astropy.io import fits
import glob
import os
import scipy.interpolate as interpolate

#next steps:
#finish implementing functions
#test functions independently
#better way to convert coordinates (find examples)
#upload to github


#should work
#if select == 2, galactic to RA / DEC
#fk4 == 1 means convert to B1950
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

#not needed
#implemented in IDL
def mprecess():
    return True

#works
def nan2undef(l, undef):
    l[np.isnan(l)] = undef

#implemented in IDL
def get_iris(num, direc, band, silent):
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
        imap = fits.open(result[0])[0].data[0]
        header = fits.open(result[0])[0].header
        print(header)
        #sxaddpar?
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

    print("shape: ", x.shape)
    
    #check indices
    f = interpolate.interp2d(np.arange(array.shape[0]), np.arange(array.shape[0]), array)
    return np.array([f(i, j) for i, j in zip(x, y)]).reshape(len(x))

def mosaique_iris(alpha, delta, direc):

    band = 4
    silent = 0
    equinox = 2000

    #the IDL code assumes square arrays?
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
        euler(alpha, delta, select=2, fk4=fk4) #how convert?
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
    weight = np.copy(result) #copy?

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
        mapi, hi = get_iris(inum[ind[i]], direc=direc, band=band, silent=silent)
        #astri = extast(hi)
        xi, yi = ad2xy(alpha, delta, hi)
        tempo = mbilinear(xi, yi, mapi)
        indw = np.where(tempo != -32768)[0]
        nbindw = indw.shape[0]
        if nbindw > 0:
              weight[indw] = weight[indw]+1 #check addition
              result[indw] = result[indw] + tempo[indw]
              
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

a = np.arange(9)
b = np.arange(9)
result = mosaique_iris(a, b, '/Users/blakechellew/Documents/DustProject/IRIS/IRISNOHOLES_B4H0')
print(result)
