#copy/paste version of ad2xy from internet
'''
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
'''

#started to translate mbilinear from IDL
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
