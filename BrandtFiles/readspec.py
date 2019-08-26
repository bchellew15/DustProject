#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from astropy.io import fits

def loadspec(specfile, i):
    
    # log_10 of starting wavelength in Angstroms
    lam0 = 10**fits.open(specfile)[2].data[i]

    # flux density
    flambda = fits.open(specfile)[0].data[i]

    # inverse variance
    ivar = fits.open(specfile)[1].data[i]

    # Definition of the wavelength array: 1e-4 spacing in log_10(lambda)
    lam = lam0*10**(1e-4*np.arange(len(ivar)))

    return [lam, flambda, ivar]


def savespec(lam, flambda, ivar, outfile=None):

    specarr = np.zeros((len(lam), 3))
    specarr[:, 0] = lam
    specarr[:, 1] = flambda
    specarr[:, 2] = ivar

    if outfile is None:
        return specarr
    else:
        np.savetxt(outfile, specarr, fmt="%.5g")
        return


if __name__ == "__main__":

    specfile = 'skyfibers_nativelam.fits'

    # Save the 101st spectrum (index 100 from 0) as an example
    i = 100
    lam, flambda, ivar = loadspec(specfile, i)
    savespec(lam, flambda, ivar, 'spec_%06d.dat' % (i))
    
    # flambda is the excess sky spectrum-it should be nearly zero
    # apart from measurement errors.  Quantify this with chi squared:

    chisq_dof = np.sum(flambda**2*ivar)/len(lam)
    print('chi2 per degree of freedom of the sky spectrum: %.2f' % chisq_dof)
