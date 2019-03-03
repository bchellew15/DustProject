import numpy as np
from scipy.interpolate import interp1d

def f99(wave, ebv=1, r_v=3.1):

    """

    Returns optical depth at an input wavelength or wavelengths as a
    function of E(B-V) (i.e. dust column density) and R_V (a
    parametrization of dust absorption properties).

    """
    
    x = 1.e4 / np.ravel(wave)
    model = 'f99'

    if np.any(x < 0.167) or np.any(x > 11.):
        raise ValueError('Wavelength(s) must be between 910 A and 6 um')
    if model == 'fm07' and abs(r_v - 3.1) > 0.001:
        raise ValueError('fm07 model not implementend for r_v != 3.1')

    k = np.zeros_like(x)
    uv_region = (x >= 1.e4 / 2700.)
    oir_region = ~uv_region

    # UV region
    y = x[uv_region]
    if model == 'f99':
        x0, gamma = 4.596, 0.99
        c3, c4, c5 = 3.23, 0.41, 5.9
        c2 = -0.824 + 4.717 / r_v
        c1 = 2.030 - 3.007 * c2
        d = y**2 / ((y**2 - x0**2)**2 + y**2 * gamma**2)
        f = np.zeros_like(y)
        valid = (y >= c5)
        f[valid] = 0.5392 * (y[valid] - c5)**2 + 0.05644 * (y[valid] - c5)**3
        k_uv = c1 + c2 * y + c3 * d + c4 * f
    if model == 'fm07':
        x0, gamma = 4.592, 0.922
        c1, c2, c3, c4, c5 = -0.175, 0.807, 2.991, 0.319, 6.097
        D = y**2 / ((y**2-x0**2)**2 + y**2 * gamma**2)
        k_uv = np.zeros_like(y)
        valid = (y <= c5)
        k_uv[valid] = c1 + c2*y[valid] + c3*D[valid]
        valid = (y > c5)
        k_uv[valid] = c1 + c2*y[valid] + c3*D[valid] + c4*(y[valid] - c5)**2
    k[uv_region] = k_uv

    # Calculate values for UV spline points to anchor OIR fit
    x_uv_spline = 1.e4 / np.array([2700., 2600.])
    d = (x_uv_spline**2 /
         ((x_uv_spline**2 - x0**2)**2 + x_uv_spline**2 * gamma**2))
    k_uv_spline = c1 + c2 * x_uv_spline + c3 * d

    # Optical / IR region
    y = x[oir_region]
    if model == 'f99':
        anchors_x = 1.e4 / np.array([np.inf, 26500., 12200., 6000., 5470.,
                                     4670., 4110.])

        # The OIR anchors are from IDL astrolib, not F99.
        anchors_extinction = np.array(
            [0.,
             0.26469 * r_v / 3.1,  # IR
             0.82925 * r_v / 3.1,  # IR
             -0.422809 + 1.00270 * r_v + 2.13572e-04 * r_v**2,  # optical
             -5.13540e-02 + 1.00216 * r_v - 7.35778e-05 * r_v**2,
             0.700127 + 1.00184 * r_v - 3.32598e-05 * r_v**2,
             (1.19456 + 1.01707 * r_v - 5.46959e-03 * r_v**2 +
              7.97809e-04 * r_v**3 - 4.45636e-05 * r_v**4)]
            )

        anchors_x = np.append(anchors_x, x_uv_spline)
        anchors_k = np.append(anchors_extinction - r_v, k_uv_spline)

    oir_spline = interp1d(anchors_x, anchors_k, kind='cubic')
    k[oir_region] = oir_spline(y)

    # notes: E(B-V) = 1 for I_100 = 50 MJy/sr at 18 K
    # E(B-V) = 1 for I_100 = 31 MJy/sr at 17 K

    return ebv*(k + r_v)
