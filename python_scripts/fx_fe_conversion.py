# Conversion of x-ray form factors to electron scattering factors
# Constants and conventions adopted from Lobato et al. https://doi.org/10.1107/S205327331401643X


import numpy as np


def mott_bethe(x, f_x, Z):

    a_0 = 0.529177210544
    x = np.asarray(x, dtype=float)
    f_x = np.asarray(f_x, dtype=float)
    f_e = (1.0 / (2.0 * np.pi**2 * a_0)) * ((Z - f_x) / (x**2))
    
    return f_e
































