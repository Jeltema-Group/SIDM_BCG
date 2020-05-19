
from astropy.table import Table

import numpy as np

import simulate_tools as st

from astropy import units as u
from astropy.coordinates import SkyCoord

def find_error(separations, z):

    median = np.median(separations)
    lower_one_sigma = np.abs(median - separations[15] )
    upper_one_sigma = np.abs(separations[83] - median)
    return lower_one_sigma, upper_one_sigma 

