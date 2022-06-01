
from astropy.table import Table

import numpy as np

import simulate_tools as st

from astropy import units as u
from astropy.coordinates import SkyCoord

def find_error(separations, z):

    median = np.median(separations)
    16perc = int(len(separations)*0.16)-1
    84perc = int(len(separations)*0.84)-1
    lower_one_sigma = np.abs(median - separations[16perc])
    upper_one_sigma = np.abs(separations[84perc] - median)
    return lower_one_sigma, upper_one_sigma 

