
import numpy as np
import copy
from astropy.table import Table
from astropy.io import fits
from astropy import wcs

from astropy.coordinates import Distance
from astropy import units as u


'''
LIBRARY TOOLS DOCUMENTATION:

    Noise Generation Functions:
        - random_poisson
            - arguments: 
                - lam: expected value of the poisson distribution (lambda)
                - N: 
            - output: random number modulated by the corresponding poisson function
        - add_noise
            - arguments: 
                - data: input image data as a 2-D array
                - lam: expected value of the poisson distribution (lambda)
            - output: an array with noise added to the array

    Measurement Functions:
        - centroid
            - arguments: 
                - data: 2-D array of the original image
            - output: 
                - new 2-D array with the same dimensions but noise has been added


    Data Preparation Functions:
        - circular_cut
            - arguments: 
                - arr: 2-D array of the data
                - center: center of the desired circle
                - radius: radius of the desired circle
            - output: 
                - 2-D array of dimensions 2*radius x 2*radius which contains a cut out
                  circle of the data around the center coordinates
        - convert_to_pixels
            - arguments:
                - RA: RA of the BCG
                - DEC: Dec of the BCG
                - radius: R500 radius 
                - redshift: redshift of the cluster
                - img_path: path to the image data
            - output:
                - the pixel position of the BCG center (ra, dec)
                - the pixel radius of the circular cut
        - convert_to_radec
            - arguments:
                - ra: RA position in pixels
                - dec: DEC position in pixels
                - w: World Coordinate System object for image
            - output: 
                - the cannonical RA and DEC coordinates in degrees
        - distance_to_theta
            - arguments:
                - redshift: redshift of the cluster
                - radius: physical separation (kpc)
            - output: 
                - angular separation (radians)
        - theta_to_distance
            - arguments:
                - redshift: redshift of the cluster
                - theta: angular separation (radians)
            - output: 
                - physical separation (in kpc)
        - degrees_to_pixels
            - arguments:
                - theta: angular separation (radians)
            - output: 
                - number of pixels corresponding to this separation
        - pixels_to_degrees
            - arguments:
                - pixels: number of pixels for separation
            - output: 
                - angular separation (radians)
'''



'''Noise Generation Functions'''

def add_noise(data):  
    noisey_array = copy.deepcopy(data)
    noisey_array += np.random.poisson(np.int_(noisey_array))
    return noisey_array


'''Centroid Measurement Function'''

def centroid(data):
    w = np.sum(data)

    x_pos = np.reshape(np.int_(np.linspace(1, len(data), len(data))), (len(data),1))
    y_pos = np.int_(np.linspace(1, len(data[0]), len(data[0])))
    
    wx = np.sum(x_pos*data)
    wy = np.sum(y_pos*data)
    
    if w == 0:
        return 0., 0.
    return wx/w, wy/w


'''Data Preparation Functions'''

def circular_cut(arr, center, radius):
    x0=center[0]-radius ; y0=center[1]-radius
    #The square of the image we're looking at.
    cut = np.transpose(arr[y0:y0+2*radius])[x0:x0+2*radius]
    
    #A circle nestled neatly in the square 
    x_pos = np.reshape(np.int_(np.linspace(0, len(cut)-1, len(cut))), (len(cut),1))
    y_pos = np.int_(np.linspace(0, len(cut[0])-1, len(cut[0])))
    circ = (x_pos-radius)**2 + (y_pos-radius)**2
    circ_true=circ < radius**2

    #Where the NaNs are
    nan_true = cut != np.nan
    
    #Zeroing out the NaNs and stuff outside the circle in the cut.
    return nan_true*circ_true*cut, x0, y0

from astropy.coordinates import SkyCoord

def convert_to_pixels(RA, DEC, radius, redshift, img_path):

    with fits.open(img_path) as hdul:  # open a FITS file
        w = wcs.WCS(hdul[0].header)

    ra_pix, dec_pix = w.all_world2pix(RA, DEC, 1)

    ra_pix = int(ra_pix)
    dec_pix = int(dec_pix)
    
    theta = distance_to_theta(redshift, radius)
    pix_radius = degrees_to_pixels(theta)

    return ra_pix, dec_pix, pix_radius

#returns theta in radians
def distance_to_theta(redshift, radius):
    R = Distance(z=redshift).to(u.kpc)
    try:
        radius = radius.to(u.kpc)
    except AttributeError:
        radius = radius*u.kpc
    return 2*np.arcsin(radius/(2*R)) # RADIUS MUST HAVE ASTROPY UNITS
    
#returns distance in kpc
def theta_to_distance(redshift, theta): 
    R = Distance(z=redshift).to(u.kpc)
    return R*2*np.sin(theta/2)
    
#converts angle of separation to a distance using Chandra's pixel to angle conversion 
def degrees_to_pixels(theta): #conversion is 0.492 arc seconds per pixel
    return int((theta.to(u.arcsec)/(0.492*(u.arcsec))))

def pixels_to_degrees(pixels):
    return (((0.492)*pixels)*u.arcsec).to(u.degree)


def pixel_distance(coords1, coords2):
    x1, y1 = coords1
    x2, y2 = coords2
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)
    




