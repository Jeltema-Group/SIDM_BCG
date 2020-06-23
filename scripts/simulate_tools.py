
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
    for row in noisey_array:
        for pixel in range(len(row)): 
            noise = np.random.poisson(int(row[pixel]))
            row[pixel] = row[pixel] + noise
    return noisey_array




'''Centroid Measurement Function'''

def centroid(data):
    wx = 0 ; wy = 0 ; w = 0
    for i in range(len(data)):
        for j in range(len(data[0])):
            wx += (i+1)*data[i][j]
            wy += (j+1)*data[i][j]
            w += data[i][j]
    
    if w == 0:
        return 0., 0.
    return wx/w, wy/w



'''Data Preparation Functions'''

def circular_cut(arr, center, radius):
    x0=center[0]-radius ; y0=center[1]-radius
    cut = np.zeros((int(2*radius), int(2*radius)))
    for x in range(2*radius):
        for y in range(2*radius):
            if (x-radius)**2 + (y-radius)**2 < radius**2:
                if arr[y+y0][x+x0] is np.nan:
                    print(arr[y+y0][x+x0])
                cut[x][y] = 0 if arr[y+y0][x+x0]==np.nan else arr[y+y0][x+x0]
    return np.array(cut), x0, y0

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
    




