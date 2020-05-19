
import sys
import simulate_tools as st 

from astropy.io import fits
from astropy import wcs
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

def NSim(ra, dec, radius, Z, img_path, N, 
        save_sample=False, ret_pix_vals=True, debug=False):
    if debug:
        print("\nAnalyzing image", img_path, '\n') 

    # open image data
    with fits.open(img_path) as hdul:  # open the FITS file
        data = hdul[0].data  # assume the first extension is an image
        w = wcs.WCS(hdul[0].header)

    # convert nans to 0.0:
    where_are_NaNs = np.isnan(data)
    data[where_are_NaNs] = 0.0

    # get pixel information 
    pix_ra, pix_dec, pix_radius = st.convert_to_pixels(ra, dec, radius, Z, img_path)

    # create circular cut
    cut, xoff, yoff = st.circular_cut(data, (pix_ra, pix_dec), pix_radius)

    new_centroids = np.zeros((2, N))
    pix_centroids = np.zeros((2, N))

    for sim in range(N):
        if sim%(int(N/10)) == 0 and debug:
            print("\t\t proccessing simulation", sim)
        new_sim = st.add_noise(cut)
        coords = st.centroid(new_sim)
        pix_centroids[0][sim], pix_centroids[1][sim] = coords[0]+xoff, coords[1]+yoff
        new_centroids[0][sim], new_centroids[1][sim] = w.all_pix2world(coords[0]+xoff, coords[1]+yoff, 1)

    if save_sample:
            np.save("example_simulation", new_sim, allow_pickle=True)

    if ret_pix_vals:
        return new_centroids, pix_centroids

    return new_centroids



def separations(BCG_ra, BCG_dec, results, z):
    
    seps = []
    a = SkyCoord(BCG_ra, BCG_dec, frame='icrs')

    units=[]
    for ra, dec in zip(results[0], results[1]):
        simul = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
        sep = a.separation(simul)
        dist = st.theta_to_distance(z, sep).to(u.kpc)
        seps += [dist.value]
    
    return seps





























