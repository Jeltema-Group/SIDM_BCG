import numpy as np
import sys
from astropy.table import Table
from astropy import wcs
from astropy.io import fits

import simulate as sim
import simulate_tools as st
import os

from astropy import units as u


from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm


#input checking
if len(sys.argv)<3 or len(sys.argv)>4: 
    print("get_centroid.py takes the following arguments\n \
            \t(1) number of trials\n \
            \t(2) simulation name \n \
            \t(optional) convergence radius in pixels")
    exit()


N = int(sys.argv[1])
sim_name = sys.argv[2]
info_path = "../prepared_data.fits"
convergence = 2 #the center of the R500 cutout converges to 2 pixels
if len(sys.argv)==4: #unless otherwise specified from the user
    convergence = int(sys.argv[3])
    print("this is it:", convergence)

cluster_info = Table.read(info_path)

if not os.path.exists("../simulations/"+ sim_name + "/"):
    os.mkdir("../simulations/" + sim_name)

for i in range(len(cluster_info)):
    catnum = cluster_info['Name'][i][10:]
    print("\n\ncluster catalog_"+catnum)
    if os.path.exists("../simulations/"+ sim_name + "/simulated_centroids_" + catnum + ".fits"):
        continue
    elif cluster_info['Name'][i] == "catalogue_169":
        continue
    ra = cluster_info['redMaPPer_ra'][i]*u.degree
    dec = cluster_info['redMaPPer_dec'][i]*u.degree
    x_ra = cluster_info['x_ray_peak_ra'][i]*u.degree
    x_dec = cluster_info['x_ray_peak_dec'][i]*u.degree
    radius = 0.15*cluster_info['R500_radius'][i]*u.kpc
    Z = cluster_info['Z'][i]
    img_path = cluster_info['Data_Path'][i]

    with fits.open(img_path) as hdul:  # open the FITS file
        data = hdul[0].data  # assume the first extension is an image
        w = wcs.WCS(hdul[0].header)

    pix_ra, pix_dec, pix_radius = st.convert_to_pixels(x_ra, x_dec, radius, Z, img_path)
    cut, xoff, yoff = st.circular_cut(data, (pix_ra, pix_dec), pix_radius)
    
    current_peak = [int(pix_ra-xoff), int(pix_dec-yoff)]
    next_peak = list(st.centroid(cut))
    next_peak = [int(next_peak[1]), int(next_peak[0])]

    while st.pixel_distance(current_peak, next_peak) > convergence: 
        current_peak = next_peak ; i+=1 
        cut, xoff, yoff = st.circular_cut(data, (current_peak[1]+xoff, current_peak[0]+yoff), pix_radius)
        next_peak = list(st.centroid(cut))
        next_peak = [int(next_peak[1]), int(next_peak[0])]

    next_peak[0] += yoff ; next_peak[1] += xoff
    temp = next_peak[0] ; next_peak[0] = next_peak[1] ; next_peak[1] = temp
    print("old peak: ({}, {}) \
            \nnew peak: ({}, {})".format(pix_ra, pix_dec, next_peak[0], next_peak[1]))

    new_x_ra, new_x_dec = w.all_pix2world(int(next_peak[0]), int(next_peak[1]), 1)
    x, pix = sim.NSim(new_x_ra, new_x_dec, radius, Z, img_path, N, debug=True)

    sim_coords = Table()
    sim_coords["RA"] = x[0]
    sim_coords["DEC"] = x[1]
    sim_coords['pix_RA'] = pix[0]
    sim_coords['pix_DEC'] = pix[1]
    sim_coords["separations_kpc"] = sim.separations(ra, dec, x, Z)

    sim_coords.write("../simulations/" + sim_name + "/simulated_centroids_" + catnum + ".fits")













