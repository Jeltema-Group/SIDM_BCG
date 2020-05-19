


import numpy as np
import sys
from astropy.table import Table
from astropy import wcs
from astropy.io import fits

import simulate as sim
import simulate_tools as st
import os

from astropy import units as u


N = 100
info_path = "prepared_data.fits"
convergence = sys.argv[1]
print("convergence to:", convergence, "kiloparsecs")

cluster_info = Table.read(info_path)


if not os.path.exists("convergence_sims/"):
    os.mkdir("./convergence_sims")

for i in range(len(cluster_info)):
    catnum = cluster_info['Name'][i][10:]
    print("\n\ncluster catalog_", catnum)
    if os.path.exists("convergence_sims/simulated_centroids_" + catnum + ".fits"):
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
    
    current_peak = (pix_ra, pix_dec)
    next_peak = list(st.centroid(cut))
    next_peak[0] += yoff ; next_peak[1] += xoff
    next_peak = (next_peak[1], next_peak[0])
    sep = st.distance_to_theta(Z, convergence)
    pix_sep = st.degrees_to_pixels(sep)
    while st.pixel_distance(current_peak, next_peak) > pix_sep:
        current_peak = next_peak
        cut, xoff, yoff = st.circular_cut(data, (pix_ra, pix_dec), pix_radius)
        next_peak = list(st.centroid(cut))
        next_peak[0] += yoff ; next_peak[1] += xoff
        next_peak = (next_peak[1], next_peak[0])

    print("\told:", (pix_ra, pix_dec), 
          "\n\tnew:", (int(next_peak[0]), int(next_peak[1])))

    new_x_ra, new_x_dec = w.all_pix2world(int(next_peak[0]), int(next_peak[1]), 1)
    x, pix = sim.NSim(new_x_ra, new_x_dec, radius, Z, img_path, N, debug=True)

    sim_coords = Table()
    sim_coords["RA"] = x[0]
    sim_coords["DEC"] = x[1]
    sim_coords['pix_RA'] = pix[0]
    sim_coords['pix_DEC'] = pix[1]
    sim_coords["separations_kpc"] = sim.separations(ra, dec, x, Z)

    sim_coords.write("./convergence_sims/simulated_centroids_" + catnum + ".fits")













