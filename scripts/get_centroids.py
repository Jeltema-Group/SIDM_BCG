import numpy as np
import sys
from astropy.table import Table

import simulate as sim
import simulate_tools as st
import os

from astropy import units as u

#input checking
if len(sys.argv) != 3:
    print("get_centroid.py takes two arguments: number of trials and simulation name")
    exit()

table_path = "../prepared_data.fits"
N = int(sys.argv[1])

cluster_info = Table.read(table_path)
header = ['Name', 'redMaPPer_ra', 'redMaPPer_dec', 'R500_radius', 'Z', 'x_ray_peak_ra', 'x_ray_peak_dec', 'Data_Path']
if cluster_info.colnames != header: 
    print("TABLE READ ERROR: table header must have the following format:", header) 
    print("current input:", cluster_info.colnames)
    exit()

if os.path.exists("simulations/"+sys.argv[2]):
    print("FILE NAME ERROR:", "../simulations/"+sys.argv[2], "already exists")
else:
    print("new directory:", "../simulations/"+sys.argv[2])
    os.mkdir("../simulations/" + sys.argv[2])

#--- end input checking----

for i in range(len(cluster_info)):
    catnum = cluster_info['Name'][i][10:]
    if cluster_info['Name'][i] == "catalogue_169":
        continue
    ra = cluster_info['redMaPPer_ra'][i]*u.degree
    dec = cluster_info['redMaPPer_dec'][i]*u.degree
    x_ra = cluster_info['x_ray_peak_ra'][i]*u.degree
    x_dec = cluster_info['x_ray_peak_dec'][i]*u.degree
    radius = 0.15*cluster_info['R500_radius'][i]*u.kpc
    Z = cluster_info['Z'][i]
    img_path = cluster_info['Data_Path'][i]

    x, pix = sim.NSim(x_ra, x_dec, radius, Z, img_path, N, debug=True)

    sim_coords = Table()
    sim_coords["RA"] = x[0]
    sim_coords["DEC"] = x[1]
    sim_coords['pix_RA'] = pix[0]
    sim_coords['pix_DEC'] = pix[1]
    sim_coords["separations_kpc"] = sim.separations(ra, dec, x, Z)

    sim_coords.write("../simulations/" + sys.argv[2] + "/simulated_centroids_" + catnum + ".fits")




















