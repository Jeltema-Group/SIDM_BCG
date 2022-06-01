import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib as mpl
#mpl.rc('font',family='Georgia', size=45)
#mpl.rc('lines', markersize=50)

import os
import sys
import numpy as np

from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy import wcs

import simulate_tools as st


if len(sys.argv) != 3:
    print("INPUT ERROR: sim_illustrations.py takes two arguments:")
    print("\t(1): the catalog number of the cluster to illustrate")
    print("\t(2): the name of the simulation under the \'simulations\' directory")
    exit()


goal_catnum = sys.argv[1]
print("illustrating catalogue_"+goal_catnum)

cluster_info = Table.read("../simulations/"+sys.argv[2]+"_prepared_data.fits")

print(cluster_info.colnames)

sims_path = "../simulations/"+sys.argv[2]+"/"
simulations_dir = os.listdir(sims_path)

for sim in simulations_dir:
    catnum = sim[20:-5]
    if catnum != goal_catnum:
        continue
    mask = (cluster_info['Name'] == "catalogue_"+catnum)
    img_path = cluster_info[mask]['Data_Path'][0]
    img_data = fits.open(img_path)[0].data
    
    with fits.open(img_path) as hdul:  # open a FITS file
        w = wcs.WCS(hdul[0].header)

    radius= 0.15*cluster_info[mask]['R500_radius'][0]*u.kpc
    Z = cluster_info[mask]['Z'][0]
    
    sim_t = Table.read(sims_path+"simulated_centroids_"+ catnum +".fits")
    sim_t.sort('separations_kpc') 
    sim_ra, sim_dec = sim_t['pix_RA']*u.degree, sim_t['pix_DEC']*u.degree
    
    theta = st.distance_to_theta(Z, radius)
    pix_radius = st.degrees_to_pixels(theta)
   
    median_coords = (sim_t['pix_RA'][49], sim_t['pix_DEC'][49])

    ax = plt.gca() ; ax.set_facecolor('whitesmoke')
    fig = plt.gcf()
    #fig.set_size_inches(30, 30)

    #create circle
    x = [pix_radius*np.cos(t/100) + median_coords[0] for t in range(2*314)]
    y = [pix_radius*np.sin(t/100) + median_coords[1] for t in range(2*314)]

    #data
    ax.scatter(sim_ra, sim_dec, alpha=0.4, label="simulation centroids", 
                    color='mediumseagreen', edgecolor='green')#, s=6000)
    ax.plot([median_coords[0]], [median_coords[1]], 'X', color='navy', label="median")
   
    #redmapper coordinates
    red_ra, red_dec = w.all_world2pix(cluster_info[mask]['redMaPPer_ra'][0], 
                                        cluster_info[mask]['redMaPPer_dec'][0], 1)
    ax.plot([red_ra], [red_dec], 'X', color='orangered', label="redmapper coordinates")

    #r500 circle
    #plt.plot(x, y, label="0.15 * r500")
  
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    print(img_data)
    im = ax.imshow(img_data, norm=LogNorm(), cmap=plt.cm.get_cmap('binary'), vmax=10)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.5)

    cbar=plt.colorbar(im, cax=cax)

    ax.legend(loc=1)#, fontsize=70)

    #ax.set_xlim(median_coords[0]-0.02*pix_radius-0.5, median_coords[0]+0.02*pix_radius-0.5)
    #ax.set_ylim(median_coords[1]-0.02*pix_radius+0.5, median_coords[1]+0.02*pix_radius+0.5)

    xlocs = ax.get_xticks()
    ylocs = ax.get_yticks()

    print(xlocs, ylocs)

    new_labels_x = []
    new_labels_y = []
    for x, y in zip(xlocs, ylocs):
        coords = w.all_pix2world(x, y, 1)
        new_labels_x += [str(coords[0])[0:8]]
        new_labels_y += [str(coords[1])[0:8]]

    new_labels_x.reverse()

    ax.set_xticklabels(new_labels_x, rotation=45)
    ax.set_yticklabels(new_labels_y, rotation=45)

    ticks=['0', '15']
    cbar.ax.set_yticklabels(ticks)

    ax.set_xlabel("RA (degrees)")
    ax.set_ylabel("DEC (degrees)")

    #fig.savefig('catalog_'+catnum+'_illustration_zoomed.png')
    plt.show()

    exit()

print("catalogue_"+goal_catnum, "does not exist")





























