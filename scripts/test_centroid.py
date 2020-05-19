import numpy as np

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import wcs

import simulate as sim
import simulate_tools as st
import centroid_statistics as cs

import sys
if len(sys.argv) != 2:
    print("test_centroid.py takes one argument: the catalog number you want to illustrate")
    exit()

catnum = int(sys.argv[1])

import pprint
pp = pprint.PrettyPrinter(indent=4)

table_path = "../prepared_data.fits"
t = Table.read(table_path)

mask = (t['Name'] == 'catalogue_'+str(catnum))

ra, dec = t[mask]["redMaPPer_ra"][0]*u.degree, t[mask]["redMaPPer_dec"][0]*u.degree
x_ra, x_dec = t[mask]['x_ray_peak_ra'][0]*u.degree, t[mask]['x_ray_peak_dec'][0]*u.degree
radius = 0.15*t[mask]["R500_radius"][0]*u.kpc
z = t[mask]["Z"][0]
img_path = t[mask]['Data_Path'][0]

x = {
        "BCG ra"        : ra, 
        "BCG dec"       : dec, 
        "x-ray ra"      : x_ra, 
        "x-ray dec"     : x_dec,
        "radius"        : radius, 
        "redshift"      : z, 
        "data path"     : img_path 
    }

pp.pprint(x)

res, pix_res = sim.NSim(x_ra, x_dec, radius, z, img_path, 100, True, True)

xx = SkyCoord(x_ra, x_dec, frame='icrs')
a = SkyCoord(ra, dec, frame='icrs')

s = a.separation(xx)
print("BCG to x-ray sep:", st.theta_to_distance(z, s))

for ra, dec in zip(res[0], res[1]):
    simul = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
    sep = a.separation(simul)
#    print(ra, dec, st.theta_to_distance(z, sep))


ra, dec = t[mask]["redMaPPer_ra"][0]*u.degree, t[mask]["redMaPPer_dec"][0]*u.degree

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib as mpl 

#mpl.rc('font',family='Georgia', size=60)

img_path = t[mask]['Data_Path'][0]

with fits.open(img_path) as hdul:  # open a FITS file
    w = wcs.WCS(hdul[0].header)
    data = hdul[0].data
rra, rdec = w.all_world2pix(ra, dec, 1)
rx_ra, rx_dec = w.all_world2pix(x_ra, x_dec, 1)

pix_ra, pix_dec, pix_radius = st.convert_to_pixels(x_ra, x_dec, radius, z, img_path)
cut, xoff, yoff = st.circular_cut(data, (pix_ra, pix_dec), pix_radius)

fig, ax1 = plt.subplots(1, 1)#, figsize=(40, 30))

im1=ax1.imshow(cut, norm=LogNorm(), cmap=plt.cm.get_cmap('viridis'), vmax=2.5)
#ax2.imshow(cut, norm=LogNorm(), cmap=plt.cm.get_cmap('Greens'), vmin=0.1)

#plt.imshow(data_y, data_x, norm=LogNorm())
ax1.plot(pix_res[1][50]-yoff, pix_res[0][50]-xoff, 'd', color='orangered', 
        markersize=10, label="Median Simulation")
ax1.plot(rx_dec-yoff, rx_ra-xoff, 'X', color= 'orangered', 
        markersize=10, label="First Order X-Ray Center")
ax1.plot(rdec-yoff, rra-xoff, 'o', color='orangered', 
        markersize=10, label="BCG")

ax1.set_ylabel("RA\n")
ax1.set_xlabel("\nDEC")

xlocs = ax1.get_xticks()
ylocs = ax1.get_yticks()

print(xlocs, ylocs)

new_labels_x = []
new_labels_y = []
for x, y in zip(xlocs, ylocs):
    coords = w.all_pix2world(x+xoff, y+yoff, 1)
    new_labels_x += [str(coords[0])[0:7]]
    new_labels_y += [str(coords[1])[0:7]]

new_labels_x.reverse()

ax1.set_xticklabels(new_labels_x)
ax1.set_yticklabels(new_labels_y)

ax1.set_facecolor('whitesmoke')

#cax = plt.axes([0.95, 0.05, 0.05,0.9 ])
tick_nums = [1, 1.2, 1.4, 1.6, 1.8, 2, 2.0, 2.2, 2.4]
ticks = ["1", ' ', ' ' , ' ', "2",'', '',  "3", "4"]

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.08)

cbar = fig.colorbar(im1, ticks=tick_nums, cax=cax)
cbar.ax.set_yticklabels(ticks)

ax1.legend()

#fig.savefig("199_illustration.png")

plt.show()




















