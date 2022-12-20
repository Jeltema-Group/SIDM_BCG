
import numpy as np
import os
import sys

from astropy.table import Table
from astropy.table import Column

import centroid_statistics as cs
import simulate_tools as st

from astropy import units as u

if len(sys.argv) != 2:
    print("INPUT ERROR: final_results.py requires one argument: \n \
                the name of the simulation")
    exit()
results_path = "../simulations/" + sys.argv[1] + "/"

#determine centroid and associated error
cluster_info = Table.read("../simulations/"+sys.argv[1]+"_prepared_data.fits")

name = np.zeros(len(os.listdir(results_path)))
center_to_BCG = np.zeros(len(os.listdir(results_path)))
up_errors = np.zeros(len(os.listdir(results_path))) ; i = 0
lo_errors = np.zeros(len(os.listdir(results_path)))

print("Calculating median separtation and errors ...")
cluster_info.sort(['Name'])
for cluster in os.listdir(results_path):
    name[i] = int(cluster[20:-5]) 
    path = results_path+cluster

    simulations = Table.read(path)
    separations = simulations['separations_kpc']
    separations.sort()
    center_to_BCG[i] = np.median(separations)

    #print(cluster, len(separations))

    up_errors[i], lo_errors[i] = cs.find_error(separations, cluster_info['Z'][i])

    i+=1

#make a final table
print("Writing results_" + sys.argv[1] + ".fits ...")
n = Column(name, name="catalog_number", dtype=np.int32)
o = Column(center_to_BCG, name="BCG_offset_kpc")
u = Column(up_errors, name="Upper_Error_kpc")
l = Column(lo_errors, name="Lower_Error_kpc")
z = cluster_info['Z']
results = Table([n, o, u, l, z])
result_table_path = "../results/results_" + sys.argv[1] + ".fits"
if os.path.exists(result_table_path):
    os.remove(result_table_path)
results.write(result_table_path)

#plot results
print("Plotting the results to result_" + sys.argv[1] + '.png')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
mpl.rc('font',family='Georgia', size=90)
mpl.rcParams.update({'errorbar.capsize': 20})


fig = mpl.pyplot.gcf()
fig.set_size_inches(80, 40)

l=len(results)
errs = np.array([results['Lower_Error_kpc'],results['Upper_Error_kpc']])
plt.errorbar(np.arange(l), results['BCG_offset_kpc'], fmt='o', color="g", 
        yerr=errs, linewidth=5, markersize=20, elinewidth=5,
        markeredgewidth=5)

plt.plot([-2, l+2], [0, 0], linestyle=':', linewidth=5, color='tomato')
plt.xlim(-0.5, l+0.5)

plt.xticks(np.arange(l), results['catalog_number'], rotation='vertical')
#plt.yticks(fontsize=70)
plt.xlabel("Catalog Numbers")
plt.ylabel("Separations (kpc)\n")

ax = plt.gca()
ax.set_facecolor('whitesmoke')

fig.savefig("../results/results_images/result_" + sys.argv[1] + '.png')
































