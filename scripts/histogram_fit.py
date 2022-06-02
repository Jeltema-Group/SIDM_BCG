
from scipy.stats import lognorm
from scipy.optimize import curve_fit

import numpy as np
import matplotlib as mpl
mpl.rc('font',family='Georgia', size=70)
import matplotlib.pyplot as plt

from astropy.table import Table

import sys


if (len(sys.argv) < 2) or (len(sys.argv) >= 5):
    print('''histogram_fit.py takes up to three arguments: 
    1)  The name of the fits table of results from final_results.py.
    2)  The maximum offset to include. Optional; default is 40.
    3)  The number (integer) of bins for the histogram. Optional; default is equal to 5.''')
    exit()

if len(sys.argv) == 2:
    sysargvs=[sys.argv[1], 40, 5]
elif len(sys.argv) == 3:
    sysargvs=[sys.argv[1], float(sys.argv[2]), 5]
elif len(sys.argv) == 4:
    sysargvs=[sys.argv[1], float(sys.argv[2]), int(sys.argv[3])]

print(f'''Inputs:
    Name of fits table for results: {sysargvs[0]}.fits
    Maximum offset: {sysargvs[1]} kpc
    Number of bins: {sysargvs[2]}
    ''')

results_path = "../results/results_"+sysargvs[0]+".fits"

print("Reading results_"+sysargvs[0]+".fits ...")
t = Table.read(results_path)

print(f"Cutting out offsets greater than {sysargvs[1]} kpc...")
mask = (t['BCG_offset_kpc']<sysargvs[1])
offsets = t['BCG_offset_kpc'][mask]
p_dist = offsets/np.sum(offsets)

n, bins, p = plt.hist(offsets, bins=sysargvs[2], 
                        color="tomato", ec="tomato", alpha = 0.5, lw=10,
                        zorder=1)

for item in p:
    item.set_height(item.get_height()/sum(n))


from scipy.stats import lognorm
params = lognorm.fit(offsets, floc=0)

print(f'''Lognorm fit values:
    mu: {params[2]}
    sigma: {params[0]}
    loc: {params[1]}
    loc should be 0.0.
    ''')

if max(offsets) <= 120:
    x = np.linspace(1/100, 120, num=12000)
else:
    ceilmax = np.ceil(max(offsets))
    x = np.linspace(1/100, ceilmax, num=100*ceilmax)

y = lognorm.pdf(x, np.ones(x.shape)*params[0], params[1], params[2])

print("Plotting the histogram to histogram_fit_" + sysargvs[0] + '.png')
plt.plot(x, y, color="tomato", linewidth=10, 
            label="$\mu$="+str(params[2])[:5]+" $\pm$ "+ str(params[0])[:4]+ " kpc")

if max(n)/sum(n) > max(y):
    plt.ylim(0, np.ceil(max(n)/sum(n)*10)/10)
else:
    plt.ylim(0, np.ceil(max(y)*10)/10)

#plt.ylim(0, 0.6)

fig = plt.gcf()
fig.set_size_inches(40, 20)

plt.xlabel("BCG/Dark Matter Offset (kpc)")
plt.ylabel("p(Offset)")

plt.xticks(fontsize=60)
plt.yticks(fontsize=60)
plt.legend()

fig.savefig('../results/results_images/histogram_fit_' + sysargvs[0] + '.png')
