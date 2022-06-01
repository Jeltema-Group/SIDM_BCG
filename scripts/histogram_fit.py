
from scipy.stats import lognorm
from scipy.optimize import curve_fit

import numpy as np
import matplotlib as mpl
mpl.rc('font',family='Georgia', size=70)
import matplotlib.pyplot as plt

from astropy.table import Table

import sys


if len(sys.argv) != 2:
    print("histogram_fit.py takes one argument: the fits table of results from final_results.py")
    exit()

results_path = "../results/results_"+sys.argv[1]+".fits"

print("Reading results_"+sys.argv[1]+".fits ...")
t = Table.read(results_path)

print("Cutting out offsets greater than 40 kpc...")
mask = (t['BCG_offset_kpc']<40)
offsets = t['BCG_offset_kpc'][mask]
p_dist = offsets/np.sum(offsets)

n, bins, p = plt.hist(offsets, bins=3, 
                        color="tomato", ec="tomato", alpha = 0.5, lw=10,
                        zorder=1)

for item in p:
    item.set_height(item.get_height()/sum(n))

x = [(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]

from scipy.stats import lognorm
def f(x, s, u, a):
    A = [a*(s * xi * np.sqrt(2*np.pi))**(-1) for xi in x]
    f = [np.exp(-(np.log(xi)**2-np.log(u)**2)/(np.sqrt(2)*s)**2) for xi in x]
    return [A[i]*f[i] for i in range(len(A))]
    
sigma = []
for counts in n:
    if counts == 0:
        sigma += [4.0]
    else:
        sigma += [1/counts]

popt, pcov = curve_fit(f, x, n/np.sum(n), sigma=sigma)

x = [i/100 for i in range(1, 10000)]
print("s:", *popt)

print("Plotting the histogram to histogram_fit_" + sys.argv[1] + '.png')
plt.plot(x, f(x, *popt), color="tomato", linewidth=10, 
            label="$\mu$="+str(popt[0])[:5]+" $\pm$ "+ str(popt[1])[:4]+ " kpc")

plt.ylim(0, 0.6)

fig = plt.gcf()
fig.set_size_inches(40, 20)

plt.xlabel("BCG/Dark Matter Offset (kpc)")
plt.ylabel("p(Offset)")

plt.xticks(fontsize=60)
plt.yticks(fontsize=60)
plt.legend()

fig.savefig('../results/results_images/histogram_fit_' + sys.argv[1] + '.png')



































