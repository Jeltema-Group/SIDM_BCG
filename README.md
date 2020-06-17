# SIDM

There have been many simulations that show that observing an offset between x-ray 
centers and the Brightest Cluster Galaxy (BCG) of a galaxy cluster can give 
competitive constraints on the self-interaction of dark matter. 

In order to determine the x-ray center of a galaxy, we must determine how much of
an effect noise has on the centroid measurement. This repository holds the code that 
creates noise simulations on x-ray data, remeasures the centroid, and determines 
the best estimate of the x-ray center.

## How to run the simulations on a small Dark Energy Survey set of galaxies:

```
cd scripts/
```

### First put the data into a standard format. 

```
python cluster_info_generator.py DES_data_mergers_removed.fits ../data/xray_data/
```

First argument (in this case ``DES_data_mergers_removed.fits``) is the data from 
DES. The columns we use from this table are:
1. ``Name`` (in 'catalogue_x' format)
2. ``redMaPPer_ra`` and ``redMaPPer_dec``
3. ``reMaPPer wrong`` (which tells us if we need to correct the ra and dec)
4. ``r500_radius``
5. ``Redshift``
6. ``x_ray_peak_ra`` and ``x_ray_peak_dec``

Second argument (in this case ``../data/x-ray_data/``) is the directory that 
contains all of the x-ray data from Chandra.

This will output a file called ``prepared_data.fits``. Do not move or alter this table. 

### We are now ready to run the simulations
```
python get_centroids.py 100 first_simulation
```
First argument (``100``) is the number of simulations you want to run per cluster.

Second argument (``first_simulation``) is the name of your simulation. This will 
create a directory in ~/simulations/ that stores the results of this step.

This step will take a while, but there should be annotations on what simulation and 
cluster is being worked on.

## Analyzing the simulations

### Compile the results of the simulations into one fits file
```
python final_results.py first_simulation
```
This will create an image in the ``~/results/results_images/`` directory called ``result_first_simulation.png``. 

### Make the offset histogram and calculate the predicted cross-section
```
python histogram_fit.py first_simulation
```
This will create another image in the ``~/results/results_images/`` directory called ``histogram_fit_first_simulation.png``.

### Check individual clusters
```
python sim_illustrations.py 199 first_simulation
```
