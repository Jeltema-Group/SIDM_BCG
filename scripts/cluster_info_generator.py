from astropy.table import Table
import os
import sys

if len(sys.argv) != 3:
    print("INPUT ERROR: cluster_info_generator.py takes two arguments:")
    print("\t(1): data table with Name, redMaPPer ra and dec, R500 radius, redshift, x-ray peak ra and dec")
    print("\t(2): path where the x-ray images are stored")
    exit()

des_data_path = "../data/" + sys.argv[1]
images_path = sys.argv[2]

data = Table.read(des_data_path)

new_data = Table()
new_data['Name'] = data['Name']
new_data['redMaPPer_ra'] = data['redMaPPer_ra']
new_data['redMaPPer_dec'] = data['redMaPPer_dec']

for row in range(len(data)):
    if data['redmapper wrong'][row]:
        new_data['redMaPPer_ra'][row] = data['right RA'][row]
        new_data['redMaPPer_dec'][row] = data['right DEC'][row]

new_data['R500_radius'] = data['r500_radius']
new_data['Z'] = data['Redshift']

new_data['x_ray_peak_ra'] = data['x_ray_peak_ra']
new_data['x_ray_peak_dec'] = data['x_ray_peak_dec']

paths = ["xxx" for x in range(len(new_data))] ; i = 0

for cluster in new_data['Name']:
    num = cluster[10:]
    match = [s for s in os.listdir(images_path) if (num+'_') in s]
    match = [x for x in match if x[:len(num)+1] == (num+'_')]
    if len(match) == 1:
        paths[i] = images_path + match[0]
    else:
        print("issue with catalogue number", num)
    i += 1
new_data['Data_Path'] = paths


new_data.write("../prepared_data.fits")


