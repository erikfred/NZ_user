from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pickle
import sys, os
from pathlib import Path

sys.path.append(os.path.abspath('util'))
import Lfun
import zrfun
import zfun

for x in range(1,25):
    print(x)

# x = range(1,25)
# for n in x:
#     print(n)

# Ldir = Lfun.Lstart()
#
# out_dir = Ldir['parent'] / 'LO_output' / 'bpress_PM' # I think this assumes I'm working out of the bpress subfolder
# print(out_dir)

"""
# read layers.nc file (will eventually become loop over all dates)
dir1 = "../LO_data/cas6_v0_live/"
dir2 = dir1 + "f2016.12.15/"
dir3 = dir2 + "ocean_his_0001.nc"
ds1 = Dataset(dir3)

# print metadata
print(ds1.__dict__)
for dim in ds1.dimensions.values():
    print(dim)
for var in ds1.variables.values():
    print(var)
"""
"""
# establish file structure
workdir = os.path.dirname(os.path.realpath(__file__)); print(workdir)
workdir = workdir + '/'
datadir = workdir.removesuffix('_user/') + '_data/'; print(datadir)
outdir = workdir.removesuffix('_user/') + '_output/'; print(outdir)
print(workdir)
webdir = 'https://liveocean.apl.uw.edu/output/'
"""
"""
webdir = 'https://liveocean.apl.uw.edu/output/'

ti = datetime.strptime('2021.12.15', '%Y.%m.%d')

url_string = (webdir + 'f' + datetime.strftime(ti, '%Y.%m.%d') + '/layers.nc#mode=bytes')
ds1 = Dataset(url_string)

t=ds1['ocean_time'][:]
tmin=datetime.timestamp(datetime(2021,12,15))
tmax=datetime.timestamp(datetime(2021,12,16))
#print(tmin, tmax)

# PLOTTING
# plotting parameters
fs = 14 # primary fontsize
lw = 3 # primary linewidth
mk = 10 # primary markersize

plt.close('all')
fig = plt.figure(figsize=(6,10))
ax = fig.add_subplot(111)
ax.plot([1, 2, 3, 4, 5], [3, 7, 8, 4, 1])
test = datetime.strftime(ti, '%Y.%m.%d')
plt.title(test)
plt.show()
"""
