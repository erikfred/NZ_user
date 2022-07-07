"""
ssh_extractor.py
Generate SSH plots using history files on perigee. Goal is to do so efficiently
by only pulling the bare minimum data from the netCDF into working memory
"""

# imports
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cmocean
from datetime import datetime, date

# read history .nc files
# dir1 = "/data1/parker/LO_roms/cas6_v0_live/"
dir1 = "../LO_data/cas6_v0_live/"
dir2 = dir1 + "f2016.12.15/" # will eventually loop over these
dir3 = dir2 + "ocean_his_0001.nc" # with inner loop over these
ds1 = Dataset(dir3)

"""# print metadata
print(ds1.__dict__)
for dim in ds1.dimensions.values():
    print(dim)
for var in ds1.variables.values():
    print(var)
"""
domain_lon = ds1['lon_rho'][:]
domain_lat = ds1['lat_rho'][:]
"""
# full domain: -129.9798<lon<-122.018, 42.0067<lat<52.0099
print(domain_lon[0,0])
print(domain_lon[0,662])
print(domain_lat[0,0])
print(domain_lat[1301,0])
"""
minlon = -126
maxlon = -123.5
minlat = domain_lat[0,0]
maxlat = 49
if True:
    ilo1 = np.argmin(np.abs(domain_lon[0,:]-minlon))
    ilo2 = np.argmin(np.abs(domain_lon[0,:]-maxlon))
    ila1 = np.argmin(np.abs(domain_lat[:,0]-minlat))
    ila2 = np.argmin(np.abs(domain_lat[:,0]-maxlat))
else:
    ilo1 = 0
    ilo2 = 662
    ila1 = 0
    ila2 = 1301

box_lon=domain_lon[ila1:ila2,ilo1:ilo2]
box_lat=domain_lat[ila1:ila2,ilo1:ilo2]

t = ds1['ocean_time'][:]
ssh=ds1['zeta'][:, ila1:ila2, ilo1:ilo2]
bath=ds1['h'][ila1:ila2, ilo1:ilo2]

# PLOTTING
# plotting parameters
fs = 14 # primary fontsize
lw = 3 # primary linewidth
mk = 10 # primary markersize
cmap = cmocean.cm.thermal

XX = box_lon
YY = box_lat

# print(surface_time.shape)
# print(len(surface_time))


t = ds1['ocean_time'][:]
ssh=ds1['zeta'][:, ila1:ila2, ilo1:ilo2]
bath=ds1['h'][ila1:ila2, ilo1:ilo2]

plt.close('all')
print(datetime.fromtimestamp(t[0]))
# print(np.mean(ssh[x,:,:]))

fig = plt.figure(figsize=(6,10))
ax = fig.add_subplot(111)
cs = ax.pcolormesh(XX,YY, ssh[0,:,:]-np.mean(ssh[0,:,:]), cmap=cmap, vmin=-0.2, vmax=0.2)
bth = ax.contour(XX, YY, bath, [300, 2000], colors='black')

ax.axis('square')
ax.set_xlim((minlon,maxlon))
ax.set_ylim((minlat, maxlat))
ax.grid(True)

cb = fig.colorbar(cs)
cb.set_label('SSH (m)', fontsize=fs)
cb.ax.tick_params(labelsize=fs)

# plt.show()
if not os.path.exists('../LO_output/SSH_maps'):
    os.mkdir('../LO_output/SSH_maps')

plt.savefig('../LO_output/SSH_maps/temp' + str(x))
