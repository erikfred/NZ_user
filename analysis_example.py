"""
analysis_example.py
Use web-hosted LiveOcean output to demonstrate the kind(s) of analysis I might
perform as part of my final dissertation chapter.
"""

# imports
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cmocean
from datetime import datetime, date

# read layers.nc file (will eventually become loop over all dates)
url1 = ('https://liveocean.apl.uw.edu/output/f2021.12.15/layers.nc#mode=bytes')
ds1 = Dataset(url1)

"""# print metadata
print(ds1.__dict__)
for dim in ds1.dimensions.values():
    print(dim)
for var in ds1.variables.values():
    print(var)
"""
layer_time = ds1['ocean_time'][:]
domain_lon = ds1['lon_rho'][:]
domain_lat = ds1['lat_rho'][:]

"""# -129.9798<lon<-122.018, 42.0067<lat<52.0099
print(domain_lon[0,0])
print(domain_lon[0,662])
print(domain_lat[0,0])
print(domain_lat[1301,0])
"""
minlon=-126
maxlon=-123.5
minlat=domain_lat[0,0]
maxlat=49

ilo1=np.argmin(np.abs(domain_lon[0,:]-minlon))
ilo2=np.argmin(np.abs(domain_lon[0,:]-maxlon))
ila1=np.argmin(np.abs(domain_lat[:,0]-minlat))
ila2=np.argmin(np.abs(domain_lat[:,0]-maxlat))

box_lon=domain_lon[ila1:ila2,ilo1:ilo2]
box_lat=domain_lat[ila1:ila2,ilo1:ilo2]

url2 = ('https://liveocean.apl.uw.edu/output/f2021.12.16/surface.nc#mode=bytes')
ds2 = Dataset(url2)

"""# print metadata
print(ds2.__dict__)
for dim in ds2.dimensions.values():
    print(dim)
for var in ds2.variables.values():
    print(var)
"""
surface_time=ds2['ocean_time'][:]
# ist=np.argwhere(surface_time==datetime.date(2021, 12, 15))
ssh=ds2['zeta'][:, ila1:ila2, ilo1:ilo2]
bath=ds2['h'][ila1:ila2, ilo1:ilo2]

tmin=datetime.timestamp(datetime(2021,12,15))
tmax=datetime.timestamp(datetime(2021,12,16))

it1=np.argwhere(surface_time>=tmin)
it2=np.argwhere(surface_time<tmax)

st=surface_time[it1:it2]
print(st)

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

"""
plt.close('all')
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.plot(test, ssh[:, 10, 10], 'r')
ax.xaxis.set_major_formatter(dates.DateFormatter('%m-%d'))
plt.show()

for x in range(0, 72):
    plt.close('all')
    print(datetime.fromtimestamp(surface_time[x]))
    # print(np.mean(ssh[x,:,:]))

    fig = plt.figure(figsize=(6,10))
    ax = fig.add_subplot(111)
    cs = ax.pcolormesh(XX,YY, ssh[x,:,:]-np.mean(ssh[x,:,:]), cmap=cmap, vmin=-0.2, vmax=0.2)
    bth = ax.contour(XX, YY, bath, [300, 2000], colors='black')

    ax.axis('square')
    ax.set_xlim((minlon,maxlon))
    ax.set_ylim((minlat, maxlat))
    ax.grid(True)

    cb = fig.colorbar(cs)
    cb.set_label('SSH (m)', fontsize=fs)
    cb.ax.tick_params(labelsize=fs)

    if not os.path.exists('../LO_output/SSH_maps'):
        os.mkdir('../LO_output/SSH_maps')

    plt.savefig('../LO_output/SSH_maps/temp' + str(x))
    # plt.show()
"""
