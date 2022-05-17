from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

"""
url1 = ('https://liveocean.apl.uw.edu/output/f2021.12.15/layers.nc#mode=bytes')
ds1 = Dataset(url1)

print(ds1.__dict__)
for dim in ds1.dimensions.values():
    print(dim)
for var in ds1.variables.values():
    print(var)

model_time = ds1['ocean_time'][:]
domain_lon = ds1['lon_rho'][:]
domain_lat = ds1['lat_rho'][:]
# print(model_time)
print(domain_lon[0,0])
print(domain_lon[0,662])
print(domain_lat[0,0])
print(domain_lat[1301,0])
"""
url2 = ('https://liveocean.apl.uw.edu/output/f2021.12.15/surface.nc#mode=bytes')
ds2 = Dataset(url2)
"""
print(ds2.__dict__)
for dim in ds2.dimensions.values():
    print(dim)
for var in ds2.variables.values():
    print(var)
"""
ssh=ds2['zeta'][:]
print(np.mean(ssh[0,:,:]))
"""
plt.close('all')
fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(111)
ax.plot([1, 2, 3])
plt.show()
"""
