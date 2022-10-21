import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pickle
import sys, os
from pathlib import Path
import time

from lo_tools import Lfun, zfun, zrfun

"""
if nn==36:
    rho2 = [ds['rho'][0,:,:,:]]
if nn==72:
    rho3 = [ds['rho'][0,:,:,:]]

rholp = np.append(rholp,[rho2],axis=0)
rho2 = rho3.copy()
"""

topdir = '../LO_output/allinone/'
loadir = topdir + 'pickles_2017-18/'
tlp = pickle.load(open((loadir + 'tlp.p'), 'rb'))
win = 90
print(range(len(tlp[:])-win))

# test = datetime.fromtimestamp(tlp[0])
# print(test.strftime("%m/%d") + " - " + test.strftime("%m/%d/%Y"))

# x1 = np.ones((365,1302,664))
# t = time.time()
# for nn in range(len(x1[0,0,0:100])):
#     x2 = np.correlate(x1,x1)
# print(time.time() - t)

# x1 = np.array([1, 2, 3, 4, 5])
# x2 = np.array([6, 3, 5, 2, 2])
# x3 = np.array([1, 2, 3, 4, 5])
# # x2 = np.array([[1, 2, 3, 4, 5],[6, 3, 5, 2, 2],[4, 4, 7, 9, 0]])
# a = (x1 - np.mean(x1)) / (np.std(x1) * len(x1))
# b = (x2 - np.mean(x2)) / (np.std(x2))
# c = (x3 - np.mean(x3)) / (np.std(x3))
# print(np.correlate(a, b))
# print(np.correlate(a,c))

# ncoutdir = '../LO_output/allinone/pickles/'
#
# tlp = pickle.load(open((ncoutdir + 'tlp.p'), 'rb'))
# zetalp = pickle.load(open((ncoutdir + 'zetalp.p'), 'rb'))
# rholp = pickle.load(open((ncoutdir + 'rholp.p'), 'rb'))
# lat = pickle.load(open((ncoutdir + 'lat.p'), 'rb'))
# lon = pickle.load(open((ncoutdir + 'lon.p'), 'rb'))
#
# print(np.shape(zetalp))
# print(np.shape(rholp))
#
# zetalp = np.append(zetalp,zetalp,axis=0)
# rholp = np.append(rholp,rholp,axis=0)
#
# print(np.shape(zetalp))
# print(np.shape(rholp))

# test = [[[1, 2, 3, 4]]]
# print(test)
# test2 = np.array(test)
# print(np.shape(test2))

# zetalp=[]
# zeta = np.ones((73, 500, 200),dtype=np.float32)
# zeta2 = zfun.lowpass(zeta, f='godin')
# print(zeta2)
# zetalp.append(zfun.lowpass(zeta,f='godin')[36:-36:24]); zeta = zeta[24:72,:,:]
# print(zeta3)
# zetalp.append(zeta3)
# print(zetalp)
# zeta = zeta[24:72,:,:]
# print(np.shape(zeta))

# test = np.around(91/5)
# test2 = range(int(test))
# print(test2)

# bath = pickle.load(open(('../LO_output/allinone/pickles/bath.p'), 'rb'))
# print(bath.shape)
# isb = np.argwhere((bath<302) & (bath>298))
# print(isb.shape)
# for nn in range(len(isb)):
#     print(isb[nn])

# t1 = pickle.load(open(('../LO_output/allinone/pickles/t_arr.p'), 'rb'))
# print(type(t1))
# print(len(t1))
# print(datetime.fromtimestamp(t1[0]))
#
# tlpd = np.zeros(len(t1))
# for tt in range(len(t1)):
#     tlpd[tt] = datetime.fromtimestamp(t1[tt])
#     print(tlpd[tt])
# print(len(tlpd))

# x = [[1, 2, 3, 4], [5, 6, 7, 8]]
# print(len(x).T)

# x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
# print(x[len(x):len(x)-5:-1])
# print(x[-1])

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
