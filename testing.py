import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pickle
import sys, os
from pathlib import Path
import time
import cmocean

from lo_tools import Lfun, zfun, zrfun

"""
if nn==36:
    rho2 = [ds['rho'][0,:,:,:]]
if nn==72:
    rho3 = [ds['rho'][0,:,:,:]]

rholp = np.append(rholp,[rho2],axis=0)
rho2 = rho3.copy()
"""

# rholp = pickle.load(open(('../LO_output/allinone/pickles_2018-19/zetalp.p'), 'rb'))
# print(np.nanmin(rholp))
# print(str(np.nanmean(rholp)) + ' +/- ' + str(np.nanstd(rholp)))
# print(np.nanmax(rholp))

# make some things
g = 9.81
fn = '../LO_data/cas6_v0_live/f2016.12.17/ocean_his_0001.nc'
ds = nc.Dataset(fn)
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
h = ds['h'][:]
z = zrfun.get_z(h, 0*h, S, only_rho=True)
z0 = z[0,:,:].squeeze()
zetalp = np.asarray([ds['zeta'][0,:,:]])
# print(zetalp.shape)
rholp = np.asarray([ds['rho'][0,:,:,:]]).squeeze()
# print(rholp.shape)
ds.close()

Gh = G['h']
z_w = zrfun.get_z(Gh, zetalp.squeeze(), S, only_w=True)
DZ = np.diff(z_w, axis=0)
# print(np.sum(DZ[:,100,100]))
# print(rholp[0,:,100,100])
# print(np.sum(DZ[:,100,100]))
# print(g * 1025 * np.sum(DZ[:,100,100]))
bp_tot = (g * (1000 + rholp) * DZ).sum(axis=0)
# print(bp_tot[100,100])
# new method
ZW = zrfun.get_z(Gh, 0*zetalp.squeeze(), S, only_w=True)
DZ = np.diff(ZW, axis=0)
# print(np.sum(DZ[:,100,100]))
# calculate the baroclinic pressure
bp_bc = np.flip(np.cumsum(np.flip(g * (rholp - 0) * DZ, axis=0), axis=0), axis=0)[0,:,:]
bp_bc2 = (g * (rholp - 0) * DZ).sum(axis=0)
bp_ssh = (g * 1000 * (zetalp + DZ.sum(axis=0)))
# print(bp_ssh[0,100,100])
# print(bp_bc[100,100])
# print(bp_bc2[100,100])
print(bp_tot[100,100]-(bp_ssh[0,100,100]+bp_bc[100,100]))

"""
# total variability
print(np.nanmean(rholp))
print(np.nanstd(rholp))
# layer-wise variability
for nn in range(30):
    print(np.nanmean(rholp[0,nn,:,:]))
    print(np.nanstd(rholp[0,nn,:,:]))
    print(np.nanmean(DZ[nn,:,:]))
    print(' ')
"""


# print(range(0,3,step=0.5))
# print(list(range(0,3,step=0.5)))

# topdir = '../LO_output/allinone/'
# loadir = topdir + 'pickles_2018-19/'
# outdir = '../LO_output/mapview/'
# savedir = outdir + loadir[-8:]
# print(savedir)

# # PLOTTING
# # plotting parameters
# fs = 14 # primary fontsize
# lw = 3 # primary linewidth
# mk = 10 # primary markersize
# cmap = cmocean.cm.balance # formerly thermal
#
# lon = np.array([[1, 2, 3, 4, 5],[1, 2, 3, 4, 5],[1, 2, 3, 4, 5],[1, 2, 3, 4, 5], \
#     [1, 2, 3, 4, 5],[1, 2, 3, 4, 5],[1, 2, 3, 4, 5],[1, 2, 3, 4, 5],[1, 2, 3, 4, 5],[1, 2, 3, 4, 5]])
# lat = np.array([[1, 1, 1, 1, 1],[2, 2, 2, 2, 2],[3, 3, 3, 3, 3],[4, 4, 4, 4, 4], \
#     [5, 5, 5, 5, 5],[6, 6, 6, 6, 6],[7, 7, 7, 7, 7],[8, 8, 8, 8, 8],[9, 9, 9, 9, 9],[10, 10, 10, 10, 10]])
# c = np.random.random([10,5])
# print(c)
#
# fig1 = plt.figure(figsize=(8,8))
# ax1 = fig1.add_subplot(111)
# cs = ax1.pcolormesh(lon, lat, c, cmap=cmap, vmin=-1, vmax=1)
# plt.show()
# plt.close()
#
# c2 = np.random.random([10,5])
# print(c2)
#
# mask = np.absolute(c) < 0.2
# C2 = np.ma.masked_array(c2,mask=mask)
# fig1 = plt.figure(figsize=(8,8))
# ax1 = fig1.add_subplot(111)
# cs = ax1.pcolormesh(lon, lat, C2, cmap=cmap, vmin=-1, vmax=1)
# plt.show()
# plt.close()

# topdir = '../LO_output/allinone/'
# loadir = topdir + 'pickles_2017-18/'
# tlp = pickle.load(open((loadir + 'tlp.p'), 'rb'))
# win = 90
# print()

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
