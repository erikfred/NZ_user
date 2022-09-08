"""
ssh_extractor.py
Generate SSH plots using history files on perigee. Goal is to do so efficiently
by only pulling the bare minimum data from the netCDF into working memory
"""

# imports
import sys, os
from datetime import datetime, timedelta, date
import numpy as np
import netCDF4 as nc
import cmocean
import matplotlib.pyplot as plt
import pickle

sys.path.append(os.path.abspath('util'))
import Lfun
import zrfun
import zfun

# config
do_press = True # ?
testbatch = False # True will only process a few history files, stored locally
if testbatch:
    etag = '_test'
else:
    etag = ''
cutout = True # True will subsample from user-defined lat/lon limits

out_dir = '../LO_output/'
ostr = 'spinup/' # arbitrary label for separating runs
# in_dir = '/data1/parker/LO_roms/cas6_v0_live/' # location on perigee
in_dir = '../LO_data/cas6_v0_live/' # location locally
dstr = 'f' # naming convention for directories
fstr = 'ocean_his_' # naming convention for history files
ti = datetime.strptime('2016.12.15', '%Y.%m.%d')
tf = datetime.strptime('2016.12.19', '%Y.%m.%d')
# tf = datetime.strptime('2022.06.30', '%Y.%m.%d')
tag = 'LiveOcean'
n_layer = 0 # bottom layer

#### end setup ####

# get grid info first (things that won't change between files)
file0 = in_dir + dstr + datetime.strftime(ti, '%Y.%m.%d') + '/' + fstr + '0001.nc'
ds0 = nc.Dataset(file0) # opens netCDF file
if False: # set to True if you want to print netCDF info
    print('\nGRID INFO')
    zfun.ncd(ds0)
lon = ds0['lon_rho'][:,:]
lat = ds0['lat_rho'][:,:]
bath = ds0['h'][:,:]
if cutout: # take cutout from full model domain, if desired
    # full range is: -129.9798<lon<-122.018, 42.0067<lat<52.0099
    minlon = -126; maxlon = -123.5; minlat = 43; maxlat = 49
    ilo1 = np.argmin(np.abs(lon[0,:] - minlon)); ilo2 = np.argmin(np.abs(lon[0,:] - maxlon))
    ila1 = np.argmin(np.abs(lat[:,0] - minlat)); ila2 = np.argmin(np.abs(lat[:,0] - maxlat))
    lon = lon[ila1:ila2,ilo1:ilo2]
    lat = lat[ila1:ila2,ilo1:ilo2]
    bath = bath[ila1:ila2,ilo1:ilo2]
ds0.close()

# make file lists
file_list = []
if testbatch: # only runs on 2 day subset
    delta = (ti + timedelta(days=2)) - ti
else: # runs on entire date range defined in CONFIG
    delta = (tf + timedelta(days=1)) - ti
drange = range(0, delta.days)
for ii in drange:  # building file names
    d_ii = in_dir + dstr + datetime.strftime(ti + timedelta(days=ii), '%Y.%m.%d') + '/'
    for jj in range(1,25):
        f_jj = fstr + str(jj).zfill(4) + '.nc'
        file_list.append(d_ii + f_jj)
nf = len(file_list)

# prepare a directory for results
outdir0 = out_dir + ostr
Lfun.make_dir(outdir0, clean=False)
ncoutdir = outdir0 + 'SSH_extractions/'
Lfun.make_dir(ncoutdir, clean=False)
# output file
out_name = 'ssh_' + tag + etag + '.nc'
out_fn = ncoutdir + out_name
# get rid of the old version, if it exists
try:
    os.remove(out_fn)
except OSError:
    pass # assume error was because the file did not exist

# make SSH
ssh_arr = np.zeros((nf, np.shape(lon)[0], np.shape(lon)[1]))
t_arr = np.zeros(nf)
errors=[]
for tt in range(nf): #for each .nc file
    # if tt==60 or tt==64 or tt==71: # something is wrong with these files (at least on local)
    #     ssh_arr[tt,:,:] = ssh_arr[tt-1,:,:] # populate with previous day
    #     t_arr[tt] = t_arr[tt-1]
    #     continue
    filex = file_list[tt]
    try:
        dsx = nc.Dataset(filex)
    except:
        errors.append(filex)
        continue
    # except OSError: # skip day if file doesn't exist
    #     ssh_arr[tt,:,:] = ssh_arr[tt-1,:,:] # populate with previous day
    #     t_arr[tt] = t_arr[tt-1]
    #     print('no file ' + filex[len(filex)-28:len(filex)])
    #     continue
    # except MemoryError: # some files can't be read (unresolved)
    #     ssh_arr[tt,:,:] = ssh_arr[tt-1,:,:] # populate with previous day
    #     t_arr[tt] = t_arr[tt-1]
    #     print('unable to read ' + filex[len(filex)-28:len(filex)])
    #     continue
    if np.mod(tt,24)==0: # print update to console every 10th file
        print('tt = ' + str(tt) + '/' + str(nf) + ' ' + filex[len(filex)-28:len(filex)])
        # local_variables = [(var, sys.getsizeof(obj)) for var, obj in locals().items() if not var.startswith('__') and var not in ['sys',]]
        # current_size = sum([size for var, size in local_variables])
        # print('Total size of local variables:', current_size)
        sys.stdout.flush()
    tx = dsx['ocean_time'][:]
    if cutout:
        ssh = dsx['zeta'][:, ila1:ila2, ilo1:ilo2]
    else:
        ssh = dsx['zeta'][:,:,:]
    t_arr[tt] = tx
    ssh_arr[tt,:,:] = ssh
    dsx.close()
ssh_mean = np.mean(ssh_arr, axis=0) # mean pressure of given location for entire time range
ssh_anom = ssh_arr - ssh_mean

# SAVING
pickle.dump(t_arr, open((ncoutdir + 't_arr.p'), 'wb'))
pickle.dump(ssh_arr, open((ncoutdir + 'ssh_arr.p'), 'wb'))
pickle.dump(ssh_anom, open((ncoutdir + 'ssh_anom.p'), 'wb'))
pickle.dump(lat, open((ncoutdir + 'lat.p'), 'wb'))
pickle.dump(lon, open((ncoutdir + 'lon.p'), 'wb'))

"""
# PLOTTING
# plotting parameters
fs = 14 # primary fontsize
lw = 3 # primary linewidth
mk = 10 # primary markersize
cmap = cmocean.cm.thermal

plt.close('all')
print(datetime.fromtimestamp(t_arr[0]))

fig1 = plt.figure(figsize=(8,8))
ax1 = fig1.add_subplot(111)
pts = ax1.plot(t_arr,ssh_anom[:,30,100])

if not os.path.exists(ncoutdir + 'SSH_maps'):
    os.mkdir(ncoutdir + 'SSH_maps')

# plt.show()
plt.savefig(ncoutdir + 'SSH_maps/mooring')
plt.close()

# for sv in range(nf): # for each day
#     fig = plt.figure(figsize=(6,10))
#     ax = fig.add_subplot(111)
#     cs = ax.pcolormesh(lon, lat, ssh_anom[sv,:,:], cmap=cmap, vmin=-0.2, vmax=0.2)
#     bth = ax.contour(lon, lat, bath, [300, 2000], colors='black')
#
#     ax.axis('square')
#     # ax.set_xlim((minlon,maxlon))
#     # ax.set_ylim((minlat, maxlat))
#     ax.grid(True)
#
#     cb = fig.colorbar(cs)
#     cb.set_label('SSH (m)', fontsize=fs)
#     cb.ax.tick_params(labelsize=fs)
#
#     # plt.show()
#     if not os.path.exists(ncoutdir + 'SSH_maps'):
#         os.mkdir(ncoutdir + 'SSH_maps')
#     plt.savefig(ncoutdir + 'SSH_maps/temp' + str(sv).zfill(4))
#     plt.close()
"""
